#!/usr/bin/env python3
"""
Orchestrator for the undiagnosed rare disease diagnostic workflow.

Responsibilities:
  * own the case state and audit trail,
  * run deterministic lib stages and agentic (LLM) stages in order,
  * drive the structured adversarial debate (Stage 7) and enforce the survival gate,
  * assemble + validate report.json and render report.html.

The agentic stages are executed by an `AgentRunner` (LLM + Bash + API tools). A reference
implementation backed by the Anthropic SDK is provided; it is dependency-guarded so this file
imports and the deterministic parts run even without the SDK or a key. Each agentic stage gets
its prompt (prompts/<stage>.md + _shared_contract.md), the accumulated case context, and tool
access; it returns the JSON described in that prompt.

Usage:
  python3 orchestrator.py --case-id RD-... --snv-vcf ... --sv-vcf ... --bam ... \
      --hpo ... --notes ... --assembly GRCh38 --out results/
"""
from __future__ import annotations
import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib import vcf_filters, sv_filters, igv_runner, acmg  # noqa: E402

ROOT = Path(__file__).resolve().parent
PROMPTS = ROOT / "prompts"
SHARED_CONTRACT = (PROMPTS / "_shared_contract.md").read_text()


# --------------------------------------------------------------------------- AgentRunner
class AgentRunner:
    """
    Executes one agentic stage: feeds (system_prompt + shared_contract) and a JSON context,
    lets the model issue Bash commands and the API-client calls, and returns parsed JSON.

    Replace `complete()` with your LLM backend. The reference impl below uses the Anthropic
    SDK with tool use (Bash + python -c into lib/api_clients). Keep reviewer and proposer in
    SEPARATE runner instances so the red team does not inherit the proposer's chain of thought.
    """

    def __init__(self, model: str, audit_dir: Path, role: str = "proposer"):
        self.model = model
        self.audit_dir = audit_dir
        self.role = role

    def run_stage(self, stage_file: str, context: dict, expect: str = "json") -> dict | list:
        prompt = (PROMPTS / stage_file).read_text() + "\n\n" + SHARED_CONTRACT
        payload = {"system": prompt, "context": context}
        (self.audit_dir / f"{self.role}_{stage_file}.input.json").write_text(
            json.dumps(payload, indent=2, default=str))
        result = self.complete(prompt, context)
        (self.audit_dir / f"{self.role}_{stage_file}.output.json").write_text(
            json.dumps(result, indent=2, default=str))
        return result

    def complete(self, system_prompt: str, context: dict):
        """
        LLM call with tool use. Reference outline (pseudocode-complete):

            import anthropic
            client = anthropic.Anthropic()
            tools = [bash_tool, python_api_tool]   # Bash + `python3 -c` into lib/api_clients
            messages = [{"role": "user", "content": json.dumps(context)}]
            while True:
                resp = client.messages.create(model=self.model, system=system_prompt,
                                              tools=tools, messages=messages,
                                              max_tokens=8192, temperature=0)
                if resp.stop_reason == "tool_use":
                    messages.append({"role":"assistant","content":resp.content})
                    messages.append({"role":"user","content":run_tools(resp)})
                    continue
                return extract_json(resp)   # the stage's JSON object

        Until wired to a backend, raise so misconfiguration is loud rather than silent.
        """
        raise NotImplementedError(
            "Wire AgentRunner.complete() to your LLM backend (see docstring). "
            "Deterministic stages run without it; agentic stages require it.")


# --------------------------------------------------------------------------- helpers
def now() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def probe_tools() -> dict:
    cfg = yaml.safe_load((ROOT / "config" / "tools.yaml").read_text())
    versions = {}
    for name, spec in cfg["binaries"].items():
        try:
            out = subprocess.run(spec["version_cmd"], shell=True, capture_output=True,
                                 text=True, timeout=30)
            versions[name] = (out.stdout or out.stderr).splitlines()[0] if (out.stdout or out.stderr) else "present"
        except Exception as e:  # noqa: BLE001
            versions[name] = f"MISSING ({e})"
    return versions


def load_inputs(args) -> dict:
    hpo = []
    for line in Path(args.hpo).read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split(None, 1)
        hpo.append({"id": parts[0], "label": parts[1] if len(parts) > 1 else ""})
    return {"hpo_terms": hpo, "clinical_notes": Path(args.notes).read_text()}


# --------------------------------------------------------------------------- adversarial gate
def adversarial_debate(proposer: AgentRunner, reviewer: AgentRunner,
                       candidates: list[dict], rebuttal_rounds: int = 1) -> list[dict]:
    """
    Stage 7. For each candidate: reviewer attacks across all surfaces; proposer rebuts; reviewer
    re-grades. Candidates with an unrefuted FATAL challenge fail the survival gate.
    Returns candidates annotated with adversarial_review and survives_review.
    """
    out = []
    for cand in candidates:
        review = reviewer.run_stage("07_adversarial_review.md", {"candidate": cand})
        for _ in range(rebuttal_rounds):
            open_challenges = [c for c in review["adversarial_review"]["challenges"]
                               if c["status"] == "unrefuted"]
            if not open_challenges:
                break
            rebuttal = proposer.run_stage(
                "07_adversarial_review.md",
                {"candidate": cand, "mode": "rebut", "challenges": open_challenges})
            # reviewer adjudicates the rebuttal (re-grade)
            review = reviewer.run_stage(
                "07_adversarial_review.md",
                {"candidate": cand, "mode": "adjudicate",
                 "prior_review": review, "rebuttal": rebuttal})
        ar = review["adversarial_review"]
        fatal_open = any(c["severity"] == "fatal" and c["status"] == "unrefuted"
                         for c in ar["challenges"])
        ar["survives_review"] = not fatal_open
        cand["adversarial_review"] = ar
        if review.get("spawned_pending_evidence"):
            cand.setdefault("_pending_seeds", []).extend(review["spawned_pending_evidence"])
        out.append(cand)
    return out


# --------------------------------------------------------------------------- main pipeline
def run(args) -> None:
    out = Path(args.out)
    audit = out / "audit"
    igv_dir = out / "igv"
    for d in (out, audit, igv_dir):
        d.mkdir(parents=True, exist_ok=True)

    cfg = yaml.safe_load((ROOT / "config" / "tools.yaml").read_text())
    thresh = yaml.safe_load((ROOT / "config" / "thresholds.yaml").read_text())
    refs = cfg["reference_data"]["files"]
    proposer = AgentRunner(cfg["llm"]["proposer_model"], audit, role="proposer")
    reviewer = AgentRunner(cfg["llm"]["reviewer_model"], audit, role="reviewer")

    report: dict = {
        "schema_version": "1.0", "case_id": args.case_id, "generated_at": now(),
        "assembly": args.assembly,
        "pipeline": {"name": "rare-disease-workflow", "version": "1.0",
                     "stages_completed": [], "tool_versions": probe_tools(),
                     "api_sources": ["ClinVar", "PubMed", "Ensembl", "UCSC", "gnomAD"]},
    }
    inputs = load_inputs(args)

    # ---- Stage 0: intake & phenotype structuring (agentic) -------------------------------
    intake = proposer.run_stage("00_intake.md", inputs)
    report["patient"] = intake["patient"]
    search = intake["derived_search_strategy"]
    pending_seeds = list(intake.get("pending_evidence_seeds", []))
    report["pipeline"]["stages_completed"].append("intake")

    # ---- Stage 1: gene prioritisation (agentic) ------------------------------------------
    priors = proposer.run_stage("01_gene_prioritization.md",
                                {"patient": report["patient"], "search_strategy": search})
    tier_a = {g["symbol"] for g in priors["gene_priors"] if g["tier"] == "A"}
    report["pipeline"]["stages_completed"].append("gene_prioritization")

    # ---- Stage 2: SNV/indel triage (deterministic facts -> agent adjudication) -----------
    rows = vcf_filters.extract_rows(args.snv_vcf)
    n_in = len(rows)
    vcf_filters.apply_qc(rows)
    vcf_filters.annotate_repeat_context(rows, refs.get("segdup_bed"), refs.get("lcr_bed"))
    kept = vcf_filters.keep_by_consequence(rows, tier_a)
    comphet = vcf_filters.find_compound_hets(kept)
    snv_stage = proposer.run_stage("02_snv_triage.md", {
        "qc_table": vcf_filters.to_table(kept),
        "compound_het_groups": {g: [vcf_filters.asdict(r) for r in v]
                                for g, v in comphet.items()},
        "gene_priors": priors["gene_priors"], "search_strategy": search,
        "thresholds": thresh})
    report["pipeline"]["stages_completed"].append("snv_triage")

    # ---- Stage 3: SV/CNV triage --------------------------------------------------------
    sv_rows = sv_filters.extract_sv_rows(args.sv_vcf)
    sv_filters.apply_sv_qc(sv_rows)
    sv_filters.annotate_gene_overlap(sv_rows, refs.get("mane_summary"))  # gene bed expected
    sv_filters.annotate_population_overlap(sv_rows, refs.get("population_cnv_bed"))
    sv_filters.annotate_segdup(sv_rows, refs.get("segdup_bed"))
    sv_stage = proposer.run_stage("03_sv_triage.md", {
        "sv_table": sv_filters.to_table(sv_rows), "gene_priors": priors["gene_priors"],
        "snv_candidates": snv_stage["snv_candidates"], "thresholds": thresh})
    report["pipeline"]["stages_completed"].append("sv_triage")

    candidates = _assign_ids(snv_stage["snv_candidates"] + sv_stage["sv_candidates"])

    # ---- Stage 4: evidence gathering ---------------------------------------------------
    evidence = proposer.run_stage("04_evidence.md", {"candidates": candidates,
                                                     "patient": report["patient"]})
    candidates = _merge_by_id(candidates, evidence)
    report["pipeline"]["stages_completed"].append("evidence")

    # ---- Stage 5: IGV visual review ----------------------------------------------------
    for c in candidates:
        g = c.get("genomic", {})
        try:
            if c["variant_type"] in ("SNV", "indel", "MNV"):
                img = igv_runner.snapshot_snv(args.bam, g["chrom"], g["pos"],
                                              f"{c['candidate_id']}_{g['chrom']}_{g['pos']}",
                                              igv_dir, args.assembly)
                c.setdefault("quality", {}).setdefault("igv_assessment", {})["image"] = img
            else:
                imgs = igv_runner.snapshot_sv(args.bam, g["chrom"], g["start"], g["end"],
                                              g.get("sv_type", "SV"), c["candidate_id"],
                                              igv_dir, args.assembly)
                qa = c.setdefault("quality", {}).setdefault("igv_assessment", {})
                qa["image"], qa["extra_images"] = imgs[0], imgs[1:]
        except Exception as e:  # noqa: BLE001
            c.setdefault("quality", {}).setdefault("igv_assessment", {})["observations"] = \
                f"IGV capture failed: {e}"
    igv_review = proposer.run_stage("05_igv_review.md", {"candidates": candidates})
    candidates = _merge_by_id(candidates, igv_review)
    report["pipeline"]["stages_completed"].append("igv_review")

    # ---- Stage 6: ACMG (agent applies codes, lib computes class) -----------------------
    acmg_stage = proposer.run_stage("06_acmg.md", {"candidates": candidates})
    for item in acmg_stage:
        comp = acmg.classify(item["acmg"]["applied_criteria"])
        item["acmg"]["classification"] = comp["classification"]
        item["acmg"]["point_score"] = comp["point_score"]
        item["acmg"]["_warnings"] = comp["warnings"]
    candidates = _merge_by_id(candidates, acmg_stage)
    report["pipeline"]["stages_completed"].append("acmg")

    # ---- Stage 7: ADVERSARIAL review + survival gate -----------------------------------
    candidates = adversarial_debate(proposer, reviewer, candidates)
    survivors = [c for c in candidates if c["adversarial_review"]["survives_review"]]
    excluded = [{"locus": _locus(c), "gene": c.get("gene", {}).get("symbol", ""),
                 "stage_excluded": "adversarial_review",
                 "reason": c["adversarial_review"]["reviewer_verdict"]}
                for c in candidates if not c["adversarial_review"]["survives_review"]]
    for c in candidates:
        pending_seeds.extend(c.get("_pending_seeds", []))
    report["pipeline"]["stages_completed"].append("adversarial_review")

    # ---- Stage 8: synthesis & ranking --------------------------------------------------
    synthesis = proposer.run_stage("08_synthesis.md", {
        "survivors": survivors, "excluded": excluded, "patient": report["patient"],
        "counts": {"in": n_in, "candidates": len(candidates)}})
    report["candidate_variants"] = synthesis["candidate_variants"]
    report["excluded_variants_summary"] = synthesis.get("excluded_variants_summary", excluded)
    report["run_summary"] = synthesis.get("run_summary", {})

    # ---- Stage 9: pending evidence + plan review ---------------------------------------
    pending = proposer.run_stage("09_pending_evidence.md", {
        "candidates": report["candidate_variants"], "seeds": pending_seeds,
        "patient": report["patient"]})
    report["pending_evidence"] = pending
    # critically review every proposed experiment (separate reviewer context)
    for pe in report["pending_evidence"]:
        if pe.get("category") == "functional_experiment" and pe.get("experimental_plan"):
            pe["plan_review"] = reviewer.run_stage(
                "09_pending_evidence.md", {"mode": "review_plan", "pending_item": pe})
    report["pipeline"]["stages_completed"].append("pending_evidence")

    # ---- validate, write, render -------------------------------------------------------
    _validate(report)
    (out / "report.json").write_text(json.dumps(report, indent=2, default=str))
    from render.render_html import render
    render(report, out / "report.html")
    print(f"Done. {len(report['candidate_variants'])} candidates, "
          f"{len(report['pending_evidence'])} pending items.\n"
          f"  {out/'report.json'}\n  {out/'report.html'}")


# --------------------------------------------------------------------------- small utils
def _assign_ids(cands: list[dict]) -> list[dict]:
    for i, c in enumerate(cands, 1):
        c.setdefault("candidate_id", f"C{i:03d}")
    return cands


def _merge_by_id(base: list[dict], updates) -> list[dict]:
    upd = {u["candidate_id"]: u for u in (updates if isinstance(updates, list) else [updates])}
    for c in base:
        u = upd.get(c["candidate_id"])
        if u:
            for k, v in u.items():
                if k == "candidate_id":
                    continue
                c[k] = {**c.get(k, {}), **v} if isinstance(v, dict) and isinstance(c.get(k), dict) else v
    return base


def _locus(c: dict) -> str:
    g = c.get("genomic", {})
    if "pos" in g:
        return f"{g.get('chrom')}:{g.get('pos')} {g.get('ref')}>{g.get('alt')}"
    return f"{g.get('chrom')}:{g.get('start')}-{g.get('end')} {g.get('sv_type','')}"


def _validate(report: dict) -> None:
    try:
        import jsonschema
    except ImportError:
        print("note: jsonschema not installed; skipping schema validation", file=sys.stderr)
        return
    schema = json.loads((ROOT / "schema" / "candidate_report.schema.json").read_text())
    jsonschema.validate(report, schema)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--case-id", required=True)
    ap.add_argument("--snv-vcf", required=True)
    ap.add_argument("--sv-vcf", required=True)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--hpo", required=True)
    ap.add_argument("--notes", required=True)
    ap.add_argument("--assembly", default="GRCh38", choices=["GRCh38", "GRCh37"])
    ap.add_argument("--out", default="results/")
    run(ap.parse_args())


if __name__ == "__main__":
    main()

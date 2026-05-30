"""
Smoke test for the deterministic parts of the workflow. Exercises everything that does NOT need
the LLM backend, against the toy fixtures in examples/fixtures/.

Run:
    python3 examples/make_fixtures.py      # once, to create fixtures
    python3 tests/smoke_test.py

Covers: input validation, SV deterministic triage (bcftools/bedtools), ACMG engine,
schema validation of the example report, and HTML rendering. Exits non-zero on any failure.
"""
from __future__ import annotations
import json
import sys
from pathlib import Path
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
FIX = ROOT / "examples" / "fixtures"

PASS, FAIL = "  PASS", "  FAIL"
_failures = 0


def check(name: str, cond: bool, detail: str = "") -> None:
    global _failures
    print(f"{PASS if cond else FAIL}  {name}{(' — ' + detail) if detail else ''}")
    if not cond:
        _failures += 1


def test_fixtures_exist() -> None:
    needed = ["snv_indel.vcf.gz", "sv_cnv.vcf.gz", "proband.bam", "proband.bam.bai",
              "gene.bed.gz", "phenotype.hpo", "clinical_notes.txt"]
    missing = [f for f in needed if not (FIX / f).exists()]
    check("fixtures present", not missing,
          "missing: %s (run examples/make_fixtures.py)" % missing if missing else "")


def test_input_validation() -> None:
    from lib import validate_inputs
    args = SimpleNamespace(
        snv_vcf=str(FIX / "snv_indel.vcf.gz"), sv_vcf=str(FIX / "sv_cnv.vcf.gz"),
        bam=str(FIX / "proband.bam"), hpo=str(FIX / "phenotype.hpo"),
        notes=str(FIX / "clinical_notes.txt"), gene_bed=str(FIX / "gene.bed.gz"),
        assembly="GRCh38")
    res = validate_inputs.validate(args)
    validate_inputs.print_report(res, stream=sys.stdout)
    check("input validation passes on good fixtures", res.ok, "; ".join(res.errors))

    # negative control: missing/unindexed inputs must be reported, not crash
    bad = SimpleNamespace(snv_vcf="nope.vcf", sv_vcf=None, bam="nope.bam",
                          hpo="nope.hpo", notes="nope.txt", gene_bed=None, assembly="GRCh38")
    rbad = validate_inputs.validate(bad)
    check("input validation flags bad inputs", not rbad.ok and len(rbad.errors) >= 4,
          f"{len(rbad.errors)} errors")


def test_sv_triage_deterministic() -> None:
    from lib import sv_filters
    rows = sv_filters.extract_sv_rows(str(FIX / "sv_cnv.vcf.gz"))
    check("SV rows parsed", len(rows) == 2, f"{len(rows)} rows")
    sv_filters.apply_sv_qc(rows)
    sv_filters.annotate_gene_overlap(rows, str(FIX / "gene.bed.gz"))
    scn2a = [r for r in rows if "SCN2A" in r.genes_affected]
    check("SV gene overlap (DUP hits SCN2A)", len(scn2a) == 1,
          f"genes: {[r.genes_affected for r in rows]}")
    weak = [r for r in rows if r.sv_type == "DEL"][0]
    check("weak DEL flagged caution", weak.qc_verdict in ("caution", "fail"),
          f"verdict={weak.qc_verdict} flags={weak.flags}")


def test_acmg_engine() -> None:
    from lib import acmg
    p = acmg.classify([{"code": "PS1", "direction": "pathogenic", "strength": "strong"},
                       {"code": "PM1", "direction": "pathogenic", "strength": "moderate"},
                       {"code": "PM2", "direction": "pathogenic", "strength": "supporting"},
                       {"code": "PP3", "direction": "pathogenic", "strength": "moderate"},
                       {"code": "PM6", "direction": "pathogenic", "strength": "moderate"}])
    check("ACMG pathogenic sum", p["classification"] == "Pathogenic" and p["point_score"] == 11,
          str(p))
    b = acmg.classify([{"code": "BA1", "direction": "benign", "strength": "stand_alone"}])
    check("ACMG BA1 override", b["classification"] == "Benign")
    v = acmg.classify([{"code": "PVS1", "direction": "pathogenic", "strength": "moderate"},
                       {"code": "BP5", "direction": "benign", "strength": "supporting"}])
    check("ACMG VUS boundary", v["classification"] == "Uncertain significance", str(v))


def test_schema_and_render() -> None:
    report = json.loads((ROOT / "examples" / "example_report.json").read_text())
    try:
        import jsonschema
        schema = json.loads((ROOT / "schema" / "candidate_report.schema.json").read_text())
        jsonschema.validate(report, schema)
        check("example_report.json validates against schema", True)
    except ImportError:
        check("schema validation (jsonschema installed)", False, "pip install jsonschema")
    from render.render_html import render
    out = ROOT / "examples" / "_smoke_render.html"
    render(report, out)
    html = out.read_text()
    check("HTML render embeds data", "KCNQ2" in html and "/*__REPORT_JSON__*/null" not in html)
    out.unlink(missing_ok=True)


def main() -> None:
    print("\n===== RARE DISEASE WORKFLOW SMOKE TEST =====")
    test_fixtures_exist()
    test_input_validation()
    test_sv_triage_deterministic()
    test_acmg_engine()
    test_schema_and_render()
    print(f"\n{'='*44}\n{'ALL PASSED' if not _failures else str(_failures)+' FAILURE(S)'}\n")
    sys.exit(1 if _failures else 0)


if __name__ == "__main__":
    main()

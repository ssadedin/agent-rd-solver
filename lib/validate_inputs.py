"""
Validate workflow inputs before any analysis runs, and emit actionable, sourcing-aware errors.

Philosophy: fail fast, fail clearly. A clinician/bioinformatician should be able to read the
output and know exactly what is wrong and how to fix or obtain each input. Nothing here needs
the LLM backend.

Run standalone:
    python3 -m lib.validate_inputs --snv-vcf ... --sv-vcf ... --bam ... --hpo ... --notes ...
or import validate(args) from the orchestrator.
"""
from __future__ import annotations
import argparse
import gzip
import os
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path

HP_RE = re.compile(r"^HP:\d{7}\b")

# How to obtain / repair each input. Shown verbatim when a check fails.
SOURCING = {
    "snv_vcf": [
        "SNV/indel VCF must be VEP-annotated, bgzipped and tabix-indexed.",
        "  Annotate : vep -i in.vcf --vcf --everything --mane --hgvs --symbol \\",
        "                 --plugin REVEL --plugin SpliceAI --plugin CADD --plugin AlphaMissense \\",
        "                 --plugin LoF --af_gnomadg -o annotated.vcf",
        "  Compress : bgzip annotated.vcf",
        "  Index    : tabix -p vcf annotated.vcf.gz",
    ],
    "sv_vcf": [
        "SV/CNV VCF (e.g. Manta/GRIDSS/DELLY + a depth caller) must be bgzipped and indexed.",
        "  Merge callers (optional): bcftools merge / SURVIVOR, then:",
        "  bgzip sv.vcf && tabix -p vcf sv.vcf.gz",
        "  Records need SVTYPE and END (and SVLEN/PR/SR where the caller provides them).",
    ],
    "bam": [
        "Aligned reads (BAM/CRAM) must be coordinate-sorted and indexed:",
        "  samtools sort -o proband.bam in.bam && samtools index proband.bam",
        "  (CRAM also accepted; ensure the reference is available to samtools.)",
    ],
    "hpo": [
        "HPO terms: one per line, 'HP:0000001[ optional label]'. Obtain terms from:",
        "  - Phenotips / PhenoTagger from the clinical letter, or",
        "  - https://hpo.jax.org/  (browse/search), or the patient's referral.",
    ],
    "notes": [
        "Clinical notes: a UTF-8 text file with age of onset, family history, prior testing,",
        "consanguinity, ethnicity and exam findings. Free text is fine.",
    ],
    "gene_bed": [
        "Gene BED (chrom<TAB>start<TAB>end<TAB>symbol) for SV gene-overlap. Generate it from an",
        "authoritative source automatically:",
        "  python3 -m lib.fetch_gene_bed --assembly GRCh38 --mane-only \\",
        "      --out $WORKFLOW_DATA/tracks/gene.bed.gz",
        "  # optionally restrict to a PanelApp Australia panel:",
        "  python3 -m lib.fetch_gene_bed --assembly GRCh38 --panelapp 250 --out panel.bed.gz",
    ],
}


@dataclass
class Result:
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    info: list[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return not self.errors

    def err(self, key: str, msg: str) -> None:
        self.errors.append(msg)
        for line in SOURCING.get(key, []):
            self.errors.append("    " + line)


def _has(tool: str) -> bool:
    return shutil.which(tool) is not None


def _exists(p: str | None) -> bool:
    return bool(p) and Path(p).exists()


def _index_present(path: str, suffixes: tuple[str, ...]) -> bool:
    return any(Path(path + s).exists() for s in suffixes)


def _read_vcf_header(path: str) -> list[str]:
    """Read header lines of a (bgzipped) VCF without external tools."""
    opener = gzip.open if path.endswith(".gz") else open
    lines = []
    try:
        with opener(path, "rt", errors="replace") as fh:
            for line in fh:
                if line.startswith("#"):
                    lines.append(line.rstrip("\n"))
                else:
                    break
    except OSError as e:
        lines.append(f"__ERROR__ {e}")
    return lines


def _contigs_from_header(header: list[str]) -> list[str]:
    out = []
    for h in header:
        m = re.match(r"##contig=<.*ID=([^,>]+)", h)
        if m:
            out.append(m.group(1))
    return out


def check_vcf(path: str | None, key: str, res: Result, require_csq: bool = False) -> None:
    label = key.replace("_", " ").upper()
    if not _exists(path):
        res.err(key, f"[{label}] file not found: {path}")
        return
    if not str(path).endswith((".vcf.gz", ".bcf")):
        res.err(key, f"[{label}] must be bgzipped (.vcf.gz) — got {path}")
    if not _index_present(path, (".tbi", ".csi")):
        res.err(key, f"[{label}] missing tabix/csi index for {path}")
    header = _read_vcf_header(path)
    if header and header[0].startswith("__ERROR__"):
        res.err(key, f"[{label}] could not read (not valid bgzip?): {header[0]}")
        return
    if not any(h.startswith("#CHROM") for h in header):
        res.err(key, f"[{label}] no #CHROM line found — not a valid VCF header")
    if require_csq and not any("ID=CSQ" in h for h in header):
        res.warnings.append(
            f"[{label}] no VEP CSQ field in header — Stage 2 expects VEP annotation. "
            "Re-run VEP with --vcf (see sourcing below).")
        for line in SOURCING["snv_vcf"]:
            res.warnings.append("    " + line)
    contigs = _contigs_from_header(header)
    if contigs:
        res.info.append(f"[{label}] {len(contigs)} contigs, e.g. {contigs[:3]}")
    return


def check_bam(path: str | None, res: Result) -> list[str]:
    if not _exists(path):
        res.err("bam", f"[BAM] file not found: {path}")
        return []
    if not _index_present(path, (".bai", ".csi")) and not Path(str(path) + ".crai").exists():
        res.err("bam", f"[BAM] missing index ({path}.bai). Run: samtools index {path}")
    contigs: list[str] = []
    if _has("samtools"):
        try:
            out = subprocess.run(["samtools", "view", "-H", path],
                                 capture_output=True, text=True, timeout=60)
            contigs = [l.split("SN:")[1].split("\t")[0]
                       for l in out.stdout.splitlines() if l.startswith("@SQ") and "SN:" in l]
            res.info.append(f"[BAM] {len(contigs)} reference contigs")
            if not any(l.startswith("@HD") and "SO:coordinate" in l
                       for l in out.stdout.splitlines()):
                res.warnings.append("[BAM] header does not declare SO:coordinate — ensure sorted.")
        except Exception as e:  # noqa: BLE001
            res.warnings.append(f"[BAM] could not read header via samtools: {e}")
    else:
        res.warnings.append("[BAM] samtools not on PATH — skipped deep BAM checks.")
    return contigs


def check_hpo(path: str | None, res: Result) -> None:
    if not _exists(path):
        res.err("hpo", f"[HPO] file not found: {path}")
        return
    terms = [l for l in Path(path).read_text().splitlines()
             if l.strip() and not l.startswith("#")]
    valid = [t for t in terms if HP_RE.match(t.strip())]
    if not valid:
        res.err("hpo", f"[HPO] no valid 'HP:#######' lines found in {path}")
    else:
        res.info.append(f"[HPO] {len(valid)} valid term(s)")
    bad = [t for t in terms if not HP_RE.match(t.strip())]
    if bad:
        res.warnings.append(f"[HPO] {len(bad)} line(s) not recognised as HPO ids: {bad[:3]}")


def check_notes(path: str | None, res: Result) -> None:
    if not _exists(path):
        res.err("notes", f"[NOTES] file not found: {path}")
        return
    if Path(path).stat().st_size == 0:
        res.err("notes", f"[NOTES] file is empty: {path}")
    else:
        res.info.append(f"[NOTES] {Path(path).stat().st_size} bytes of clinical text")


def check_gene_bed(path: str | None, res: Result) -> None:
    if not path:
        res.warnings.append("[GENE BED] none supplied — SV gene-overlap annotation will be skipped.")
        for line in SOURCING["gene_bed"]:
            res.warnings.append("    " + line)
        return
    if not _exists(path):
        res.err("gene_bed", f"[GENE BED] file not found: {path}")
        return
    res.info.append(f"[GENE BED] present: {path}")


def check_assembly_concordance(vcf_contigs_key: dict[str, list[str]], assembly: str,
                               res: Result) -> None:
    """Warn on chr-prefix mismatch across inputs and against the declared assembly."""
    prefixed = {k: (any(c.startswith("chr") for c in v) if v else None)
                for k, v in vcf_contigs_key.items() if v}
    vals = set(prefixed.values())
    if len(vals) > 1:
        res.err("snv_vcf",
                f"[ASSEMBLY] contig naming differs across inputs (chr-prefix): {prefixed}. "
                "All inputs must share one convention; re-header with bcftools annotate "
                "--rename-chrs.")
    res.info.append(f"[ASSEMBLY] declared {assembly}; contig prefix style {prefixed}")


def validate(args) -> Result:
    res = Result()
    check_vcf(getattr(args, "snv_vcf", None), "snv_vcf", res, require_csq=True)
    check_vcf(getattr(args, "sv_vcf", None), "sv_vcf", res, require_csq=False)
    bam_contigs = check_bam(getattr(args, "bam", None), res)
    check_hpo(getattr(args, "hpo", None), res)
    check_notes(getattr(args, "notes", None), res)
    check_gene_bed(getattr(args, "gene_bed", None), res)

    contigs = {"bam": bam_contigs}
    for k in ("snv_vcf", "sv_vcf"):
        p = getattr(args, k, None)
        if _exists(p):
            contigs[k] = _contigs_from_header(_read_vcf_header(p))
    check_assembly_concordance(contigs, getattr(args, "assembly", "GRCh38"), res)
    return res


def print_report(res: Result, stream=sys.stderr) -> None:
    def w(s):
        print(s, file=stream)
    w("\n========== INPUT VALIDATION ==========")
    for i in res.info:
        w(f"  ok    {i}")
    for warn in res.warnings:
        w(f"  warn  {warn}" if not warn.startswith("    ") else warn)
    for e in res.errors:
        w(f"  ERROR {e}" if not e.startswith("    ") else e)
    w("--------------------------------------")
    w("RESULT: " + ("PASS — inputs look usable.\n" if res.ok
                    else f"FAIL — {len(res.errors)} blocking issue(s). Fix the items above.\n"))


def main() -> None:
    ap = argparse.ArgumentParser(description="Validate rare-disease workflow inputs.")
    for a in ("snv-vcf", "sv-vcf", "bam", "hpo", "notes", "gene-bed"):
        ap.add_argument("--" + a)
    ap.add_argument("--assembly", default="GRCh38")
    args = ap.parse_args()
    res = validate(args)
    print_report(res)
    sys.exit(0 if res.ok else 2)


if __name__ == "__main__":
    main()

"""
Deterministic SV/CNV triage. Builds the sv_table the Stage-3 agent reasons over:
type, size, breakpoint support, gene overlap, population-CNV overlap, dosage sensitivity.

Assumes a VCF with SV conventions (SVTYPE, END, SVLEN, and support fields such as PR/SR for
Manta or BC/CN for read-depth callers). Falls back gracefully when fields are absent.
Relies on bcftools + bedtools.
"""
from __future__ import annotations
import json
import shlex
import subprocess
from dataclasses import dataclass, field, asdict
from pathlib import Path

import yaml

THRESH = yaml.safe_load((Path(__file__).resolve().parent.parent /
                         "config" / "thresholds.yaml").read_text())["sv_cnv_qc"]


def _run(cmd: str) -> str:
    p = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"command failed: {cmd}\n{p.stderr}")
    return p.stdout


@dataclass
class SVRow:
    chrom: str
    start: int
    end: int
    sv_type: str
    svlen: int | None
    vcf_id: str
    filter: str
    copy_number: float | None = None
    paired_support: int | None = None
    split_support: int | None = None
    genes_affected: list[str] = field(default_factory=list)
    pop_overlap_fraction: float | None = None
    in_segdup: bool = False
    qc_verdict: str = "pass"
    flags: list[str] = field(default_factory=list)


def _header_tags(vcf: str) -> set[str]:
    """Tags actually declared in the header. bcftools query errors hard on absent tags, and
    SV vs read-depth CNV callers expose different INFO/FORMAT fields, so we query only what
    exists and treat the rest as null."""
    hdr = _run(f"bcftools view -h {shlex.quote(vcf)}")
    tags = set()
    for line in hdr.splitlines():
        m = None
        if line.startswith("##INFO=<ID="):
            m = "INFO/" + line.split("ID=", 1)[1].split(",", 1)[0]
        elif line.startswith("##FORMAT=<ID="):
            m = "FMT/" + line.split("ID=", 1)[1].split(",", 1)[0]
        if m:
            tags.add(m)
    return tags


def extract_sv_rows(vcf: str) -> list[SVRow]:
    tags = _header_tags(vcf)
    # column order is fixed; absent optional tags are emitted as "." so positions stay aligned
    cols = ["%CHROM", "%POS",
            "%INFO/END" if "INFO/END" in tags else ".",
            "%INFO/SVTYPE" if "INFO/SVTYPE" in tags else ".",
            "%INFO/SVLEN" if "INFO/SVLEN" in tags else ".",
            "%ID", "%FILTER",
            "%INFO/CN" if "INFO/CN" in tags else ".",
            "[%PR]" if "FMT/PR" in tags else ".",
            "[%SR]" if "FMT/SR" in tags else "."]
    fmt = "\\t".join(cols) + "\\n"
    out = _run(f"bcftools query -f '{fmt}' {shlex.quote(vcf)}")
    rows = []
    for line in out.splitlines():
        c = line.split("\t")
        if len(c) < 10:
            continue
        chrom, pos, end, svtype, svlen, vid, filt, cn_s, pr_s, sr_s = c[:10]
        start = int(pos)
        end_i = _int(end) or (start + abs(_int(svlen) or 0))
        rows.append(SVRow(chrom=chrom, start=start, end=end_i,
                          sv_type=(svtype if svtype != "." else "UNKNOWN"),
                          svlen=_int(svlen), vcf_id=vid, filter=filt, copy_number=_float(cn_s),
                          paired_support=_support(pr_s), split_support=_support(sr_s)))
    return rows


def _support(field_val: str):
    # PR/SR are usually ref,alt -> take the alt (second) count
    if not field_val or field_val == ".":
        return None
    parts = field_val.split(",")
    return _int(parts[-1])


def _int(x):
    try:
        return int(float(x))
    except (TypeError, ValueError):
        return None


def _float(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def apply_sv_qc(rows: list[SVRow]) -> list[SVRow]:
    for r in rows:
        size = abs(r.svlen) if r.svlen else (r.end - r.start)
        if size < THRESH["min_sv_size_bp"]:
            r.flags.append(f"too_small={size}")
            r.qc_verdict = "fail"
        if r.filter not in ("PASS", "."):
            r.flags.append(f"FILTER={r.filter}")
            r.qc_verdict = _worse(r.qc_verdict, "caution")
        # read-depth CNV callers won't have PR/SR; only enforce for breakpoint callers
        if r.paired_support is not None and r.paired_support < THRESH["min_paired_read_support"]:
            r.flags.append(f"low_PR={r.paired_support}")
            r.qc_verdict = _worse(r.qc_verdict, "caution")
        if r.split_support is not None and r.split_support < THRESH["min_split_read_support"]:
            r.flags.append(f"low_SR={r.split_support}")
            r.qc_verdict = _worse(r.qc_verdict, "caution")
    return rows


_ORDER = {"pass": 0, "caution": 1, "fail": 2}


def _worse(a, b):
    return a if _ORDER[a] >= _ORDER[b] else b


def annotate_gene_overlap(rows: list[SVRow], gene_bed: str) -> None:
    """gene_bed: chrom start end gene. Populates genes_affected via bedtools intersect."""
    if not gene_bed or not Path(gene_bed).exists():
        return
    loci = "\n".join(f"{r.chrom}\t{r.start}\t{r.end}\t{i}" for i, r in enumerate(rows))
    hit = _run(f"echo {shlex.quote(loci)} | "
               f"bedtools intersect -a - -b {shlex.quote(gene_bed)} -wa -wb")
    for line in hit.splitlines():
        f = line.split("\t")
        idx = int(f[3])
        gene = f[-1]
        if gene not in rows[idx].genes_affected:
            rows[idx].genes_affected.append(gene)


def annotate_population_overlap(rows: list[SVRow], pop_cnv_bed: str) -> None:
    """Reciprocal-overlap fraction against benign population CNV map."""
    if not pop_cnv_bed or not Path(pop_cnv_bed).exists():
        return
    loci = "\n".join(f"{r.chrom}\t{r.start}\t{r.end}\t{i}" for i, r in enumerate(rows))
    hit = _run(f"echo {shlex.quote(loci)} | bedtools intersect -a - -b {shlex.quote(pop_cnv_bed)} "
               f"-f 0.5 -r -wa -wb")
    seen = {int(l.split('\t')[3]) for l in hit.splitlines() if l}
    for i in seen:
        rows[i].pop_overlap_fraction = 0.5  # >= reciprocal 0.5; exact value computed if needed
        if rows[i].pop_overlap_fraction >= THRESH["max_benign_cnv_overlap_fraction"]:
            rows[i].flags.append("common_benign_cnv")
            rows[i].qc_verdict = _worse(rows[i].qc_verdict, "caution")


def annotate_segdup(rows: list[SVRow], segdup_bed: str) -> None:
    if not segdup_bed or not Path(segdup_bed).exists():
        return
    # flag if EITHER breakpoint sits in a segdup
    bps = []
    for i, r in enumerate(rows):
        bps.append(f"{r.chrom}\t{max(r.start-1,0)}\t{r.start}\t{i}")
        bps.append(f"{r.chrom}\t{max(r.end-1,0)}\t{r.end}\t{i}")
    hit = _run("echo " + shlex.quote("\n".join(bps)) +
               f" | bedtools intersect -a - -b {shlex.quote(segdup_bed)} -wa")
    for line in hit.splitlines():
        i = int(line.split("\t")[3])
        rows[i].in_segdup = True
        if "breakpoint_in_segdup" not in rows[i].flags:
            rows[i].flags.append("breakpoint_in_segdup")
            rows[i].qc_verdict = _worse(rows[i].qc_verdict, "caution")


def to_table(rows: list[SVRow]) -> list[dict]:
    return [asdict(r) for r in rows]

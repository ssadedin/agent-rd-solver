"""
Deterministic SNV/indel triage on a VEP-annotated VCF.

This produces the *facts* the Stage-2 agent reasons over: a QC table, frequency table,
consequence table (VEP CSQ parsed, MANE/canonical preferred) and genotype table. It does NOT
make the final keep/drop call on borderline variants — it labels them so the agent (and the IGV
stage) can adjudicate. Every hard cut comes from config/thresholds.yaml.

Relies on `bcftools` (with the +split-vep plugin) being on PATH, per config/tools.yaml.
"""
from __future__ import annotations
import json
import shlex
import subprocess
from dataclasses import dataclass, field, asdict
from pathlib import Path

import yaml

THRESH = yaml.safe_load((Path(__file__).resolve().parent.parent /
                         "config" / "thresholds.yaml").read_text())


def _run(cmd: str) -> str:
    """Run a Bash command in the sandbox and return stdout (raises on failure)."""
    p = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"command failed: {cmd}\n{p.stderr}")
    return p.stdout


def csq_fields(vcf: str) -> list[str]:
    """Parse the VEP CSQ field order from the VCF header."""
    hdr = _run(f"bcftools view -h {shlex.quote(vcf)} | grep -m1 'ID=CSQ'")
    # ...Format: Allele|Consequence|IMPACT|SYMBOL|...">
    fmt = hdr.split("Format:")[1].strip().rstrip('">').strip()
    return [f.strip() for f in fmt.split("|")]


@dataclass
class VariantRow:
    chrom: str
    pos: int
    ref: str
    alt: str
    vcf_id: str = "."
    filter: str = "."
    # genotype / quality
    gt: str = "./."
    zygosity: str = "unknown"
    depth: int | None = None
    gq: int | None = None
    ad_ref: int | None = None
    ad_alt: int | None = None
    allele_balance: float | None = None
    mq: float | None = None
    strand_bias: float | None = None
    # vep
    gene: str = ""
    transcript: str = ""
    consequence: str = ""
    impact: str = ""
    hgvsc: str = ""
    hgvsp: str = ""
    mane_select: bool = False
    # frequency / predictors (from VEP if present; else filled by api_clients)
    gnomad_af: float | None = None
    revel: float | None = None
    cadd: float | None = None
    spliceai: float | None = None
    # derived flags
    flags: list[str] = field(default_factory=list)
    qc_verdict: str = "pass"


def _zygosity(gt: str) -> str:
    a = gt.replace("|", "/").split("/")
    if len(a) != 2 or "." in a:
        return "unknown"
    if a[0] == a[1] == "0":
        return "ref"
    if a[0] == a[1]:
        return "hom"
    return "het"


def extract_rows(vcf: str, regions_bed: str | None = None) -> list[VariantRow]:
    """
    Stream the VCF through bcftools +split-vep, pulling per-variant + per-(MANE) transcript
    annotations and the proband genotype/format fields.
    """
    fields = csq_fields(vcf)
    want_csq = [f for f in ("SYMBOL", "Feature", "Consequence", "IMPACT", "HGVSc", "HGVSp",
                            "MANE_SELECT", "gnomADg_AF", "REVEL", "CADD_PHRED",
                            "SpliceAI_pred_DS_AG") if f in fields]
    csq_expr = "".join(f"%{c}|" for c in want_csq).rstrip("|")
    region = f"-R {shlex.quote(regions_bed)}" if regions_bed else ""
    # -s worst: keep the most severe consequence; prefer MANE via -p / select transcripts.
    query = (
        f"bcftools +split-vep {shlex.quote(vcf)} {region} "
        f"-f '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%FILTER\\t"
        f"[%GT]\\t[%DP]\\t[%GQ]\\t[%AD]\\t%INFO/MQ\\t%INFO/FS\\t{csq_expr}\\n' "
        f"-d -s worst"
    )
    out = _run(query)
    rows: list[VariantRow] = []
    for line in out.splitlines():
        parts = line.split("\t")
        if len(parts) < 12:
            continue
        chrom, pos, vid, ref, alt, filt, gt, dp, gq, ad, mq, fs = parts[:12]
        csq = parts[12].split("|") if len(parts) > 12 else []
        cmap = dict(zip(want_csq, csq))
        ad_ref = ad_alt = None
        if ad and ad != ".":
            ad_parts = ad.split(",")
            if len(ad_parts) >= 2 and ad_parts[0].isdigit():
                ad_ref, ad_alt = int(ad_parts[0]), int(ad_parts[1])
        ab = (ad_alt / (ad_ref + ad_alt)) if (ad_ref is not None and (ad_ref + ad_alt) > 0) else None
        rows.append(VariantRow(
            chrom=chrom, pos=int(pos), ref=ref, alt=alt, vcf_id=vid, filter=filt,
            gt=gt, zygosity=_zygosity(gt),
            depth=_int(dp), gq=_int(gq), ad_ref=ad_ref, ad_alt=ad_alt, allele_balance=ab,
            mq=_float(mq), strand_bias=_float(fs),
            gene=cmap.get("SYMBOL", ""), transcript=cmap.get("Feature", ""),
            consequence=cmap.get("Consequence", ""), impact=cmap.get("IMPACT", ""),
            hgvsc=cmap.get("HGVSc", ""), hgvsp=cmap.get("HGVSp", ""),
            mane_select=bool(cmap.get("MANE_SELECT", "").strip()),
            gnomad_af=_float(cmap.get("gnomADg_AF")), revel=_float(cmap.get("REVEL")),
            cadd=_float(cmap.get("CADD_PHRED")), spliceai=_float(cmap.get("SpliceAI_pred_DS_AG")),
        ))
    return rows


def _int(x):
    try:
        return int(x)
    except (TypeError, ValueError):
        return None


def _float(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def apply_qc(rows: list[VariantRow]) -> list[VariantRow]:
    """Annotate each row with QC flags and a verdict (pass/caution/fail). Never deletes rows."""
    q = THRESH["snv_indel_qc"]
    for r in rows:
        if q["require_filter_pass"] and r.filter not in ("PASS", "."):
            r.flags.append(f"FILTER={r.filter}")
            r.qc_verdict = "fail"
        if r.depth is not None and r.depth < q["min_depth"]:
            r.flags.append(f"low_depth={r.depth}")
            r.qc_verdict = _worse(r.qc_verdict, "fail")
        if r.gq is not None and r.gq < q["min_genotype_quality"]:
            r.flags.append(f"low_gq={r.gq}")
            r.qc_verdict = _worse(r.qc_verdict, "caution")
        if r.ad_alt is not None and r.ad_alt < q["min_alt_reads"]:
            r.flags.append(f"low_alt_reads={r.ad_alt}")
            r.qc_verdict = _worse(r.qc_verdict, "caution")
        if r.zygosity == "het" and r.allele_balance is not None:
            if not (q["het_allele_balance_min"] <= r.allele_balance <= q["het_allele_balance_max"]):
                r.flags.append(f"skewed_ab={r.allele_balance:.2f}")
                r.qc_verdict = _worse(r.qc_verdict, "caution")
        if r.mq is not None and r.mq < q["min_mapping_quality"]:
            r.flags.append(f"low_mq={r.mq}")
            r.qc_verdict = _worse(r.qc_verdict, "caution")
        if r.strand_bias is not None and r.strand_bias > q["max_strand_bias_phred"]:
            r.flags.append(f"strand_bias={r.strand_bias}")
            r.qc_verdict = _worse(r.qc_verdict, "caution")
    return rows


_ORDER = {"pass": 0, "caution": 1, "fail": 2}


def _worse(a: str, b: str) -> str:
    return a if _ORDER[a] >= _ORDER[b] else b


def keep_by_consequence(rows: list[VariantRow], tier_a_genes: set[str]) -> list[VariantRow]:
    """Keep HIGH/MODERATE impact; keep conditional consequences only in Tier-A genes or with
    a splice/predictor signal. Returns the retained subset (the rest go to 'excluded')."""
    cp = THRESH["consequence_priority"]
    keep_terms = set(cp["high_impact"]) | set(cp["moderate_impact"])
    cond = set(cp["conditional"])
    sp_thr = THRESH["in_silico"]["spliceai_supporting"]
    kept = []
    for r in rows:
        terms = set(r.consequence.split("&"))
        if terms & keep_terms:
            kept.append(r)
        elif terms & cond and (
            r.gene in tier_a_genes
            or (r.spliceai is not None and r.spliceai >= sp_thr)
        ):
            r.flags.append("conditional_consequence_retained")
            kept.append(r)
    return kept


def annotate_repeat_context(rows: list[VariantRow], segdup_bed: str, lcr_bed: str) -> None:
    """Flag variants in segdup / low-complexity regions (mandatory IGV)."""
    for bed, label in ((segdup_bed, "segdup"), (lcr_bed, "low_complexity")):
        if not bed or not Path(bed).exists():
            continue
        # write loci to a temp bed and intersect
        loci = "\n".join(f"{r.chrom}\t{r.pos-1}\t{r.pos}\t{i}" for i, r in enumerate(rows))
        hit = _run(f"echo {shlex.quote(loci)} | bedtools intersect -a - -b {shlex.quote(bed)} -wa")
        idx = {int(l.split('\t')[3]) for l in hit.splitlines() if l}
        for i in idx:
            rows[i].flags.append(f"in_{label}")
            rows[i].qc_verdict = _worse(rows[i].qc_verdict, "caution")


def find_compound_hets(rows: list[VariantRow]) -> dict[str, list[VariantRow]]:
    """Group rare het variants by gene as candidate compound-het pairs (phasing done later via BAM)."""
    by_gene: dict[str, list[VariantRow]] = {}
    rec_max = THRESH["frequency"]["recessive_max_af"]
    for r in rows:
        if r.zygosity == "het" and (r.gnomad_af is None or r.gnomad_af <= rec_max):
            by_gene.setdefault(r.gene, []).append(r)
    return {g: v for g, v in by_gene.items() if len(v) >= 2}


def to_table(rows: list[VariantRow]) -> list[dict]:
    return [asdict(r) for r in rows]


if __name__ == "__main__":
    import sys
    rows = extract_rows(sys.argv[1])
    apply_qc(rows)
    print(json.dumps(to_table(rows)[:5], indent=2))

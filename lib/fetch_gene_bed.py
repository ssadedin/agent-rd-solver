"""
Fetch an authoritative gene BED (chrom, start, end, gene_symbol) used by the SV stage for
gene-overlap annotation, with no manual curation.

Sources (authoritative, automatable):
  --source gencode  (default)  GENCODE basic annotation GTF from EBI. Versioned, stable URLs,
                               the reference gene model for GRCh38/GRCh37. Gene-level features.
  --source ucsc                UCSC 'mane'/'ncbiRefSeqSelect' track via the UCSC REST API.

Optional clinical restriction:
  --panelapp <panel_id>        Restrict to genes on a PanelApp Australia panel
                               (https://panelapp-aus.org). Coordinates still come from
                               GENCODE/UCSC; PanelApp supplies the gene symbol whitelist and
                               confidence (green = level 3).

  --mane-only                  Keep only MANE Select transcripts (GENCODE 'tag "MANE_Select"').

Output is BED4, coordinate-sorted, and bgzip+tabix-indexed when those tools are available.

Examples:
  python3 -m lib.fetch_gene_bed --assembly GRCh38 --mane-only --out tracks/gene.bed.gz
  python3 -m lib.fetch_gene_bed --assembly GRCh38 --panelapp 250 --out panel.bed.gz
"""
from __future__ import annotations
import argparse
import gzip
import io
import re
import shutil
import subprocess
import sys
from pathlib import Path

import requests

# GENCODE release pinned for reproducibility; bump deliberately.
GENCODE_RELEASE = "46"
GENCODE_URL = {
    "GRCh38": f"https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/"
              f"release_{GENCODE_RELEASE}/gencode.v{GENCODE_RELEASE}.basic.annotation.gtf.gz",
    "GRCh37": f"https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/"
              f"release_{GENCODE_RELEASE}/GRCh37_mapping/"
              f"gencode.v{GENCODE_RELEASE}lift37.basic.annotation.gtf.gz",
}
UCSC_API = "https://api.genome.ucsc.edu"
UCSC_GENOME = {"GRCh38": "hg38", "GRCh37": "hg19"}
PANELAPP_AU = "https://panelapp-aus.org/api/v1"

_ATTR_GENE = re.compile(r'gene_name "([^"]+)"')


# --------------------------------------------------------------------------- GENCODE
def from_gencode(assembly: str, mane_only: bool, max_records: int | None = None) -> list[tuple]:
    url = GENCODE_URL[assembly]
    print(f"[gencode] streaming {url}", file=sys.stderr)
    feature = "transcript" if mane_only else "gene"
    rows: list[tuple] = []
    with requests.get(url, stream=True, timeout=120) as r:
        r.raise_for_status()
        gz = gzip.GzipFile(fileobj=io.BufferedReader(r.raw))  # type: ignore[arg-type]
        for raw in gz:
            line = raw.decode("utf-8", "replace")
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != feature:
                continue
            attrs = f[8]
            if mane_only and 'tag "MANE_Select"' not in attrs:
                continue
            m = _ATTR_GENE.search(attrs)
            if not m:
                continue
            rows.append((f[0], int(f[3]) - 1, int(f[4]), m.group(1)))
            if max_records and len(rows) >= max_records:
                break
    return _collapse_genes(rows)


def _collapse_genes(rows: list[tuple]) -> list[tuple]:
    """For mane-only/transcript input, collapse multiple transcripts to one span per gene."""
    span: dict[tuple, list] = {}
    for chrom, start, end, name in rows:
        key = (chrom, name)
        if key not in span:
            span[key] = [start, end]
        else:
            span[key][0] = min(span[key][0], start)
            span[key][1] = max(span[key][1], end)
    return [(c, s, e, n) for (c, n), (s, e) in span.items()]


# --------------------------------------------------------------------------- UCSC
def from_ucsc(assembly: str, max_records: int | None = None) -> list[tuple]:
    genome = UCSC_GENOME[assembly]
    track = "mane"  # MANE Select on hg38; falls back to ncbiRefSeqSelect on hg19
    print(f"[ucsc] genome={genome} track={track}", file=sys.stderr)
    chroms = _ucsc_chroms(genome)
    rows: list[tuple] = []
    for chrom in chroms:
        data = requests.get(f"{UCSC_API}/getData/track",
                            params={"genome": genome, "track": track, "chrom": chrom},
                            timeout=120).json()
        items = data.get(track) or data.get(chrom) or []
        if isinstance(items, dict):
            items = items.get(chrom, [])
        for it in items:
            name = it.get("geneName") or it.get("name2") or it.get("name")
            if name and "txStart" in it:
                rows.append((chrom, int(it["txStart"]), int(it["txEnd"]), name))
        if max_records and len(rows) >= max_records:
            break
    return _collapse_genes(rows)


def _ucsc_chroms(genome: str) -> list[str]:
    data = requests.get(f"{UCSC_API}/list/chromosomes",
                        params={"genome": genome}, timeout=60).json()
    chroms = list((data.get("chromosomes") or {}).keys())
    # main chromosomes only (skip alts/randoms/Un)
    return [c for c in chroms if re.fullmatch(r"chr([0-9]{1,2}|X|Y|M|MT)", c)]


# --------------------------------------------------------------------------- PanelApp AU
def panelapp_genes(panel_id: str, min_confidence: int = 3) -> set[str]:
    """Return green/confident gene symbols for a PanelApp Australia panel."""
    url = f"{PANELAPP_AU}/panels/{panel_id}/"
    print(f"[panelapp] {url}", file=sys.stderr)
    try:
        r = requests.get(url, timeout=60)
        r.raise_for_status()
    except requests.RequestException as e:
        raise SystemExit(
            f"ERROR: could not reach PanelApp Australia ({e}).\n"
            f"  Check connectivity to {PANELAPP_AU} or browse panels at "
            f"https://panelapp-aus.org/panels/ to find a panel id.\n"
            f"  You can also export a panel's gene list manually and pass --gene-list file.txt")
    genes = r.json().get("genes", [])
    return {g["gene_data"]["gene_symbol"] for g in genes
            if int(g.get("confidence_level", 0)) >= min_confidence and g.get("gene_data")}


# --------------------------------------------------------------------------- write
def write_bed(rows: list[tuple], out: str, restrict: set[str] | None = None) -> Path:
    if restrict is not None:
        rows = [r for r in rows if r[3] in restrict]
        missing = restrict - {r[3] for r in rows}
        if missing:
            print(f"[warn] {len(missing)} panel gene(s) had no coordinates "
                  f"(symbol/alias mismatch): {sorted(missing)[:10]}", file=sys.stderr)
    rows = sorted(set(rows), key=lambda r: (_chrom_key(r[0]), r[1]))
    out_path = Path(out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plain = out_path.with_suffix("") if out_path.suffix == ".gz" else out_path
    with open(plain, "w") as fh:
        for chrom, start, end, name in rows:
            fh.write(f"{chrom}\t{start}\t{end}\t{name}\n")
    print(f"[ok] wrote {len(rows)} gene intervals -> {plain}", file=sys.stderr)
    if out.endswith(".gz") and shutil.which("bgzip") and shutil.which("tabix"):
        subprocess.run(["bgzip", "-f", str(plain)], check=True)
        subprocess.run(["tabix", "-p", "bed", str(out_path)], check=True)
        print(f"[ok] bgzipped + tabixed -> {out_path}", file=sys.stderr)
        return out_path
    if out.endswith(".gz"):
        print("[warn] bgzip/tabix not found; left uncompressed BED.", file=sys.stderr)
    return plain


def _chrom_key(c: str):
    c2 = c.replace("chr", "")
    return (0, int(c2)) if c2.isdigit() else (1, c2)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--assembly", default="GRCh38", choices=["GRCh38", "GRCh37"])
    ap.add_argument("--source", default="gencode", choices=["gencode", "ucsc"])
    ap.add_argument("--mane-only", action="store_true",
                    help="GENCODE: keep only MANE Select transcripts")
    ap.add_argument("--panelapp", help="PanelApp Australia panel id to restrict genes to")
    ap.add_argument("--panelapp-min-confidence", type=int, default=3,
                    help="3=green (default), 2=amber, 1=red")
    ap.add_argument("--gene-list", help="local file of gene symbols (one per line) to restrict to")
    ap.add_argument("--max-records", type=int, help="dev/testing: stop after N source records")
    ap.add_argument("--out", required=True, help="output .bed or .bed.gz")
    args = ap.parse_args()

    restrict: set[str] | None = None
    if args.panelapp:
        restrict = panelapp_genes(args.panelapp, args.panelapp_min_confidence)
        print(f"[panelapp] {len(restrict)} genes at confidence>={args.panelapp_min_confidence}",
              file=sys.stderr)
    if args.gene_list:
        gl = {l.strip() for l in Path(args.gene_list).read_text().splitlines() if l.strip()}
        restrict = gl if restrict is None else (restrict | gl)

    if args.source == "gencode":
        rows = from_gencode(args.assembly, args.mane_only, args.max_records)
    else:
        rows = from_ucsc(args.assembly, args.max_records)
    write_bed(rows, args.out, restrict)


if __name__ == "__main__":
    main()

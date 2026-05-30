"""
Generate IGV screenshots for candidate loci by writing an IGV batch script and running it
headless (xvfb-run igv -b batch.txt). Snapshots are saved under <out>/igv/ and the relative
path is returned for embedding in the report and the evidence_ref.

For SNV/indel: a single padded window.
For SV/CNV: a snapshot at each breakpoint plus a zoomed-out view (and a coverage view).
"""
from __future__ import annotations
import shlex
import subprocess
from pathlib import Path

import yaml

_CFG = yaml.safe_load((Path(__file__).resolve().parent.parent /
                       "config" / "tools.yaml").read_text())["igv"]


def _genome_id(assembly: str) -> str:
    return _CFG["genome_id_grch38"] if assembly == "GRCh38" else _CFG["genome_id_grch37"]


def _batch_header(bam: str, genome: str, out_dir: Path) -> list[str]:
    lines = ["new", f"genome {genome}", f"load {bam}", f"snapshotDirectory {out_dir}"]
    for k, v in _CFG.get("preferences", {}).items():
        lines.append(f"preference {k} {v}")
    return lines


def snapshot_snv(bam: str, chrom: str, pos: int, name: str, out_dir: Path,
                 assembly: str = "GRCh38") -> str:
    flank = _CFG["flank_bp"]
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    png = f"{name}.png"
    batch = _batch_header(bam, _genome_id(assembly), out_dir)
    batch += [f"goto {chrom}:{max(pos-flank,1)}-{pos+flank}",
              "sort base", "collapse", f"snapshot {png}", "exit"]
    _execute(batch, out_dir)
    return f"igv/{png}"


def snapshot_sv(bam: str, chrom: str, start: int, end: int, sv_type: str, name: str,
                out_dir: Path, assembly: str = "GRCh38") -> list[str]:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    genome = _genome_id(assembly)
    flank = _CFG["flank_bp"]
    images: list[str] = []
    batch = _batch_header(bam, genome, out_dir)
    span = end - start
    if sv_type in ("BND", "INV") or span > _CFG["sv_view_max_bp"]:
        # snapshot each breakpoint separately
        for label, p in (("bp1", start), ("bp2", end)):
            png = f"{name}_{label}.png"
            batch += [f"goto {chrom}:{max(p-flank*4,1)}-{p+flank*4}",
                      "viewaspairs", "squish", f"snapshot {png}"]
            images.append(f"igv/{png}")
    else:
        png = f"{name}.png"
        batch += [f"goto {chrom}:{max(start-flank*4,1)}-{end+flank*4}",
                  "viewaspairs", "squish", f"snapshot {png}"]
        images.append(f"igv/{png}")
    # zoomed-out coverage context
    png_ctx = f"{name}_context.png"
    pad = max(span, 2000)
    batch += [f"goto {chrom}:{max(start-pad,1)}-{end+pad}", "collapse",
              f"snapshot {png_ctx}", "exit"]
    images.append(f"igv/{png_ctx}")
    _execute(batch, out_dir)
    return images


def _execute(batch_lines: list[str], out_dir: Path) -> None:
    batch_file = out_dir / "_igv_batch.txt"
    batch_file.write_text("\n".join(batch_lines) + "\n")
    wrapper = _wrapper()
    cmd = f"{wrapper} igv -b {shlex.quote(str(batch_file))}"
    p = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=600)
    if p.returncode != 0:
        raise RuntimeError(f"IGV batch failed: {p.stderr[-2000:]}")


def _wrapper() -> str:
    try:
        from . import tools_cfg  # not present; intentional fallthrough
    except Exception:
        pass
    full = yaml.safe_load((Path(__file__).resolve().parent.parent /
                           "config" / "tools.yaml").read_text())
    return full["binaries"]["igv"].get("headless_wrapper", "xvfb-run -a")

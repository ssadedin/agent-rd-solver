"""
Render report.json into a single self-contained report.html.

The JSON is embedded into template.html (no external assets, no server needed) and rendered
client-side. IGV PNGs are referenced by their relative paths under igv/.
"""
from __future__ import annotations
import json
from pathlib import Path

TEMPLATE = Path(__file__).resolve().parent / "template.html"


def render(report: dict, out_path: Path | str) -> Path:
    html = TEMPLATE.read_text()
    blob = json.dumps(report, default=str)
    # Embed safely: close-tag escape so a string in the JSON can't break out of <script>.
    blob = blob.replace("</", "<\\/")
    html = html.replace("/*__REPORT_JSON__*/null", blob)
    out_path = Path(out_path)
    out_path.write_text(html)
    return out_path


if __name__ == "__main__":
    import sys
    rep = json.loads(Path(sys.argv[1]).read_text())
    print(render(rep, sys.argv[2] if len(sys.argv) > 2 else "report.html"))

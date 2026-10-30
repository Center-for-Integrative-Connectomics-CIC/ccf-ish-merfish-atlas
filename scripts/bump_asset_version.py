"""Stamp one version on every local JS/CSS reference in docs/ so browsers never mix old and new files.

Run after editing anything in docs/assets/:  python scripts/bump_asset_version.py
"""

from __future__ import annotations

import re
import time
from pathlib import Path

DOCS = Path(__file__).resolve().parents[1] / "docs"
REF = re.compile(r"""(['"])((?:\.{1,2}/|assets/)[^'"?]*?\.(?:js|css))(?:\?v=[\w.]+)?\1""")


def main() -> None:
    v = time.strftime("%Y%m%d%H%M%S")
    for f in [*DOCS.rglob("*.html"), *(DOCS / "assets" / "js").glob("*.js")]:
        s = f.read_text()
        t = REF.sub(lambda m: f"{m[1]}{m[2]}?v={v}{m[1]}", s)
        if t != s:
            f.write_text(t)
            print(f"{f.relative_to(DOCS)}: v={v}")


if __name__ == "__main__":
    main()

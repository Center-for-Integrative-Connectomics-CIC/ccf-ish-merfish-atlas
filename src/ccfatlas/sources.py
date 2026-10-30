"""Download helpers for the original sources (resumable, atomic)."""

from __future__ import annotations

import json
import shutil
import urllib.parse
import urllib.request
from pathlib import Path

from ccfatlas import config as C


def fetch(url: str, dest: Path, chunk: int = 1 << 22) -> Path:
    """Download url -> dest (skips existing files; resumes partial downloads with Range)."""
    if dest.exists():
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_suffix(dest.suffix + ".part")
    have = part.stat().st_size if part.exists() else 0
    req = urllib.request.Request(url, headers={"Range": f"bytes={have}-"} if have else {})
    with urllib.request.urlopen(req, timeout=300) as r, open(part, "ab" if have else "wb") as f:
        shutil.copyfileobj(r, f, chunk)
    part.rename(dest)
    return dest


def abc(rel: str) -> Path:
    """An Allen Brain Cell Atlas file, mirrored under RAW/abc/ with its S3 key."""
    return fetch(f"{C.ABC_S3}/{rel}", C.RAW / "abc" / rel)


def ish_gene_list() -> list[dict]:
    """All ISH genes in the EBRAINS bucket: [{gene, bytes, last_modified}]."""
    out, marker = [], None
    while True:
        q = {"limit": 1000, "prefix": C.ISH_PREFIX}
        if marker:
            q["marker"] = marker
        objs = json.load(urllib.request.urlopen(f"{C.EBRAINS_BUCKET}?{urllib.parse.urlencode(q)}", timeout=120))["objects"]
        if not objs:
            return out
        for o in objs:
            if o["name"].endswith(C.ISH_SUFFIX):
                out.append({"gene": o["name"][len(C.ISH_PREFIX):-len(C.ISH_SUFFIX)], "bytes": o["bytes"],
                            "last_modified": o["last_modified"]})
        marker = objs[-1]["name"]


def ish_url(gene: str) -> str:
    return f"{C.ISH_DIRECT}{C.ISH_PREFIX}{urllib.parse.quote(gene)}{C.ISH_SUFFIX}"


def ish_volume(gene: str, cache: Path | None = None) -> bytes:
    path = (cache or C.RAW / "ish") / f"{gene}{C.ISH_SUFFIX}"
    return fetch(ish_url(gene), path).read_bytes()

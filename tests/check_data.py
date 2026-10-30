"""Integrity checks for the baked data (run after any rebuild):  python tests/check_data.py

Checks manifests against the arrays they describe, cell positions inside the grid, cell-type and region
indices, region tables, the ISH orientation (Prox1 peaks in the dentate gyrus) and the MERFISH expression
files (one value per cell, expected gene set) for each brain.
"""

from __future__ import annotations

import gzip
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from ccfatlas import config as C  # noqa: E402

D = C.SITE_DATA
fails = []


def check(ok, msg):
    print(("ok   " if ok else "FAIL ") + msg)
    if not ok:
        fails.append(msg)


def arr(path, dtype):
    return np.frombuffer(gzip.decompress(Path(path).read_bytes()), dtype)


def main() -> None:
    onto = json.loads((D / "ref/ontology.json").read_text())
    n_reg = len(onto) + 1
    tax = json.loads((D / "merfish/taxonomy.json").read_text())
    n_cl = len(tax["cluster"])
    ann = arr(D / "ref/annotation.u16.gz", "<u2")
    check(ann.size == np.prod(C.GRID), f"annotation size {ann.size}")
    check(ann.max() < n_reg, "annotation indices within the ontology")
    for lv, par in (("cluster", "supertype"), ("supertype", "subclass"), ("subclass", "class")):
        check(all(0 <= e["parent"] < len(tax[par]) for e in tax[lv]), f"taxonomy {lv} -> {par} parents valid")
    ds = json.loads((D / "datasets.json").read_text())
    check(len(ds["merfish"]) == 5, "five MERFISH brains listed")
    lim = np.array(C.GRID) * C.VOXEL_UM
    for d in ds["merfish"]:
        k = d["key"]
        man = json.loads((D / f"merfish/{k}/manifest.json").read_text())
        n = man["cells"]
        pos = arr(D / f"merfish/{k}/pos.u16.gz", "<u2").reshape(-1, 3)
        check(len(pos) == n, f"{k}: {n:,} positions")
        check(bool((pos < lim).all()), f"{k}: positions inside the CCFv3BBP grid")
        check(bool((np.diff(pos[:, 0].astype(int)) >= 0).all()), f"{k}: cells sorted by AP")
        cl = arr(D / f"merfish/{k}/cluster.u16.gz", "<u2")
        check(len(cl) == n and cl.max() < n_cl, f"{k}: cluster indices valid")
        for src in ("ccf", "bbp"):
            r = arr(D / f"merfish/{k}/region_{src}.u16.gz", "<u2")
            check(len(r) == n and r.max() < n_reg, f"{k}: region_{src} valid ({(r > 0).mean():.1%} labelled)")
            rc = man["files"][f"region_cluster_{src}"]
            reg = arr(D / f"merfish/{k}/{rc['region']['file']}", "<u2")
            cnt = arr(D / f"merfish/{k}/{rc['count']['file']}", "<u4")
            root = [o["i"] for o in onto if o["id"] == 997][0]
            check(int(cnt[reg == root].sum()) == int((r > 0).sum()), f"{k}: region_cluster_{src} root total = labelled cells")
            check(bool((np.diff(reg.astype(int)) >= 0).all()), f"{k}: region_cluster_{src} sorted by region")
        genes = man.get("genes", [])
        check(len(genes) == d["genes"], f"{k}: {len(genes)} genes in the manifest")
        exp = C.HF_STAGE / "merfish" / k / "expr"
        files = {f.name[:-6] for f in exp.glob("*.u8.gz")}
        check(files == {g["symbol"] for g in genes}, f"{k}: expression files match the gene list")
        g0 = genes[0]["symbol"]
        check(arr(exp / f"{g0}.u8.gz", "u1").size == n, f"{k}: expression file has one value per cell")
    ish = json.loads((D / "ish/genes.json").read_text())
    G = len(ish["genes"])
    check(G == 4085, f"ISH genes {G}")
    mean = arr(D / "ish/region_mean.f16.gz", "<f2").reshape(n_reg, G).astype(np.float32)
    names = [g["gene"] for g in ish["genes"]]
    if "Prox1" in names:
        j = names.index("Prox1")
        acr = {o["i"]: o["acronym"] for o in onto}
        top = [acr.get(i) for i in np.argsort(-np.nan_to_num(mean[:, j], nan=-1))[:3]]
        check("DG-sg" in top or "DG" in top, f"ISH orientation: Prox1 top regions {top}")
    conc = json.loads((D / "qc/concordance.json").read_text())
    check(len(conc["ish_vs_merfish"]) == 5, "concordance for all five brains")
    print(f"\n{len(fails)} failures" if fails else "\nall checks passed")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()

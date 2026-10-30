"""Step 2: how the MERFISH space (Allen CCFv3, CCF-2020 parcellation) maps onto the ISH space (CCFv3BBP).

The MERFISH cells carry Allen CCFv3 coordinates (ABC atlas, x = AP, y = DV, z = LR in mm on the CCFv3
grid); the ISH volumes live on CCFv3BBP (25 µm, 566 x 320 x 456), which extends CCFv3 anteriorly. This
step compares the two annotations whole-brain to fix and document the mapping:
  1. the CCF-2020 annotation the cells were registered with (10 µm) -> structure ids -> 25 µm grid;
  2. integer shift scan AP 8..20, DV/LR -2..2 voxels: brain-mask Dice and label agreement at
     several ontology depths;
  3. per-structure Dice at the best shift (where CCFv3BBP redraws or adds regions, cell region
     assignments differ between the two spaces -- the site shows both labels).
Output: docs/data/qc/alignment.json, data_work/ccfv3_on_bbp_grid.npy
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
import tifffile
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from ccfatlas import config as C  # noqa: E402
from ccfatlas.bake import write_json  # noqa: E402


def parcellation_to_structure() -> pd.Series:
    """ABC parcellation_index -> CCF structure id (the 'substructure' level label)."""
    mem = pd.read_csv(C.RAW / "abc" / C.CCF2020 / "parcellation_to_parcellation_term_membership.csv")
    lab = mem.drop_duplicates("parcellation_index").set_index("parcellation_index").parcellation_label
    return lab.str.rsplit("-", n=1).str[1].astype(int)


def ccfv3_on_25(lut: pd.Series) -> np.ndarray:
    img = nib.load(C.RAW / "abc/image_volumes/Allen-CCF-2020/20230630/annotation_10.nii.gz")
    a10 = np.asarray(img.dataobj)                                    # (AP, DV, LR), P I R codes = ASR index order
    l = np.zeros(int(a10.max()) + 1, np.uint32)
    l[lut.index.values[lut.index.values <= a10.max()]] = lut.values[lut.index.values <= a10.max()]
    shape25 = tuple(int(np.floor((n - 1) * 0.4)) + 1 for n in a10.shape)   # centres at 25 µm steps
    idx = [np.round(np.arange(n) * 2.5).astype(int) for n in shape25]
    return l[a10[np.ix_(*idx)]]


def main() -> None:
    bbp_dir = C.RAW / "bbp" / C.BBP_PACKAGE
    bbp = tifffile.imread(bbp_dir / "annotation.tiff").astype(np.int64)
    structs = json.loads((bbp_dir / "structures.json").read_text())
    path = {s["id"]: s["structure_id_path"] for s in structs}
    lut = parcellation_to_structure()
    cc = ccfv3_on_25(lut).astype(np.int64)
    C.WORK.mkdir(parents=True, exist_ok=True)
    print("CCFv3 (25 µm)", cc.shape, "CCFv3BBP", bbp.shape)

    ids_cc = np.unique(cc)
    missing = [int(i) for i in ids_cc if i and i not in path]
    print(f"CCFv3 structure ids absent from the CCFv3BBP ontology: {len(missing)} {missing[:12]}")

    def ancestor_map(ids, depth):
        out = np.zeros(int(max(ids.max(), bbp.max())) + 1, np.int64)
        for i in ids:
            if i in path:
                p = path[i]
                out[i] = p[min(depth, len(p) - 1)]
        return out

    allids = np.union1d(ids_cc, np.unique(bbp))
    depth_maps = {d: ancestor_map(allids, d) for d in (3, 5, 7)}

    def score(shift):
        da, dd, dl = shift
        # place the CCFv3 grid on the CCFv3BBP grid: BBP index = CCFv3 index + shift
        out = np.zeros_like(bbp)
        src = [slice(max(0, -s), min(n, N - s)) for s, n, N in zip(shift, cc.shape, bbp.shape)]
        dst = [slice(max(0, s), max(0, s) + (sl.stop - sl.start)) for s, sl in zip(shift, src)]
        out[tuple(dst)] = cc[tuple(src)]
        a, b = out > 0, bbp > 0
        res = {"shift": list(shift), "brain_dice": float(2 * (a & b).sum() / (a.sum() + b.sum()))}
        both = a & b
        res["exact"] = float((out[both] == bbp[both]).mean())
        for d, m in depth_maps.items():
            x, y = m[out[both]], m[bbp[both]]
            ok = (x > 0) & (y > 0)
            res[f"depth{d}"] = float((x[ok] == y[ok]).mean())
        return res, out

    rows = []
    for da in range(8, 21, 2):
        rows.append(score((da, 0, 0))[0])
    best_ap = max(rows, key=lambda r: r["depth5"])["shift"][0]
    for da in range(best_ap - 1, best_ap + 2):
        for dd in (-2, -1, 0, 1, 2):
            for dl in (-2, -1, 0, 1, 2):
                if (da, dd, dl) != (best_ap, 0, 0) or True:
                    rows.append(score((da, dd, dl))[0])
    scan = pd.DataFrame(rows).drop_duplicates(subset=None if False else ["brain_dice", "exact", "depth5"])
    best = max(rows, key=lambda r: (r["depth5"], r["brain_dice"]))
    print(pd.DataFrame(rows).sort_values("depth5", ascending=False).head(8).to_string(index=False))
    print("best:", best)
    res, placed = score(tuple(best["shift"]))
    np.save(C.WORK / "ccfv3_on_bbp_grid.npy", placed.astype(np.uint32))

    # per-structure Dice (structure = itself + descendants) at the best shift
    desc = {}
    for s in structs:
        for anc in s["structure_id_path"]:
            desc.setdefault(anc, []).append(s["id"])
    per = []
    for s in structs:
        d = np.array(desc[s["id"]])
        a = np.isin(placed, d)
        b = np.isin(bbp, d)
        nb, na = int(b.sum()), int(a.sum())
        if nb + na < 200:
            continue
        per.append({"id": s["id"], "acronym": s["acronym"], "name": s["name"], "depth": len(s["structure_id_path"]) - 1,
                    "vox_bbp": nb, "vox_ccfv3": na, "dice": float(2 * (a & b).sum() / (na + nb))})
    per = pd.DataFrame(per)
    C.WORK.joinpath("alignment_per_structure.csv").write_text(per.to_csv(index=False))
    summ = {"best_shift_vox": best["shift"], "best": best, "scan": rows,
            "per_structure": per.to_dict(orient="records"),
            "ccfv3_ids_absent_from_bbp": missing,
            "note": "CCFv3BBP index = CCFv3 index + shift (AP, DV, LR) at 25 µm"}
    write_json(C.SITE_DATA / "qc" / "alignment.json", summ)
    big = per[per.vox_bbp > 5000]
    print(f"structures with >5000 voxels: {len(big)}; Dice median {big.dice.median():.3f}; "
          f"< 0.8: {(big.dice < 0.8).sum()}; < 0.5: {(big.dice < 0.5).sum()}")
    print(big.sort_values("dice").head(15)[["acronym", "name", "vox_bbp", "vox_ccfv3", "dice"]].to_string(index=False))


if __name__ == "__main__":
    main()

"""Step 8: MERFISH cells in the 3D ARA, for the HGEA-ARA tab, with an accuracy check.

The HGEA-ARA atlas (the HGEA plates on their own slices of the 3D Allen Reference Atlas Nissl brain, CCFv2 2011)
and its 3D ARA -> CCFv3BBP transform come from the HGEA 3D ARA atlas
(https://github.com/Center-for-Integrative-Connectomics-CIC/hgea-ara3d-atals), placed in data_raw/hgea_ara3d/:
    transforms/   ara_to_ccfv3bbp_0.nii.gz, ara_to_ccfv3bbp_1.mat, transforms.json
    site/         its web data: nissl.u8.gz, ara.u16.gz, ontology.json, hgea_full.u8.gz, hgea_plates.u8.gz, meta.json, svg/

Cells: CCFv3BBP µm (this site's cell positions) -> 3D ARA µm by the same transform used for the atlas (for points
ANTs applies the forward list fixed -> moving; checked: CCFv3BBP DG/CA1/SUB voxels land in the ARA's DG/CA1/SUB
more often with it than with the AP offset alone or the inverse list). Only cells inside the region where the
non-rigid part was fitted (the hippocampal crop) are kept.

Check, per brain:
    region   each cell's own CCF field (CA1, CA2, CA3, DG, SUB+ProS; the data owners' registration) vs the 2011 ARA
             field at its mapped position, and vs the CCFv3BBP field at its CCF position for reference
    types    DG Glut, CA1-ProS Glut, CA3 Glut cells: share landing in the ARA's DG, CA1, CA3

Output: docs/data/hgea_ara/  the atlas files and the HGEA SVG plates (copied), merfish/<brain>/{pos.u16.gz (3D ARA µm: AP, DV, LR),
        cluster.u16.gz, section.u16.gz, hgea.u8.gz, idx.u32.gz (row in the brain's full cell list, for gene lookups)},
        cells.json, accuracy.json, ish_lookup.{ara,ccf}.u32.gz (3D ARA voxels around the hippocampus and the CCFv3BBP
        voxel each maps to, for showing the Blue Brain ISH in the ARA; ARA -> CCF uses the inverse list, round trip 0.14 µm)
        docs/data/hgea_ccf/  the HGEA and the ARA Nissl in CCFv3BBP (from the hgea-ara3d-atlas CCFv3BBP edition: site/ccfv3bbp/),
        for the overlay and side-by-side toggle
"""

from __future__ import annotations

import gzip
import json
import shutil
import sys
import time
from pathlib import Path

import ants
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from ccfatlas import config as C  # noqa: E402

SRC = C.RAW / "hgea_ara3d"
OUT = C.SITE_DATA / "hgea_ara"
FIELD_IDS = {"CA1": 382, "CA2": 423, "CA3": 463, "DG": 726, "SUB": 502, "ProS": 484682470}
FIELD_OF_ACR = {"CA1": "CA1", "CA2": "CA2", "CA3": "CA3", "DG": "DG", "SUB": "SUB"}   # 2011 acronym prefixes
TYPES = {"037 DG Glut": "DG", "016 CA1-ProS Glut": "CA1", "017 CA3 Glut": "CA3"}
V = 25.0
CHUNK = 400_000


def read_gz(path, dtype, shape=None):
    a = np.frombuffer(gzip.open(path).read(), dtype)
    return a.reshape(shape) if shape else a


def write_gz(path, a) -> int:
    with gzip.open(path, "wb", compresslevel=9) as f:
        f.write(np.ascontiguousarray(a).tobytes())
    return path.stat().st_size


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    T = json.loads((SRC / "transforms" / "transforms.json").read_text())
    tx = [str(SRC / "transforms" / n) for n in T["forward"]]
    lo, hi = np.array(T["crop_lo_vox"]) * V, np.array(T["crop_hi_vox"]) * V
    ap_origin = T["moving_origin_um"][0]

    # the atlas files, unchanged
    for f in ("nissl.u8.gz", "ara.u16.gz", "ontology.json", "hgea_full.u8.gz", "hgea_plates.u8.gz", "meta.json"):
        shutil.copy(SRC / "site" / f, OUT / f)
    shutil.copytree(SRC / "site" / "svg", OUT / "svg", dirs_exist_ok=True)
    meta = json.loads((OUT / "meta.json").read_text())
    meta.update(title="HGEA-ARA · MERFISH", subtitle="MERFISH cells mapped into the 3D Allen Reference Atlas with the HGEA domains (Bienkowski et al. 2018)")
    (OUT / "meta.json").write_text(json.dumps(meta, indent=1))
    shape = meta["grid"]["shape"]
    ara_lab = read_gz(OUT / "ara.u16.gz", np.uint16, shape)
    hgea = read_gz(OUT / "hgea_full.u8.gz", np.uint8, shape)
    ara_onto = json.loads((OUT / "ontology.json").read_text())
    ara_field = np.array([next((f for p, f in FIELD_OF_ACR.items() if e["acronym"].startswith(p)), "") for e in ara_onto], dtype=object)

    onto = json.loads((C.SITE_DATA / "ref" / "ontology.json").read_text())
    parent = {o["i"]: o["parent"] for o in onto}
    fid = {o["i"]: o["id"] for o in onto}
    n_onto = max(parent) + 1

    def field_of_index():
        out = np.full(n_onto, "", dtype=object)
        for i in range(n_onto):
            a = i
            while a in parent and a:
                f = next((k for k, v in FIELD_IDS.items() if fid.get(a) == v), None)
                if f:
                    out[i] = "SUB" if f == "ProS" else f
                    break
                if parent[a] == a:
                    break
                a = parent[a]
        return out
    ccf_field = field_of_index()
    tax = json.loads((C.SITE_DATA / "merfish" / "taxonomy.json").read_text())
    sub_of_cluster = np.array([tax["supertype"][c["parent"]]["parent"] for c in tax["cluster"]])
    sub_names = [s["name"] for s in tax["subclass"]]
    bbp_ann = None

    cells, acc = {}, {}
    for d in json.loads((C.SITE_DATA / "datasets.json").read_text())["merfish"]:
        k = d["key"]
        man = json.loads((C.SITE_DATA / "merfish" / k / "manifest.json").read_text())
        f, n = man["files"], man["cells"]
        base = C.SITE_DATA / "merfish" / k
        pos = read_gz(base / f["pos"]["file"], np.uint16, (n, 3)).astype(np.float64)
        keep = np.all((pos >= lo) & (pos <= hi), axis=1)
        idx = np.nonzero(keep)[0]
        w = np.empty((len(idx), 3))
        for s in range(0, len(idx), CHUNK):
            p = pd.DataFrame(pos[idx[s:s + CHUNK]], columns=["x", "y", "z"])
            w[s:s + CHUNK] = ants.apply_transforms_to_points(3, p, tx).to_numpy()
        ara_um = np.column_stack([w[:, 0] - ap_origin, w[:, 1], w[:, 2]])
        vox = np.rint(ara_um / V).astype(int)
        inside = np.all((vox >= 0) & (vox < np.array(shape)), axis=1)
        idx, ara_um, vox = idx[inside], ara_um[inside], vox[inside]
        cluster = read_gz(base / f["cluster"]["file"], np.uint16)[idx]
        section = read_gz(base / f["section"]["file"], np.uint16)[idx]
        rccf = read_gz(base / f["region_ccf"]["file"], np.uint16)[idx]
        dom = hgea[tuple(vox.T)]
        own = ccf_field[rccf]
        landed = ara_field[ara_lab[tuple(vox.T)]]
        hp = own != ""
        res = dict(cells_in_crop=int(len(idx)), hippocampal_cells=int(hp.sum()),
                   region_agreement=float((own[hp] == landed[hp]).mean()),
                   per_field={fl: float((landed[own == fl] == fl).mean()) for fl in ("CA1", "CA2", "CA3", "DG", "SUB") if (own == fl).any()},
                   cell_types={})
        sub = sub_of_cluster[cluster]
        for tname, fl in TYPES.items():
            m = sub == sub_names.index(tname)
            if m.any():
                res["cell_types"][tname] = dict(cells=int(m.sum()), in_ara_field=float((landed[m] == fl).mean()),
                                               own_ccf_label_in_field=float((own[m] == fl).mean()))
        acc[k] = res
        bd = OUT / "merfish" / k
        bd.mkdir(parents=True, exist_ok=True)
        files = dict(pos=write_gz(bd / "pos.u16.gz", np.clip(np.rint(ara_um), 0, 65535).astype(np.uint16)),
                     cluster=write_gz(bd / "cluster.u16.gz", cluster), section=write_gz(bd / "section.u16.gz", section),
                     hgea=write_gz(bd / "hgea.u8.gz", dom), idx=write_gz(bd / "idx.u32.gz", idx.astype(np.uint32)))
        cells[k] = dict(cells=int(len(idx)), label=d["label"], short=d["short"], bytes_gz=files)
        print(f"{k}: {len(idx):,} cells mapped; region agreement {res['region_agreement']:.3f} "
              f"{ {a: round(b, 3) for a, b in res['per_field'].items()} }; types "
              f"{ {a: round(b['in_ara_field'], 3) for a, b in res['cell_types'].items()} } ({time.time() - t0:.0f}s)", flush=True)
    (OUT / "cells.json").write_text(json.dumps(cells, indent=1))

    # ISH lookup: ARA voxels around the hippocampus -> CCFv3BBP voxel (nearest), ARA -> CCF with the inverse list
    from scipy import ndimage
    region = ndimage.binary_dilation((hgea > 0) | np.isin(ara_lab, [i for i, e in enumerate(ara_onto) if ara_field[i]]), iterations=8)
    av = np.argwhere(region)
    inv = [str(SRC / "transforms" / n) for n in T["inverse"]]
    cc = np.empty(len(av), np.uint32)
    for s in range(0, len(av), CHUNK):
        a = av[s:s + CHUNK].astype(np.float64) * V
        a[:, 0] += ap_origin
        w = ants.apply_transforms_to_points(3, pd.DataFrame(a, columns=["x", "y", "z"]), inv, whichtoinvert=[True, False]).to_numpy()
        v = np.clip(np.rint(w / V).astype(int), 0, np.array(C.GRID) - 1)
        cc[s:s + CHUNK] = np.ravel_multi_index(tuple(v.T), C.GRID)
    write_gz(OUT / "ish_lookup.ara.u32.gz", np.ravel_multi_index(tuple(av.T), shape).astype(np.uint32))
    write_gz(OUT / "ish_lookup.ccf.u32.gz", cc)
    print(f"ISH lookup: {len(av):,} ARA voxels", flush=True)

    # HGEA + ARA Nissl in CCFv3BBP, for the overlay / side-by-side toggle
    HC = C.SITE_DATA / "hgea_ccf"
    HC.mkdir(exist_ok=True)
    for f in ("nissl.u8.gz", "hgea_full.u8.gz", "meta.json"):
        shutil.copy(SRC / "site" / "ccfv3bbp" / f, HC / f)
    (OUT / "accuracy.json").write_text(json.dumps(acc, indent=1))


if __name__ == "__main__":
    main()

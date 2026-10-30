"""Step 3: bake each MERFISH brain separately for the site and the Hugging Face data store.

  python scripts/03_bake_merfish.py [--datasets zhuang4,...] [--no-expr]

Per brain, into docs/data/merfish/<key>/ (GitHub Pages):
  pos.u16.gz          (n, 3) CCFv3BBP µm (AP, DV, LR), cells sorted by AP (fast coronal slab queries)
  cluster.u16.gz      WMB cluster index (taxonomy.json), class/subclass/supertype follow from it
  region_ccf.u16.gz   the cell's own Allen CCF-2020 parcellation, as a CCFv3BBP ontology index
  region_bbp.u16.gz   the CCFv3BBP structure at the cell's position
  section.u16.gz      original brain section (index into manifest.sections)
  quality.u8.gz       mapping quality x 255 (Yao: average_correlation_score; Zhuang: cluster_confidence_score)
  region_cluster.*.gz cell counts per (region incl. ancestors, cluster), for each label source
  region_gene_{mean,frac}.f16.gz  per region (incl. ancestors, own labels) x gene: mean log2 expression and
                      fraction of cells with expression > 0 -- full precision float32 sums, stored float16
  manifest.json
and into hf_dataset/merfish/<key>/expr/<gene>.u8.gz the per-cell log2 expression of every panel gene, same
cell order, quantised to uint8 with value = q / 255 * vmax (vmax per gene in the manifest; 0 stays 0).
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from ccfatlas import config as C  # noqa: E402
from ccfatlas import merfish as M  # noqa: E402
from ccfatlas.bake import write_array, write_json  # noqa: E402

BAND = 50_000


class Bands:
    """Row bands of an h5ad X as dense float32, for dense datasets and CSR groups alike."""

    def __init__(self, X):
        self.X = X
        self.sparse = isinstance(X, h5py.Group)
        if self.sparse:
            self.shape = tuple(int(v) for v in X.attrs["shape"])
            self.indptr = X["indptr"][:]
        else:
            self.shape = X.shape

    def __call__(self, r0: int, r1: int) -> np.ndarray:
        r1 = min(r1, self.shape[0])
        if not self.sparse:
            return self.X[r0:r1]
        a, b = int(self.indptr[r0]), int(self.indptr[r1])
        m = sparse.csr_matrix((self.X["data"][a:b], self.X["indices"][a:b], self.indptr[r0:r1 + 1] - a),
                              shape=(r1 - r0, self.shape[1]))
        return m.toarray().astype(np.float32)


def ancestors_matrix(onto: list[dict]) -> sparse.csr_matrix:
    """A[a, r] = 1 if a is r or an ancestor of r (ontology indices, 0 = none)."""
    n = len(onto) + 1
    parent = np.zeros(n, int)
    for o in onto:
        parent[o["i"]] = o["parent"]
    rows, cols = [], []
    for r in range(1, n):
        a = r
        while a:
            rows.append(a)
            cols.append(r)
            a = parent[a]
    return sparse.csr_matrix((np.ones(len(rows), np.float64), (rows, cols)), shape=(n, n))


def depth_ancestor(onto, depth):
    n = len(onto) + 1
    parent = np.zeros(n, int)
    dep = np.zeros(n, int)
    for o in onto:
        parent[o["i"]], dep[o["i"]] = o["parent"], o["depth"]
    out = np.arange(n)
    for r in range(1, n):
        a = r
        while dep[a] > depth and parent[a]:
            a = parent[a]
        out[r] = a
    return out


def bake(m: C.Merfish, tax: dict, onto: list[dict], ann: np.ndarray, A, expr: bool) -> dict:
    t0 = time.time()
    out = C.SITE_DATA / "merfish" / m.key
    cells = M.load_cells(m)
    n_all = len(cells)
    xyz = cells[["x", "y", "z"]].to_numpy(float)
    um = np.column_stack([xyz[:, 0] * 1000 + C.CCFV3_AP_OFFSET_VOX * C.VOXEL_UM, xyz[:, 1] * 1000, xyz[:, 2] * 1000])
    lim = np.array(C.GRID) * C.VOXEL_UM
    ok = np.isfinite(um).all(1) & (um >= 0).all(1) & (um < lim - 0.5).all(1)
    alias_ok = cells.cluster_alias.isin(tax["alias_to_cluster"].keys()).to_numpy()
    keep = ok & alias_ok
    cells, um = cells[keep].reset_index(drop=True), um[keep]
    order = np.argsort(um[:, 0], kind="stable")
    cells, um = cells.iloc[order].reset_index(drop=True), um[order]
    n = len(cells)

    id_to_i = {o["id"]: o["i"] for o in onto}
    lut = M.parcellation_to_structure()
    sid = cells.parcellation_index.map(lut).fillna(0).astype(int).to_numpy()
    region_ccf = np.array([id_to_i.get(s, 0) for s in sid], np.uint16)
    vox = np.round(um / C.VOXEL_UM).astype(int)
    vox = np.minimum(vox, np.array(C.GRID) - 1)
    region_bbp = ann[tuple(vox.T)]
    cluster = cells.cluster_alias.map(tax["alias_to_cluster"]).astype(np.uint16).to_numpy()
    sections = sorted(cells.brain_section_label.astype(str).unique())
    section = cells.brain_section_label.astype(str).map({s: i for i, s in enumerate(sections)}).astype(np.uint16).to_numpy()
    qcol = "average_correlation_score" if "average_correlation_score" in cells else "cluster_confidence_score"
    quality = np.clip(np.nan_to_num(cells[qcol].to_numpy(float)) * 255, 0, 255).astype(np.uint8)

    both = (region_ccf > 0) & (region_bbp > 0)
    d5 = depth_ancestor(onto, 5)
    agree = {"exact": float((region_ccf[both] == region_bbp[both]).mean()),
             "depth5": float((d5[region_ccf[both]] == d5[region_bbp[both]]).mean()),
             "cells_with_ccf_label": int((region_ccf > 0).sum()), "cells_with_bbp_label": int((region_bbp > 0).sum())}

    man = {"key": m.key, "dataset": m.name, "label": m.label, "reference": m.reference,
           "source": {k: m.url(getattr(m, k)) for k in ("expr", "meta", "genes", "ccf")},
           "cells": n, "cells_in_source": n_all, "dropped_outside_grid": int((~ok).sum()),
           "dropped_unknown_cluster": int((ok & ~alias_ok).sum()), "sections": sections,
           "quality_field": qcol, "region_label_agreement": agree,
           "position": "CCFv3BBP µm (AP, DV, LR) = (x*1000 + 350, y*1000, z*1000) from ABC CCFv3 coordinates"}
    man["files"] = {"pos": write_array(out / "pos.u16.gz", np.round(um).astype(np.uint16)),
                    "cluster": write_array(out / "cluster.u16.gz", cluster),
                    "region_ccf": write_array(out / "region_ccf.u16.gz", region_ccf),
                    "region_bbp": write_array(out / "region_bbp.u16.gz", region_bbp.astype(np.uint16)),
                    "section": write_array(out / "section.u16.gz", section),
                    "quality": write_array(out / "quality.u8.gz", quality)}
    if "high_quality_transfer" in cells:
        man["files"]["high_quality"] = write_array(out / "high_quality.u8.gz",
                                                    cells.high_quality_transfer.astype(bool).to_numpy().astype(np.uint8))

    n_reg, n_cl = len(onto) + 1, len(tax["cluster"])
    for src, reg in (("ccf", region_ccf), ("bbp", region_bbp)):
        leaf = sparse.csr_matrix((np.ones(n), (reg.astype(int), cluster.astype(int))), shape=(n_reg, n_cl))
        full = (A @ leaf).tocoo()
        man["files"][f"region_cluster_{src}"] = {
            "region": write_array(out / f"region_cluster_{src}.region.u16.gz", full.row.astype(np.uint16)),
            "cluster": write_array(out / f"region_cluster_{src}.cluster.u16.gz", full.col.astype(np.uint16)),
            "count": write_array(out / f"region_cluster_{src}.count.u32.gz", full.data.astype(np.uint32))}
    print(f"  {m.key}: {n:,} cells ({n_all - n:,} dropped), {len(sections)} sections, "
          f"label agreement exact {agree['exact']:.3f} / depth-5 {agree['depth5']:.3f} [{time.time() - t0:.0f}s]", flush=True)

    if expr:
        f = h5py.File(C.RAW / "abc" / m.expr, "r")
        X = Bands(f["X"])
        obs = pd.Index(f["obs"]["cell_label"][:].astype(str))
        sym = f["var"]["gene_symbol"][:].astype(str)
        ens = f["var"]["gene_identifier"][:].astype(str)
        rows = obs.get_indexer(cells.cell_label.astype(str))
        has = rows >= 0
        man["cells_without_expression"] = int((~has).sum())
        n_g = X.shape[1]
        # pass 1: per-gene max and region sums (float64) over own CCF labels
        vmax = np.zeros(n_g, np.float64)
        sums = np.zeros((n_reg, n_g))
        npos = np.zeros((n_reg, n_g))
        ncell = np.bincount(region_ccf.astype(int), minlength=n_reg).astype(np.float64)
        cell_of_row = np.full(X.shape[0], -1, np.int64)
        cell_of_row[rows[has]] = np.nonzero(has)[0]
        Q = np.zeros((n_g, n), np.uint8)
        bands = list(range(0, X.shape[0], BAND))
        for pass_ in (1, 2):
            for b0 in bands:
                blk = X(b0, b0 + BAND)
                ci = cell_of_row[b0:b0 + len(blk)]
                sel = ci >= 0
                blk, ci = blk[sel], ci[sel]
                if pass_ == 1:
                    vmax = np.maximum(vmax, blk.max(0, initial=0))
                    oh = sparse.csr_matrix((np.ones(len(ci)), (region_ccf[ci].astype(int), np.arange(len(ci)))), shape=(n_reg, len(ci)))
                    sums += oh @ blk.astype(np.float64)
                    npos += oh @ (blk > 0).astype(np.float64)
                else:
                    q = np.ceil(blk / np.maximum(vmax, 1e-9)[None, :] * 255.0)
                    Q[:, ci] = np.where(blk > 0, np.clip(q, 1, 255), 0).astype(np.uint8).T   # quantize_u8, vectorised
            print(f"    expression pass {pass_} done [{time.time() - t0:.0f}s]", flush=True)
        f.close()
        S, P, N = A @ sums, A @ npos, A @ ncell
        mean = np.where(N[:, None] > 0, S / np.maximum(N[:, None], 1), np.nan).astype(np.float16)
        frac = np.where(N[:, None] > 0, P / np.maximum(N[:, None], 1), np.nan).astype(np.float16)
        man["files"]["region_gene_mean"] = write_array(out / "region_gene_mean.f16.gz", mean)
        man["files"]["region_gene_frac"] = write_array(out / "region_gene_frac.f16.gz", frac)
        man["region_cells"] = N.astype(int).tolist()
        hf = C.HF_STAGE / "merfish" / m.key / "expr"
        hf.mkdir(parents=True, exist_ok=True)
        genes = []
        seen = {}
        for g in range(n_g):
            name = sym[g] if sym[g] not in seen else f"{sym[g]}__{ens[g]}"
            seen[sym[g]] = True
            info = write_array(hf / f"{name}.u8.gz", Q[g])
            genes.append({"symbol": name, "ensembl": ens[g], "vmax": float(vmax[g]), "bytes_gz": info["bytes_gz"],
                          "detected_frac": float((Q[g] > 0).mean())})
        man["genes"] = genes
        man["expr_encoding"] = "uint8, value = q / 255 * vmax (log2 expression, ABC 'log2' matrix); 0 = not detected"
        print(f"    {n_g} genes written, {sum(x['bytes_gz'] for x in genes) / 1e6:.0f} MB [{time.time() - t0:.0f}s]", flush=True)
    write_json(out / "manifest.json", man, compact=False)
    return man


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--datasets", default=",".join(m.key for m in C.MERFISH))
    ap.add_argument("--no-expr", action="store_true")
    args = ap.parse_args()
    import tifffile
    onto = json.loads((C.SITE_DATA / "ref" / "ontology.json").read_text())
    import gzip
    ann = np.frombuffer(gzip.decompress((C.SITE_DATA / "ref" / "annotation.u16.gz").read_bytes()), "<u2").reshape(C.GRID)
    tax = M.taxonomy()
    write_json(C.SITE_DATA / "merfish" / "taxonomy.json", {k: v for k, v in tax.items() if k != "alias_to_cluster"})
    A = ancestors_matrix(onto)
    for key in args.datasets.split(","):
        bake(C.MERFISH_BY_KEY[key], tax, onto, ann, A, not args.no_expr)


if __name__ == "__main__":
    main()

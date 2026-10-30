"""MERFISH datasets (Allen Brain Cell Atlas): cells, cell types, regions and expression, one brain at a time.

Positions: the ABC atlas gives each cell Allen CCFv3 coordinates (mm; x = AP, y = DV, z = LR, on the
CCFv3 grid). On the CCFv3BBP grid: AP_um = x * 1000 + 14 * 25, DV_um = y * 1000, LR_um = z * 1000
(step 2 verified the 14-voxel AP offset whole-brain). Positions are kept at 1 µm (uint16).

Regions per cell, both kept:
  region_ccf  the cell's own Allen CCF-2020 parcellation from the ABC atlas (the data owners' registration)
  region_bbp  the CCFv3BBP structure at the cell's position (the ISH atlas' own labels)
Cell types: WMB taxonomy (class > subclass > supertype > cluster), shared by all five brains.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from ccfatlas import config as C

LEVELS = ["class", "subclass", "supertype", "cluster"]


def taxonomy() -> dict:
    base = C.RAW / "abc" / C.TAXONOMY / "views"
    piv = pd.read_csv(base / "cluster_to_cluster_annotation_membership_pivoted.csv")
    col = pd.read_csv(base / "cluster_to_cluster_annotation_membership_color.csv").set_index("cluster_alias")
    piv = piv.sort_values("cluster_alias").reset_index(drop=True)
    out = {"levels": LEVELS, "source": f"Allen Brain Cell Atlas {C.TAXONOMY} (Yao et al. 2023)"}
    index = {}
    for lv in LEVELS:
        names = sorted(piv[lv].dropna().unique())
        index[lv] = {n: i for i, n in enumerate(names)}
        colors = piv.drop_duplicates(lv).set_index(lv)["cluster_alias"].map(col[f"{lv}_color"])
        out[lv] = [{"name": n, "color": colors.get(n, "#888888")} for n in names]
    for lv, parent in zip(LEVELS[1:], LEVELS[:-1]):
        par = piv.drop_duplicates(lv).set_index(lv)[parent]
        for e in out[lv]:
            e["parent"] = index[parent][par[e["name"]]]
    nt = piv.drop_duplicates("cluster").set_index("cluster")["neurotransmitter"].fillna("")
    for e in out["cluster"]:
        e["nt"] = nt[e["name"]]
    out["alias_to_cluster"] = {int(a): index["cluster"][c] for a, c in zip(piv.cluster_alias, piv.cluster)}
    return out


def parcellation_to_structure() -> pd.Series:
    mem = pd.read_csv(C.RAW / "abc" / C.CCF2020 / "parcellation_to_parcellation_term_membership.csv")
    lab = mem.drop_duplicates("parcellation_index").set_index("parcellation_index").parcellation_label
    return lab.str.rsplit("-", n=1).str[1].astype(int)


def load_cells(m: C.Merfish) -> pd.DataFrame:
    ccf = pd.read_csv(C.RAW / "abc" / m.ccf)
    use = ["cell_label", "brain_section_label", "cluster_alias", "average_correlation_score",
           "cluster_confidence_score", "subclass_confidence_score", "high_quality_transfer"]
    meta = pd.read_csv(C.RAW / "abc" / m.meta, usecols=lambda c: c in use, dtype={"cell_label": str})
    ccf["cell_label"] = ccf.cell_label.astype(str)
    cells = ccf.merge(meta, on="cell_label", how="inner")
    return cells

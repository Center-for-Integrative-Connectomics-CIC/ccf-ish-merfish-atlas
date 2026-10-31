"""Step 5: site index + ISH-vs-MERFISH concordance (per gene, per brain).

datasets.json   what the viewer loads: the ISH atlas and the five MERFISH brains (labels, counts, panels)
qc/concordance.json
  For every gene measured by ISH and on a brain's MERFISH panel: Spearman correlation between the ISH
  region means and that brain's MERFISH region means (own CCF labels) across CCF regions with enough
  data in both (depth 5-8 structures, >= 50 cells and >= 500 voxels); Pearson on the same values (driven
  by where the gene is high, so it suits sparse markers whose many near-zero regions make ranks noisy);
  and peak agreement: the percentile rank, in each atlas, of the other atlas' top region (1.0 = same top
  region). As a reference ceiling, the same
  correlation between pairs of MERFISH brains for their shared genes (two independent animals, same
  modality). Region means are on different scales (ISH voxel value vs mean log2 per cell), so only rank
  agreement is compared.
"""

from __future__ import annotations

import gzip
import itertools
import json
import sys
from pathlib import Path

import numpy as np
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from ccfatlas import config as C  # noqa: E402
from ccfatlas.bake import write_json  # noqa: E402

SHORT = {"yao638850": "Allen 638850", "zhuang1": "Zhuang 1", "zhuang2": "Zhuang 2", "zhuang3": "Zhuang 3", "zhuang4": "Zhuang 4"}
MIN_CELLS, MIN_VOX = 50, 500


def load(path, shape):
    return np.frombuffer(gzip.decompress(Path(path).read_bytes()), "<f2").reshape(shape).astype(np.float32)


def main() -> None:
    onto = json.loads((C.SITE_DATA / "ref" / "ontology.json").read_text())
    depth = np.zeros(len(onto) + 1, int)
    for o in onto:
        depth[o["i"]] = o["depth"]
    ish = json.loads((C.SITE_DATA / "ish" / "genes.json").read_text())
    ish_genes = [g["gene"] for g in ish["genes"]]
    ish_mean = load(C.SITE_DATA / "ish" / ish["files"]["region_mean"]["file"], ish["files"]["region_mean"]["shape"])
    vox = np.array(ish["region_voxels"])

    mer = {}
    for m in C.MERFISH:
        d = C.SITE_DATA / "merfish" / m.key
        if not (d / "manifest.json").exists():
            continue
        man = json.loads((d / "manifest.json").read_text())
        mer[m.key] = {"man": man, "mean": load(d / man["files"]["region_gene_mean"]["file"], man["files"]["region_gene_mean"]["shape"]),
                      "cells": np.array(man["region_cells"]), "genes": [g["symbol"] for g in man["genes"]]}

    ds = {"ish": {"label": "Blue Brain ISH gene atlas (Carey, Piluso et al. 2026)", "genes": len(ish_genes),
                  "space": "CCFv3BBP 25 µm", "url_template": ish["url_template"]},
          "merfish": [{"key": k, "short": SHORT[k], "label": v["man"]["label"], "dataset": v["man"]["dataset"],
                       "reference": v["man"]["reference"], "cells": v["man"]["cells"], "genes": sum(not C.is_control_probe(g) for g in v["genes"]),
                       "sections": len(v["man"]["sections"]), "label_agreement": v["man"]["region_label_agreement"]}
                      for k, v in mer.items()]}
    write_json(C.SITE_DATA / "datasets.json", ds, compact=False)

    regs_ok = (depth >= 5) & (depth <= 8) & (vox >= MIN_VOX)
    gi_ish = {g: i for i, g in enumerate(ish_genes)}
    conc = {"min_cells": MIN_CELLS, "min_voxels": MIN_VOX, "depths": [5, 8], "ish_vs_merfish": {}, "merfish_vs_merfish": {}}
    for k, v in mer.items():
        rows = {}
        sel = regs_ok & (v["cells"] >= MIN_CELLS)
        for j, g in enumerate(v["genes"]):
            if g not in gi_ish:
                continue
            a, b = ish_mean[sel, gi_ish[g]], v["mean"][sel, j]
            ok = np.isfinite(a) & np.isfinite(b)
            if ok.sum() < 20 or np.ptp(a[ok]) == 0 or np.ptp(b[ok]) == 0:
                continue
            a, b = a[ok], b[ok]
            pr = lambda x, i: float((x < x[i]).mean() + 0.5 * (x == x[i]).mean())   # percentile rank of element i
            peak_ish_in_mer = pr(b, int(np.argmax(a)))
            peak_mer_in_ish = pr(a, int(np.argmax(b)))
            rows[g] = [round(float(stats.spearmanr(a, b)[0]), 4), int(ok.sum()), round(float(stats.pearsonr(a, b)[0]), 4),
                       round(peak_ish_in_mer, 4), round(peak_mer_in_ish, 4)]
        r = np.array([x[0] for x in rows.values()])
        pe = np.array([x[2] for x in rows.values()])
        pk = np.array([min(x[3], x[4]) for x in rows.values()])
        conc["ish_vs_merfish"][k] = {"genes": rows, "columns": ["spearman", "n_regions", "pearson", "peak_ish_in_mer", "peak_mer_in_ish"],
                                     "n_genes": len(rows), "median_rho": float(np.median(r)) if len(r) else None,
                                     "median_pearson": float(np.median(pe)) if len(pe) else None,
                                     "frac_rho_gt_0.5": float((r > 0.5).mean()) if len(r) else None,
                                     "frac_peaks_top5pct": float((pk >= 0.95).mean()) if len(pk) else None}
        print(f"ISH vs {k}: {len(rows)} shared genes, median Spearman {np.median(r):.3f}, Pearson {np.median(pe):.3f}, "
              f"both peaks in the other's top 5%: {(pk >= 0.95).mean():.0%}")
    for k1, k2 in itertools.combinations(mer, 2):
        v1, v2 = mer[k1], mer[k2]
        sel = regs_ok & (v1["cells"] >= MIN_CELLS) & (v2["cells"] >= MIN_CELLS)
        i2 = {g: i for i, g in enumerate(v2["genes"])}
        rs = []
        for j, g in enumerate(v1["genes"]):
            if g in i2:
                a, b = v1["mean"][sel, j], v2["mean"][sel, i2[g]]
                ok = np.isfinite(a) & np.isfinite(b)
                if ok.sum() >= 20 and np.ptp(a[ok]) > 0 and np.ptp(b[ok]) > 0:
                    rs.append(stats.spearmanr(a[ok], b[ok])[0])
        rs = np.array(rs)
        conc["merfish_vs_merfish"][f"{k1}|{k2}"] = {"n_genes": len(rs), "median_rho": float(np.median(rs)) if len(rs) else None}
        print(f"{k1} vs {k2}: {len(rs)} shared genes, median rho {np.median(rs):.3f}")
    write_json(C.SITE_DATA / "qc" / "concordance.json", conc)
    print("wrote datasets.json and qc/concordance.json")


if __name__ == "__main__":
    main()

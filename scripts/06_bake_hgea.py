"""Step 6: the HGEA gene set (Bienkowski et al. 2018) on the Blue Brain ISH atlas -> docs/data/hgea/.

Sources (data_raw/hgea/, copied from the paper's Supplementary Information, doi:10.1038/s41593-018-0241-y):
  supp_table2_gene_annotation.xlsx   Supplementary Table 2: binary annotation of 258 genes x 19 HGEA
                                     domains / SUB layers + "Putative HPF Interneurons" (the paper's definition)
  supp_table7_hgea_nomenclature.xlsx Supplementary Table 7: per HGEA domain (22, incl. connectivity-defined
                                     CA3id/ic, CA1dr/dc and the SUB subregions): genes, inputs, outputs
  bienkowski2018_natneurosci.pdf, supp_text_and_figures.pdf   for gene mentions per figure legend / text
Everything is kept exactly as published: names (with their asterisks), domain names, presence marks.
The paper and its supplement do not define the asterisk on 14 gene names; it is shown verbatim.

Derived here:
  - ISH volume for each gene: exact symbol, asterisk stripped, or a documented alias (old / renamed symbols);
    Epha5 and Fosl2 have no coronal ISH volume (Allen imaged them sagittally only).
  - Table 7 vs Table 2 consistency per domain.
  - figure/legend mentions of each gene in the main paper and the Supplementary Information.
  - per Table 2 column, two signature volumes on the hippocampus (CA, DG, SUB, ProS, FC, IG + 6 voxels):
      signature = mean over the column's marker genes of ISH / p99.5 (clipped to 1)
      contrast  = signature - mean over the genes the paper annotates in OTHER columns of the same group
                  but NOT in this one (absent marks), clipped at 0
    stored as uint8 (x 255) in a crop box (manifest gives the box). Interneuron genes form their own column.
"""

from __future__ import annotations

import concurrent.futures as cf
import gzip
import json
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from ccfatlas import config as C  # noqa: E402
from ccfatlas.bake import write_array, write_json  # noqa: E402

import importlib.util  # noqa: E402

_s = importlib.util.spec_from_file_location("ish4", Path(__file__).resolve().parent / "04_bake_ish.py")
_ish = importlib.util.module_from_spec(_s)
_s.loader.exec_module(_ish)

SRC = C.RAW / "hgea"
OUT = C.SITE_DATA / "hgea"
CITATION = ("Bienkowski M.S., Bowman I., Song M.Y. et al. Integration of gene expression and brain-wide connectivity "
            "reveals the multiscale organization of mouse hippocampal networks. Nat. Neurosci. 21, 1628-1643 (2018).")
DOI = "10.1038/s41593-018-0241-y"
# old or alternative names in the HGEA tables -> current symbol used by the Allen / Blue Brain volumes
ALIASES = {"Dscr1": "Rcan1", "Fam40b": "Strip2", "Igsf4d": "Cadm2", "Pcdh21": "Cdhr1", "Neuromedin B": "Nmb",
           "D930040M24Rik": "L3mbtl4", "Gm1399": "Gm1399",
           # lower-case l in the HGEA table; the Blue Brain volume is named "A730096I18Rik*" (capital I, asterisk)
           "A730096l18Rik": "A730096I18Rik*"}
GROUPS = [
    ("Dentate gyrus", ["DGd", "DGi", "DGv", "DG_pod", "DG_pov"]),
    ("CA3", ["CA3dd", "CA3d", "CA3i", "CA3v", "CA3vv"]),
    ("CA2", ["CA2"]),
    ("CA1", ["CA1d", "CA1i", "CA1v", "CA1vv"]),
    ("Subiculum layers", ["SUB_1", "SUB_2", "SUB_3", "SUB_4"]),
    ("Interneurons", ["Putative HPF Interneurons"]),
]
T7_TO_T2 = {"DGpod": "DG_pod", "DGpov": "DG_pov", "CA3id": "CA3i", "CA3ic": "CA3i", "CA1dr": "CA1d", "CA1dc": "CA1d"}
SUBREGION_LAYERS = {"SUBdd": [1, 4], "SUBdv": [1, 4], "ProSUB": [3, 4], "SUBv": [2, 3, 4], "SUBvv": [2, 3, 4]}


def clean(name: str) -> str:
    n = re.sub(r"\s*\(.*\)\s*", "", str(name)).strip().rstrip("*").strip()
    return n


def paren_alias(name: str) -> str | None:
    m = re.search(r"\(([^)]+)\)", str(name))
    return m.group(1).strip() if m else None


def table7():
    """One entry per domain label (22). A domain is labelled over its Inputs column; its genes are the gene
    columns before it. Neighbouring SUB subregions share gene columns in the table (SUBdd/SUBdv: layers 1+4;
    SUBv/SUBvv: layers 2-4), so those lists are shared exactly as published."""
    t7 = pd.read_excel(SRC / "supp_table7_hgea_nomenclature.xlsx", header=None)
    out, genes, after_outputs = [], {}, False
    for j in range(t7.shape[1]):
        head = str(t7.iat[1, j]).strip()
        vals = [str(x).strip() for x in t7.iloc[2:, j].dropna()]
        if "Genes" in head:
            if after_outputs:
                genes, after_outputs = {}, False
            genes[head] = vals
        elif head == "Inputs":
            label = t7.iat[0, j]
            out.append({"name": str(label).strip(), "genes": dict(genes), "inputs": vals, "outputs": [],
                        "shares_genes_with": None})
        elif head == "Outputs":
            out[-1]["outputs"] = vals
            after_outputs = True
    for k in range(1, len(out)):
        if out[k]["genes"] == out[k - 1]["genes"]:
            out[k]["shares_genes_with"] = out[k - 1]["name"]
    return out


def figure_mentions(genes: list[str]) -> dict:
    """gene -> {'figures': [...], 'text': n} from the paper and the Supplementary Information."""
    def text(pdf):
        return subprocess.run(["pdftotext", str(pdf), "-"], capture_output=True, text=True).stdout
    main, supp = text(SRC / "bienkowski2018_natneurosci.pdf"), text(SRC / "supp_text_and_figures.pdf")
    out = {g: {"figures": [], "text": 0} for g in genes}

    def legends(t, prefix, pat):
        parts = re.split(pat, t)
        # re.split with a capturing group: [pre, num, body, num, body, ...]
        for k in range(1, len(parts) - 1, 2):
            yield f"{prefix} {parts[k]}", parts[k + 1][:2500]
    for g in genes:
        w = re.compile(rf"(?<![A-Za-z0-9]){re.escape(g)}(?![A-Za-z0-9])")
        out[g]["text"] = len(w.findall(main))
        for lab, body in list(legends(main, "Fig.", r"\nFig\. (\d) \|")) + list(legends(supp, "Supplementary Fig.", r"\nSupplementary Figure (\d+)\n")):
            if w.search(body) and lab not in out[g]["figures"]:
                out[g]["figures"].append(lab)
    return out


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    t2 = pd.read_excel(SRC / "supp_table2_gene_annotation.xlsx", sheet_name="Sheet1")
    cols = [c for c in t2.columns if c != "Gene"]
    ish = json.loads((C.SITE_DATA / "ish" / "genes.json").read_text())
    ish_names = {g["gene"]: i for i, g in enumerate(ish["genes"])}

    genes = []
    for _, row in t2.iterrows():
        name = str(row.Gene).strip()
        sym = clean(name)
        alias = paren_alias(name) or ALIASES.get(sym)
        cand = [sym, alias, ALIASES.get(sym)] + [c for c in ish_names if c.lower() == sym.lower()]
        hit = next((c for c in cand if c and c in ish_names), None)
        how = None if hit is None else ("exact" if hit == name else "asterisk removed" if hit == sym and name.endswith("*")
                                        else "same symbol" if hit == sym else f"alias ({sym} -> {hit})")
        doms = [k for k, c in enumerate(cols) if pd.notna(row[c]) and str(row[c]).strip()]
        genes.append({"name": name, "symbol": sym, "asterisk": name.endswith("*"), "ish": hit, "match": how or "no coronal ISH volume",
                      "domains": doms})
    ment = figure_mentions([g["symbol"] for g in genes])
    for g in genes:
        g["figures"], g["text_mentions"] = ment[g["symbol"]]["figures"], ment[g["symbol"]]["text"]

    # Table 7 and its consistency with Table 2
    t7 = table7()
    t2sets = {c: {genes[i]["symbol"] for i in range(len(genes)) if k in genes[i]["domains"]} for k, c in enumerate(cols)}
    for d in t7:
        col = T7_TO_T2.get(d["name"], d["name"])
        allg = {clean(x) for v in d["genes"].values() for x in v}
        if col in t2sets:
            ref = t2sets[col]
            d["table2_column"] = col
            d["vs_table2"] = {"same": len(allg & ref), "only_table7": sorted(allg - ref), "only_table2": sorted(ref - allg)}
        if d["name"] in SUBREGION_LAYERS:
            d["layers"] = SUBREGION_LAYERS[d["name"]]
            per = {}
            for head, v in d["genes"].items():
                m = re.match(r"Layer (\d) Genes", head)
                if m:
                    lc = f"SUB_{m.group(1)}"
                    per[head] = {"same_as_table2": {clean(x) for x in v} == t2sets.get(lc, set())}
            d["vs_table2_layers"] = per

    # signature volumes on the hippocampus
    onto = json.loads((C.SITE_DATA / "ref" / "ontology.json").read_text())
    ann = np.frombuffer(gzip.decompress((C.SITE_DATA / "ref" / "annotation.u16.gz").read_bytes()), "<u2").reshape(C.GRID)
    keep_ids = {o["i"] for o in onto if o["acronym"] in ("HIP", "SUB", "ProS")}
    par = {o["i"]: o["parent"] for o in onto}

    def under(i):
        while i:
            if i in keep_ids:
                return True
            i = par.get(i, 0)
        return False
    lut = np.zeros(len(onto) + 1, bool)
    lut[[o["i"] for o in onto if under(o["i"])]] = True
    from scipy import ndimage
    mask = ndimage.binary_dilation(lut[ann], iterations=6)
    nz = np.argwhere(mask)
    lo, hi = nz.min(0), nz.max(0) + 1
    box = tuple(slice(a, b) for a, b in zip(lo, hi))
    mbox = mask[box]
    print(f"hippocampal crop {list(lo)} - {list(hi)} ({np.prod(hi - lo) / 1e6:.1f} M voxels)", flush=True)

    with_vol = sorted({g["ish"] for g in genes if g["ish"]})
    vol = {}

    def load(sym):
        v = _ish.decode(_ish.fetch(sym))[box].astype(np.float32)
        p = max(1.0, ish["genes"][ish_names[sym]]["p995"])
        return sym, np.clip(v / p, 0, 1).astype(np.float16)
    with cf.ThreadPoolExecutor(8) as ex:
        for k, (sym, v) in enumerate(ex.map(load, with_vol), 1):
            vol[sym] = v
            if k % 50 == 0 or k == len(with_vol):
                print(f"  ISH volumes {k}/{len(with_vol)}", flush=True)

    group_of = {c: grp for grp, cs in GROUPS for c in cs}
    sig_files = {}
    for k, c in enumerate(cols):
        markers = [g["ish"] for g in genes if k in g["domains"] and g["ish"]]
        same_group = [j for j, cc in enumerate(cols) if group_of.get(cc) == group_of.get(c) and cc != c]
        absent = sorted({g["ish"] for g in genes if g["ish"] and k not in g["domains"] and any(j in g["domains"] for j in same_group)})
        if not markers:
            continue
        s = np.mean([vol[x].astype(np.float32) for x in markers], axis=0)
        cst = np.clip(s - (np.mean([vol[x].astype(np.float32) for x in absent], axis=0) if absent else 0), 0, 1)
        s[~mbox] = 0
        cst[~mbox] = 0
        safe = re.sub(r"[^A-Za-z0-9_]+", "_", c)
        sig_files[c] = {"signature": write_array(OUT / "signatures" / f"{safe}.sig.u8.gz", np.round(s * 255).astype(np.uint8)),
                        "contrast": write_array(OUT / "signatures" / f"{safe}.con.u8.gz", np.round(cst * 255).astype(np.uint8)),
                        "markers": len(markers), "absent_reference": len(absent), "p995_sig": float(np.percentile(s[mbox], 99.5)),
                        "p995_con": float(np.percentile(cst[mbox], 99.5))}
        print(f"  {c}: {len(markers)} markers, {len(absent)} absent-reference genes", flush=True)

    # all HGEA genes together: a 50 µm stack over the hippocampus (CA, DG, SUB, ProS, FC, IG; no margin),
    # each gene = mean of its 2x2x2 normalised 25 µm voxels, x255 -> one uint8 row per gene
    core = lut[ann][box]
    sh = [(n // 2) * 2 for n in core.shape]
    core2 = core[:sh[0], :sh[1], :sh[2]].reshape(sh[0] // 2, 2, sh[1] // 2, 2, sh[2] // 2, 2).any(axis=(1, 3, 5))
    idx = np.argwhere(core2)
    stack_genes = [g["ish"] for g in genes if g["ish"]]
    rows = np.zeros((len(stack_genes), len(idx)), np.uint8)
    for k, sym in enumerate(stack_genes):
        v = vol[sym][:sh[0], :sh[1], :sh[2]].astype(np.float32).reshape(sh[0] // 2, 2, sh[1] // 2, 2, sh[2] // 2, 2).mean(axis=(1, 3, 5))
        rows[k] = np.round(v[tuple(idx.T)] * 255).astype(np.uint8)
    stack = {"step_um": 50, "origin_vox25": [int(x) for x in lo], "shape50": [int(x) // 2 for x in sh],
             "genes": stack_genes, "n_voxels": int(len(idx)),
             "voxels": write_array(OUT / "stack50.vox.u16.gz", idx.astype(np.uint16)),
             "values": write_array(OUT / "stack50.val.u8.gz", rows),
             "definition": "ISH / p99.5 (clip 1), mean over 2x2x2 voxels -> 50 µm, inside CA/DG/SUB/ProS/FC/IG, x255; "
                           "row k = genes[k] (ISH volume names), column j = voxels[j] (50 µm indices from origin)"}
    print(f"  stack: {len(stack_genes)} genes x {len(idx):,} voxels at 50 µm ({stack['values']['bytes_gz'] / 1e6:.1f} MB gz)", flush=True)

    man = {"citation": CITATION, "doi": DOI, "source_files": ["Supplementary Table 2", "Supplementary Table 7"],
           "asterisk_note": "Shown as in Supplementary Table 2; the paper and its supplement do not state what the asterisk marks.",
           "columns": cols, "groups": [{"group": g, "columns": cs} for g, cs in GROUPS], "genes": genes, "table7": t7,
           "subregion_layers": SUBREGION_LAYERS,
           "signatures": {"box_lo": [int(x) for x in lo], "box_hi": [int(x) for x in hi],
                          "definition": "signature = mean over the column's marker genes of ISH/p99.5 (clip 1); contrast = signature "
                                        "minus the mean over genes annotated in other columns of the same group but not this one (clip 0); "
                                        "x255, inside CA/DG/SUB/ProS/FC/IG + 6 voxels",
                          "files": sig_files}, "stack50": stack}
    write_json(OUT / "hgea.json", man)
    n_ish = sum(1 for g in genes if g["ish"])
    print(f"wrote {OUT}: {len(genes)} genes ({n_ish} with ISH), {len(t7)} Table 7 domains, {len(sig_files)} signature columns")
    bad = [(d["name"], d["vs_table2"]) for d in t7 if "vs_table2" in d and (d["vs_table2"]["only_table7"] or d["vs_table2"]["only_table2"])]
    print("Table 7 vs Table 2 differences:", bad if bad else "none")


if __name__ == "__main__":
    main()

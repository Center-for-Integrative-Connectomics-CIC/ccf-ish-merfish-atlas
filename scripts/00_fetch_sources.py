"""Step 0: download the original sources into data_raw/ (resumable; existing files are kept).

  python scripts/00_fetch_sources.py [--expr] [--datasets yao638850,zhuang1] [--bbp-from DIR]

- WMB taxonomy (cluster annotation terms, colours), Allen CCF 2020 parcellation tables + 10 µm annotation
- per MERFISH dataset: cell metadata with cluster annotation, gene panel, CCFv3 coordinates
  (+ the log2 expression matrix with --expr; 0.1-7.6 GB each)
- the CCFv3BBP BrainGlobe package (Zenodo 14034334) unless --bbp-from points to an unpacked copy
- the ISH gene list from the EBRAINS bucket (volumes themselves are streamed by the site / step 4)
"""

from __future__ import annotations

import argparse
import json
import sys
import zipfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from ccfatlas import config as C  # noqa: E402
from ccfatlas import sources as S  # noqa: E402

TAXONOMY_FILES = ["cluster.csv", "cluster_annotation_term.csv", "cluster_annotation_term_set.csv",
                  "cluster_to_cluster_annotation_membership.csv",
                  "views/cluster_to_cluster_annotation_membership_pivoted.csv",
                  "views/cluster_to_cluster_annotation_membership_color.csv"]
CCF2020_FILES = ["parcellation.csv", "parcellation_term.csv", "parcellation_term_set.csv",
                 "parcellation_to_parcellation_term_membership.csv",
                 "views/parcellation_to_parcellation_term_membership_acronym.csv",
                 "views/parcellation_to_parcellation_term_membership_name.csv",
                 "views/parcellation_to_parcellation_term_membership_color.csv"]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--expr", action="store_true", help="also download the log2 expression matrices")
    ap.add_argument("--datasets", default=",".join(m.key for m in C.MERFISH))
    ap.add_argument("--bbp-from", default=None, help="existing unpacked CCFv3BBP BrainGlobe folder to link")
    args = ap.parse_args()

    for f in TAXONOMY_FILES:
        print("taxonomy", f, S.abc(f"{C.TAXONOMY}/{f}").stat().st_size)
    for f in CCF2020_FILES:
        S.abc(f"{C.CCF2020}/{f}")
    S.abc("image_volumes/Allen-CCF-2020/20230630/annotation_10.nii.gz")      # the annotation the MERFISH cells used

    for key in args.datasets.split(","):
        m = C.MERFISH_BY_KEY[key]
        for rel in (m.genes, m.ccf, m.meta) + ((m.expr,) if args.expr else ()):
            p = S.abc(rel)
            print(f"{key}: {rel} ({p.stat().st_size / 1e6:.0f} MB)", flush=True)

    bbp = C.RAW / "bbp" / C.BBP_PACKAGE
    if not bbp.exists():
        if args.bbp_from:
            bbp.parent.mkdir(parents=True, exist_ok=True)
            bbp.symlink_to(Path(args.bbp_from).resolve())
        else:
            z = S.fetch("https://zenodo.org/api/records/14034334/files-archive", C.RAW / "bbp" / "zenodo_14034334.zip")
            with zipfile.ZipFile(z) as zz:
                zz.extractall(C.RAW / "bbp")
    print("CCFv3BBP:", bbp, "ok" if (bbp / "annotation.tiff").exists() else "MISSING annotation.tiff")

    genes = S.ish_gene_list()
    (C.RAW / "ish_genes.json").write_text(json.dumps(genes))
    print(f"ISH: {len(genes)} gene volumes, {sum(g['bytes'] for g in genes) / 1e9:.1f} GB in total")


if __name__ == "__main__":
    main()

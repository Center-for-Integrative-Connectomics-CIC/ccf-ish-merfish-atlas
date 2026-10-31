"""Paths, original data sources and dataset definitions.

Every dataset is read from its original publisher:
- ISH gene volumes: Carey, Piluso et al. (2026), Blue Brain Project / EBRAINS bucket "deepslice",
  4,085 genes as 25 µm NIfTI volumes in CCFv3BBP (BrainGlobe ASR grid 566 x 320 x 456).
- CCFv3BBP reference space: Piluso et al. (BrainGlobe package "ccfv3augmented_mouse_25um", same bucket /
  Zenodo 14034334).
- MERFISH: Allen Brain Cell Atlas (S3 "allen-brain-cell-atlas"), five datasets, each kept
  separate: MERFISH-C57BL6J-638850 (Yao et al. 2023) and Zhuang-ABCA-1..4 (Zhang et al. 2023),
  with the shared WMB taxonomy (Yao et al. 2023) and the Allen CCFv3 coordinates provided per cell.

Environment overrides: CCFATLAS_RAW (raw downloads), CCFATLAS_HF (Hugging Face staging folder).
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RAW = Path(os.environ.get("CCFATLAS_RAW", ROOT / "data_raw"))
WORK = ROOT / "data_work"                       # intermediate (not committed)
SITE = ROOT / "docs"                            # GitHub Pages root
SITE_DATA = SITE / "data"                       # small baked data served by Pages
HF_STAGE = Path(os.environ.get("CCFATLAS_HF", ROOT / "hf_dataset"))   # large baked data -> Hugging Face

EBRAINS_BUCKET = "https://data-proxy.ebrains.eu/api/v1/buckets/deepslice"
ISH_PREFIX = "expression_volumes/"
ISH_SUFFIX = "_interp_25um.nii.gz"
ISH_DIRECT = "https://rgw.cscs.ch/ebrains:deepslice/"          # where the proxy redirects; CORS enabled
BBP_PACKAGE = "ccfv3augmented_mouse_25um_v1.0"

ABC_S3 = "https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com"
TAXONOMY = "metadata/WMB-taxonomy/20231215"
CCF2020 = "metadata/Allen-CCF-2020/20230630"

GRID = (566, 320, 456)          # CCFv3BBP, ASR, 25 µm
VOXEL_UM = 25.0
CCFV3_AP_OFFSET_VOX = 14        # CCFv3BBP AP index = CCFv3 AP index + 14 (verified in scripts/02_align_spaces.py)


@dataclass(frozen=True)
class Merfish:
    key: str            # short id used on the site
    name: str           # ABC dataset name
    label: str
    reference: str
    expr: str           # log2 expression matrix (h5ad)
    meta: str           # cell metadata with cluster annotation
    genes: str
    ccf: str            # per-cell CCFv3 coordinates + parcellation

    def url(self, rel: str) -> str:
        return f"{ABC_S3}/{rel}"


MERFISH = [
    Merfish("yao638850", "MERFISH-C57BL6J-638850", "Allen MERFISH, C57BL/6J 638850 (Yao et al. 2023)",
            "Yao Z. et al. Nature 624, 317-332 (2023)",
            "expression_matrices/MERFISH-C57BL6J-638850/20230830/C57BL6J-638850-log2.h5ad",
            "metadata/MERFISH-C57BL6J-638850/20241115/views/cell_metadata_with_cluster_annotation.csv",
            "metadata/MERFISH-C57BL6J-638850/20241115/gene.csv",
            "metadata/MERFISH-C57BL6J-638850-CCF/20231215/ccf_coordinates.csv"),
    *[Merfish(f"zhuang{k}", f"Zhuang-ABCA-{k}", f"Zhuang lab MERFISH, ABCA-{k} (Zhang et al. 2023)",
              "Zhang M. et al. Nature 624, 343-354 (2023)",
              f"expression_matrices/Zhuang-ABCA-{k}/20230830/Zhuang-ABCA-{k}-log2.h5ad",
              f"metadata/Zhuang-ABCA-{k}/20241115/views/cell_metadata_with_cluster_annotation.csv",
              f"metadata/Zhuang-ABCA-{k}/20241115/gene.csv",
              f"metadata/Zhuang-ABCA-{k}-CCF/20230830/ccf_coordinates.csv") for k in (1, 2, 3, 4)],
]
MERFISH_BY_KEY = {m.key: m for m in MERFISH}


def is_control_probe(symbol: str) -> bool:
    """Blank control probes on the Allen 638850 panel; kept in the data, not counted or offered as genes."""
    return symbol.startswith("Blank-")

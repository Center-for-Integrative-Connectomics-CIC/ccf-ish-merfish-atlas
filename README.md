# CCF Gene Atlas · ISH × MERFISH

An interactive 3D browser for comparing mouse brain gene expression across two atlases in one common space: the Blue Brain in situ hybridisation (ISH) gene atlas and five MERFISH brains, all on the Allen Common Coordinate Framework (CCFv3BBP, 25 µm).

**Open the atlas:** https://center-for-integrative-connectomics-cic.github.io/ccf-ish-merfish-atlas/

It runs in the browser with nothing to install. Data are loaded from their original sources on demand.

## Data

| Dataset | Content | Source |
|---|---|---|
| Blue Brain ISH gene atlas | 4,085 genes as 25 µm volumes in CCFv3BBP | Carey, Piluso et al., bioRxiv 2026 (EBRAINS) |
| Allen MERFISH C57BL6J-638850 | 3.7 M cells, 500 genes | Yao et al., Nature 2023 (Allen Brain Cell Atlas) |
| Zhuang ABCA-1 to -4 | 5.4 M cells across 4 brains, 1,122 genes | Zhang et al., Nature 2023 (Allen Brain Cell Atlas) |
| Hippocampus Gene Expression Atlas (HGEA) | 258 domain-defining genes, gene × domain annotation | Bienkowski et al., Nat. Neurosci. 2018 |

Each MERFISH brain is shown separately. Only Allen 638850 samples both hemispheres; the Zhuang brains sample about one, and Zhuang-ABCA-4 has three sections. Every cell carries its WMB cell type (class → cluster), its own CCF region and the CCFv3BBP region at its position.

## Features

- **Five views of the data:** ISH, MERFISH, overlay, side by side and HGEA-ISH, each in 3D or as coronal, sagittal or horizontal slices.
- **Gene search across both atlases.** Show up to 50 genes at once, each in its own colour and toggleable; paste a list to add many together.
- **Region profiles.** For any CCF region, see the cell-type composition per brain, a gene's level in each dataset and the most enriched genes.
- **HGEA-ISH.** The genes Bienkowski et al. (2018) used to define hippocampal domains, shown on the ISH volumes exactly as published:
  - the full Supplementary Table 2 gene × domain table;
  - per-domain gene lists, signatures and Supplementary Table 7 connectivity;
  - all 256 available genes as one layer, toggleable by gene or by domain;
  - hippocampal CCF meshes.
- **Hover details:** region, coordinates, ISH value, and the nearest cell's type, labels and expression.
- **Shareable links:** the URL keeps the current view, e.g. `?gene=Prox1&brain=zhuang1&region=DG&view=coronal&mode=overlay`.

## Accuracy

The site's **Accuracy** page documents every check. In brief:

- **Space alignment.** MERFISH coordinates (CCFv3) map to CCFv3BBP with a 350 µm AP shift. Over the whole brain, region overlap is Dice 0.990, and region labels agree 97.8% at depth 5.
- **ISH vs MERFISH concordance.** Median per-gene Spearman over regions is 0.29–0.57, depending on the brain. For comparison, between two MERFISH brains it is 0.55–0.74.

## Data hosting

| Data | Hosted on |
|---|---|
| Site, reference space, cells, region tables, HGEA | GitHub Pages (`docs/`, ~220 MB) |
| ISH volumes | Streamed from the original EBRAINS bucket |
| Per-cell MERFISH expression | Hugging Face dataset [`ybb2297/ccf-ish-merfish-data`](https://huggingface.co/datasets/ybb2297/ccf-ish-merfish-data) |

## Run locally

```sh
git clone https://github.com/Center-for-Integrative-Connectomics-CIC/ccf-ish-merfish-atlas.git
cd ccf-ish-merfish-atlas
python3 scripts/serve.py          # http://127.0.0.1:8777/docs/
```

## Rebuild the data

Rebuilding needs Python ≥ 3.9 and about 30 GB of free disk. Each step reads the original sources:

```sh
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt
.venv/bin/python scripts/00_fetch_sources.py --expr     # download sources (~17 GB)
.venv/bin/python scripts/01_bake_reference.py           # CCFv3BBP ontology, annotation, template, meshes
.venv/bin/python scripts/02_align_spaces.py             # CCFv3 -> CCFv3BBP mapping and check
.venv/bin/python scripts/03_bake_merfish.py             # cells and region tables per brain
.venv/bin/python scripts/04_bake_ish.py                 # ISH region statistics (all 4,085 genes)
.venv/bin/python scripts/05_index_and_concordance.py    # dataset index and concordance
.venv/bin/python scripts/06_bake_hgea.py                # HGEA gene set (see note below)
.venv/bin/python scripts/07_publish_hf.py --repo <user>/<dataset>   # upload per-cell expression
.venv/bin/python tests/check_data.py                    # verify the baked data
```

Notes:
- **HGEA step:** step 6 needs the paper's Supplementary Information files in `data_raw/hgea/`.
- **Viewing your own build:** open the site with `?expr=local` to read expression from a locally built `hf_dataset/`, or change `HF_DATASET` in `docs/assets/js/config.js` to point at your own dataset.
- **After editing `docs/assets/`:** run `python3 scripts/bump_asset_version.py` so browsers load the new files.

## Repository layout

```
docs/          the website (GitHub Pages root; static, no build step)
  data/        baked data served with the site
src/ccfatlas/  build package: paths, sources, MERFISH readers, file formats
scripts/       build pipeline (00–07), local server, asset versioning
tests/         data checks and a headless-browser viewer test
```

Coordinates are CCFv3BBP, BrainGlobe ASR orientation, 566 × 320 × 456 voxels at 25 µm. The binary formats are described in `src/ccfatlas/bake.py`.

## Citation and licence

The code is released under the MIT licence (`LICENSE`). The data belong to their creators. If you use this atlas, please cite the original datasets; see `DATA_LICENSES.md` for sources, terms and references.

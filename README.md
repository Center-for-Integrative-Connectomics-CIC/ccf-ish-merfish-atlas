# CCF Gene Atlas · ISH × MERFISH

An interactive 3D reference for looking up mouse brain genes in two atlases at once, in one CCF space:

- **Blue Brain ISH gene atlas.** Carey, Piluso et al. 2026: 4,085 genes as 25 µm volumes on CCFv3BBP. They are streamed from the original EBRAINS bucket.
- **Five whole-brain MERFISH datasets.** From the Allen Brain Cell Atlas: Allen MERFISH C57BL6J-638850 (Yao et al. 2023) and Zhuang-ABCA-1 to -4 (Zhang et al. 2023). That is 9.1 M cells, each with its WMB cell type, its own CCF label and the CCFv3BBP label at its position. Each brain is kept separate.

**Live site:** https://center-for-integrative-connectomics-cic.github.io/ccf-ish-merfish-atlas/  
GitHub Pages serves the `docs/` folder of this repository.

## Using the site

- **Hub** (`docs/index.html`): gene finder, datasets, accuracy summary, sources.
- **3D atlas** (`docs/atlas/index.html`):
  - **Modes:** *ISH*, *MERFISH* (one brain), *overlay* (both in CCFv3BBP), *side by side* (synced views), and *HGEA-ISH* (see below).
  - **Views:** 3D (cells as points, ISH voxels above a threshold, region meshes) and coronal / sagittal / horizontal slices. The slices show the template, region outlines, the ISH heat map and the cells within a ±slab.
  - **Gene search:** covers the ISH atlas and every MERFISH panel. Cells can be coloured by gene (per brain) or by cell type at any taxonomy level.
  - **Regions:** search any CCFv3BBP region, or click a slice. The profile shows each brain's cell-type composition separately, the gene's value there in each dataset, and the most enriched genes.
  - **Hover (slices):** region, coordinates, ISH value, and the nearest cell (type at all levels, both region labels, quality, expression).
  - **Sharing:** the URL carries the state (`?gene=Prox1&brain=zhuang1&region=DG&view=coronal&mode=overlay`).
- **HGEA-ISH mode:** the Hippocampus Gene Expression Atlas gene set (Bienkowski et al. 2018) on the Blue Brain ISH volumes, exactly as published. Above the tabs:
  - **Supplementary Table 2:** 258 genes × 20 columns (19 domains / SUB layers + putative interneurons), shown as a clickable gene × domain table with 1,631 marks. Names keep their published asterisks; the paper does not define them.
  - **Gene coverage:** 256 of the 258 genes have an ISH volume, matched by exact symbol, asterisk removed or a documented alias. Epha5 and Fosl2 have none (Allen imaged them sagittally only).

  The panel has three tabs:
  - **Domains:** grouped as in the paper (DG, CA3, CA2, CA1, SUB layers, interneurons). Each domain has its gene list (genes named in the paper's figure legends are starred) and its Supplementary Table 7 entry (inputs and outputs for all 22 domains, including CA3id/ic, CA1dr/dc and the SUB subregions). Each domain also has a *signature* (mean ISH of its marker genes) and a *contrast* (that mean minus the genes the paper marks in sibling domains only). Both are derived and labelled as such.
  - **All genes:** all 256 genes as one layer at 50 µm inside the hippocampus, combined as a count (how many genes are expressed at each spot), max or mean. Genes switch on and off one by one or by domain, and hovering a slice lists the strongest genes at that spot.
  - **View:** restricts the ISH to the hippocampal formation and shows the CCF meshes of its regions, each toggleable.
- Up to four genes or derived layers can be shown together, each in its own colour.
- **Accuracy** (`docs/qc/index.html`): space alignment, cell label agreement, ISH-vs-MERFISH concordance, data handling, limits.

## Repository layout

```
docs/                      GitHub Pages root (static; three.js from a CDN, no build step)
  index.html  atlas/  qc/  assets/{css,js}
  data/                    small baked data served with the site (~220 MB)
    ref/                   CCFv3BBP ontology, 25 µm annotation + template, 887 region meshes
    merfish/<brain>/       cells (1 µm positions, cluster, both region labels, section, quality), region tables
    ish/                   ISH gene index + region means / fractions
    hgea/                  HGEA gene annotation, domain signatures, 50 µm all-genes stack
    qc/                    alignment.json, concordance.json
    datasets.json
src/ccfatlas/              build package (config, sources, merfish, bake formats)
scripts/                   build pipeline (steps 00-07), serve.py (local server), bump_asset_version.py
tests/                     smoke_viewer.py (headless-Chrome test of the viewer), check_data.py (baked-data checks)
hf_dataset/                (not committed) per-cell MERFISH expression -> Hugging Face dataset
data_raw/, data_work/      (not committed) downloads and intermediates
```

## Data hosting

| Data | Where | Size |
|---|---|---|
| Site + reference + cells + region tables + HGEA | GitHub Pages (`docs/`) | ~220 MB |
| ISH volumes | streamed from EBRAINS (`rgw.cscs.ch/ebrains:deepslice/expression_volumes/`, CORS-enabled) | 24 GB, ~6 MB per gene, loaded on demand |
| Per-cell MERFISH expression (uint8 per gene and brain) | Hugging Face dataset | ~1.8 GB, 4,388 files |

Set your Hugging Face dataset id in `docs/assets/js/config.js` (`HF_DATASET`). On `localhost` the site reads `hf_dataset/` instead. Any page also accepts `?expr=<base-url>`.

## Rebuilding

Requires Python ≥ 3.9 and about 30 GB of disk while building.

```sh
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt
.venv/bin/python scripts/00_fetch_sources.py --expr        # original sources -> data_raw/ (~17 GB)
.venv/bin/python scripts/01_bake_reference.py              # CCFv3BBP ontology, annotation, template, meshes
.venv/bin/python scripts/02_align_spaces.py                # CCFv3 -> CCFv3BBP, whole-brain check
.venv/bin/python scripts/03_bake_merfish.py                # five brains, one at a time
.venv/bin/python scripts/04_bake_ish.py --workers 8        # streams all 4,085 ISH volumes once (~40 min)
.venv/bin/python scripts/05_index_and_concordance.py       # datasets.json + concordance
.venv/bin/python scripts/06_bake_hgea.py                   # HGEA gene set (needs the paper's supplementary files in data_raw/hgea/)
.venv/bin/python scripts/07_publish_hf.py --repo <hf-user>/<dataset>   # upload hf_dataset/ (after huggingface-cli login)
```

Local preview, from the repo root so that `hf_dataset/` is reachable:

```sh
python3 scripts/serve.py        # no-cache server, then open http://127.0.0.1:8777/docs/
.venv/bin/python tests/check_data.py
.venv/bin/python tests/smoke_viewer.py
```

After editing anything in `docs/assets/`, run `python3 scripts/bump_asset_version.py` so that browsers load the new files instead of cached ones.

### Publishing on GitHub Pages
**Settings → Pages → Deploy from a branch → `main` / `/docs`.** The `.nojekyll` file in `docs/` stops Jekyll from processing the site.

## Coordinates and formats

- **Space:** CCFv3BBP, BrainGlobe ASR, 566 × 320 × 456 at 25 µm. Every layer uses µm = voxel index × 25 (voxel centres), axis order AP, DV, LR.
- **MERFISH positions:** ABC CCFv3 mm → CCFv3BBP µm as (x·1000 + 350, y·1000, z·1000). The 14-voxel AP shift is verified whole-brain in step 2.
- **Binary formats:** gzip-compressed little-endian arrays (`*.u16.gz`, `*.u8.gz`, `*.f16.gz`), with dtype and shape in the JSON manifests. Meshes are `*.mesh.gz` (uint16 µm vertices + uint32 triangles); see `src/ccfatlas/bake.py`.

## Citation and licences

Code: MIT (`LICENSE`). The data belong to their creators; see `DATA_LICENSES.md` and cite the original papers listed on the hub.

# CCF Gene Atlas · ISH × MERFISH

An interactive 3D browser for comparing mouse brain gene expression across the Blue Brain in situ hybridisation (ISH) atlas and five MERFISH brains, all in one common space (CCFv3BBP, 25 µm).

**[Open the atlas](https://center-for-integrative-connectomics-cic.github.io/ccf-ish-merfish-atlas/)**. It runs in the browser, with nothing to install.

## Features

- **Views:** ISH, MERFISH, overlay, side by side and HGEA-ISH, in 3D or as coronal, sagittal and horizontal slices.
- **Gene search:** look up a gene across both atlases and show up to 50 at once.
- **Region profiles:** see cell-type composition and gene levels for any CCF region, brain by brain.
- **HGEA-ISH:** explore the hippocampal domain genes of Bienkowski et al. (2018) on the ISH volumes.
- **Shareable links:** the URL keeps the current view.

## Data

| Dataset | Content | Reference |
|---|---|---|
| Blue Brain ISH atlas | 4,085 genes, 25 µm volumes | Carey, Piluso et al., bioRxiv 2026 |
| Allen MERFISH C57BL6J-638850 | 3.7 M cells, 500 genes | Yao et al., Nature 2023 |
| Zhuang ABCA-1 to -4 | 5.4 M cells across 4 brains, 1,122 genes | Zhang et al., Nature 2023 |
| Hippocampus Gene Expression Atlas | 258 genes, gene × domain annotation | Bienkowski et al., Nat. Neurosci. 2018 |

Only the Allen brain samples both hemispheres; the Zhuang brains sample about one. The site's **Accuracy** page covers alignment, ISH and MERFISH agreement, and other limits to keep in mind.

## Run locally

```sh
git clone https://github.com/Center-for-Integrative-Connectomics-CIC/ccf-ish-merfish-atlas.git
cd ccf-ish-merfish-atlas
python3 scripts/serve.py    # then open http://127.0.0.1:8777/docs/
```

To rebuild the data from the original sources, install `requirements.txt` and run `scripts/00` to `07` in order.

## Citation and licence

The code is under the MIT licence (see `LICENSE`). The data belong to their creators, so please cite the original datasets; see `DATA_LICENSES.md` for sources and terms.

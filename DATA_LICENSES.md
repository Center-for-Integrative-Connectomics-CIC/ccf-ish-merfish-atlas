# Data sources and terms

This repository contains code that downloads and re-packages third-party data for viewing. The data stay
the property of their creators. **Check each source's current terms before redistributing**; the notes
below are pointers, not legal advice.

| Data | Source | Terms (check the source) | Cite |
|---|---|---|---|
| ISH gene volumes (4,085 genes, 25 µm, CCFv3BBP) | EBRAINS data-proxy bucket `deepslice`, `expression_volumes/` | as stated by the EBRAINS dataset / bucket (published as CC BY 4.0 at the time of writing) | Carey H., Piluso S. et al., bioRxiv 2026, doi:10.64898/2026.01.20.700446 |
| CCFv3BBP reference space (annotation, template, meshes, ontology) | Zenodo record 14034334; BrainGlobe package `ccfv3augmented_mouse_25um` | see the Zenodo record | Piluso S. et al., bioRxiv 2024, doi:10.1101/2024.11.06.622212 |
| MERFISH C57BL6J-638850, WMB taxonomy, CCF registration | Allen Brain Cell Atlas, `s3://allen-brain-cell-atlas` | Allen Institute terms of use / the ABC atlas data licence | Yao Z. et al., Nature 624, 317–332 (2023) |
| MERFISH Zhuang-ABCA-1 to -4 | Allen Brain Cell Atlas, `s3://allen-brain-cell-atlas` | as above | Zhang M. et al., Nature 624, 343–354 (2023) |
| HGEA gene annotation (Supplementary Tables 2 and 7) | Supplementary Information of the article (Nature Neuroscience) | publisher's terms; the gene x domain annotation is reproduced for reference with citation | Bienkowski M.S. et al., Nat. Neurosci. 21, 1628–1643 (2018), doi:10.1038/s41593-018-0241-y |
| Allen CCFv3 / CCF-2020 parcellation | Allen Brain Cell Atlas `image_volumes/Allen-CCF-2020`, `metadata/Allen-CCF-2020` | Allen Institute terms of use | Wang Q. et al., Cell 181, 936–953 (2020) |

What this repository derives and ships:
- `docs/data/hgea/`: the HGEA gene annotation tables (as published), gene-to-volume matching, domain signature volumes and the 50 µm all-genes stack, all derived from the ISH.
- `docs/data/`: region meshes (decimated), the annotation as indices, the template as 8-bit, MERFISH cell positions, cell types, region labels and quality scores, and region-level statistics.
- The Hugging Face dataset: per-cell MERFISH expression, quantised to 8 bits.

The ISH volumes are not redistributed; the site streams them from the original bucket.

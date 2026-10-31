"""Step 7: publish the per-cell MERFISH expression (hf_dataset/) as a public Hugging Face dataset.

  huggingface-cli login                      # once, with a write token
  python scripts/07_publish_hf.py --repo <hf-user>/ccf-ish-merfish-data [--dry-run]

Writes a dataset card (README.md) with the file layout and encoding, then uploads hf_dataset/ with
resumable large-folder upload. Afterwards set HF_DATASET in docs/assets/js/config.js to the same repo id.
The site reads files as  https://huggingface.co/datasets/<repo>/resolve/main/merfish/<brain>/expr/<gene>.u8.gz
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from ccfatlas import config as C  # noqa: E402

CARD = """---
license: other
pretty_name: CCF Gene Atlas - per-cell MERFISH expression
tags: [neuroscience, mouse-brain, spatial-transcriptomics, merfish, allen-ccf]
---
# CCF Gene Atlas · per-cell MERFISH expression

Per-cell log2 expression of every panel gene for five MERFISH brains from the Allen Brain Cell
Atlas, packaged for the browser viewer of the CCF Gene Atlas (ISH x MERFISH). Each brain is kept separate.

| brain | dataset | cells | genes |
|---|---|---|---|
{rows}

Layout: `merfish/<brain>/expr/<gene>.u8.gz`, gzip-compressed uint8 arrays, one value per cell. Cell order and
per-gene `vmax` are in the site's `docs/data/merfish/<brain>/manifest.json`. Decoding: log2 value = q / 255 * vmax,
and 0 means not detected.

Sources and citation:
- Yao Z. et al., Nature 624, 317-332 (2023): MERFISH C57BL6J-638850.
- Zhang M. et al., Nature 624, 343-354 (2023): Zhuang-ABCA-1 to -4.
- Allen Brain Cell Atlas (Allen Institute for Brain Science).

This is a re-packaging for visualisation; the original data and their terms are at the Allen Brain Cell Atlas.
"""


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", required=True)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    rows = []
    n_files = size = 0
    for m in C.MERFISH:
        man = json.loads((C.SITE_DATA / "merfish" / m.key / "manifest.json").read_text())
        files = list((C.HF_STAGE / "merfish" / m.key / "expr").glob("*.u8.gz"))
        n_files += len(files)
        size += sum(f.stat().st_size for f in files)
        rows.append(f"| {m.key} | {m.name} | {man['cells']:,} | {len(files)} |")
    (C.HF_STAGE / "README.md").write_text(CARD.format(rows="\n".join(rows)))
    print(f"{n_files} files, {size / 1e9:.2f} GB in {C.HF_STAGE}")
    if args.dry_run:
        return
    from huggingface_hub import HfApi
    api = HfApi()
    api.create_repo(args.repo, repo_type="dataset", exist_ok=True)
    api.upload_large_folder(repo_id=args.repo, repo_type="dataset", folder_path=str(C.HF_STAGE))
    print(f"uploaded -> https://huggingface.co/datasets/{args.repo}")


if __name__ == "__main__":
    main()

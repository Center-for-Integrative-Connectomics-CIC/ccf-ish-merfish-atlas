"""Step 4: index and summarise the Blue Brain ISH gene volumes (Carey, Piluso et al. 2026).

Streams every 25 µm volume once from the EBRAINS bucket (CORS-enabled; the site streams the same files
on demand), decodes it onto the CCFv3BBP grid and computes, per gene:
  - per structure (own voxels + descendants): mean ISH value and fraction of voxels > 0
  - display range: max and 99.5th percentile inside the brain
Volumes are not kept (24 GB). Orientation: NIfTI (LR, AP, DV) -> ASR via transpose(1, 2, 0)[::-1, ::-1, :]
(checked: Prox1 peaks in the dentate gyrus granule layer, see the printed check).
Output: docs/data/ish/genes.json, region_mean.f16.gz, region_frac.f16.gz (regions x genes)
  python scripts/04_bake_ish.py [--workers 8] [--limit N]
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import gzip
import io
import json
import sys
import time
import urllib.request
from pathlib import Path

import nibabel as nib
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from ccfatlas import config as C  # noqa: E402
from ccfatlas import sources as S  # noqa: E402
from ccfatlas.bake import write_array, write_json  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import importlib.util  # noqa: E402

_spec = importlib.util.spec_from_file_location("m3", Path(__file__).resolve().parent / "03_bake_merfish.py")
_m3 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_m3)


def decode(blob: bytes) -> np.ndarray:
    img = nib.Nifti1Image.from_bytes(gzip.decompress(blob))
    vol = np.asarray(img.dataobj).transpose(1, 2, 0)[::-1, ::-1, :]
    assert vol.shape == C.GRID, vol.shape
    return np.ascontiguousarray(vol)


def fetch(gene: str, tries: int = 4) -> bytes:
    for k in range(tries):
        try:
            with urllib.request.urlopen(S.ish_url(gene), timeout=180) as r:
                return r.read()
        except Exception:
            if k == tries - 1:
                raise
            time.sleep(2 + 3 * k)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()
    onto = json.loads((C.SITE_DATA / "ref" / "ontology.json").read_text())
    ann = np.frombuffer(gzip.decompress((C.SITE_DATA / "ref" / "annotation.u16.gz").read_bytes()), "<u2").reshape(C.GRID)
    flat = ann.ravel().astype(np.int64)
    inside = flat > 0
    n_reg = len(onto) + 1
    A = _m3.ancestors_matrix(onto)
    nvox = A @ np.bincount(flat, minlength=n_reg).astype(np.float64)
    genes = json.loads((C.RAW / "ish_genes.json").read_text())
    if args.limit:
        genes = genes[: args.limit]
    G = len(genes)
    mean = np.full((n_reg, G), np.nan, np.float32)
    frac = np.full((n_reg, G), np.nan, np.float32)
    info = [None] * G

    def work(j):
        g = genes[j]["gene"]
        v = decode(fetch(g)).ravel()
        s = np.bincount(flat, weights=v, minlength=n_reg)
        p = np.bincount(flat, weights=(v > 0), minlength=n_reg)
        vin = v[inside]
        return j, A @ s, A @ p, {"gene": g, "bytes_gz": genes[j]["bytes"], "dtype": str(v.dtype),
                                 "max": float(vin.max()), "p995": float(np.percentile(vin, 99.5)),
                                 "frac_pos": float((vin > 0).mean())}

    t0 = time.time()
    done = 0
    with cf.ThreadPoolExecutor(args.workers) as ex:
        for j, s, p, meta in ex.map(work, range(G)):
            ok = nvox > 0
            mean[ok, j] = s[ok] / nvox[ok]
            frac[ok, j] = p[ok] / nvox[ok]
            info[j] = meta
            done += 1
            if done % 100 == 0 or done == G:
                print(f"  {done}/{G} genes [{time.time() - t0:.0f}s]", flush=True)

    names = {o["i"]: o["acronym"] for o in onto}
    if any(x["gene"] == "Prox1" for x in info):
        j = [x["gene"] for x in info].index("Prox1")
        top = np.argsort(-np.nan_to_num(mean[:, j], nan=-1))[:5]
        print("orientation check, Prox1 top structures:", [names.get(i) for i in top])
    out = C.SITE_DATA / "ish"
    man = {"source": "Carey, Piluso et al. (2026), Blue Brain Project; EBRAINS bucket 'deepslice'",
           "url_template": C.ISH_DIRECT + C.ISH_PREFIX + "{gene}" + C.ISH_SUFFIX,
           "proxy_template": C.EBRAINS_BUCKET + "/" + C.ISH_PREFIX + "{gene}" + C.ISH_SUFFIX,
           "grid": list(C.GRID), "orientation": "NIfTI (LR, AP, DV) -> ASR: transpose(1,2,0)[::-1, ::-1, :]",
           "genes": info, "region_voxels": nvox.astype(int).tolist()}
    man["files"] = {"region_mean": write_array(out / "region_mean.f16.gz", mean.astype(np.float16)),
                    "region_frac": write_array(out / "region_frac.f16.gz", frac.astype(np.float16))}
    write_json(out / "genes.json", man)
    print(f"wrote {out} [{time.time() - t0:.0f}s]")


if __name__ == "__main__":
    main()

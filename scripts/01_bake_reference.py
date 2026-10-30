"""Step 1: bake the CCFv3BBP reference space for the site (docs/data/ref/).

Outputs
  ontology.json        every CCFv3BBP structure: index, id, acronym, name, parent index, colour, depth,
                       voxel count (own + descendants), whether a mesh exists
  annotation.u16.gz    25 µm annotation as uint16 structure indices (566 x 320 x 456, ASR)
  reference.u8.gz      the CCFv3BBP average template, 1st-99.5th percentile -> uint8
  meshes/<id>.mesh.gz  region surfaces, quadric-decimated (pyfqmr), µm in the same frame as the volumes
Coordinates: µm = voxel index x 25 (voxel centres), axis order AP, DV, LR -- shared by every layer.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pyfqmr
import tifffile

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from ccfatlas import config as C  # noqa: E402
from ccfatlas.bake import write_array, write_json, write_mesh  # noqa: E402

OUT = C.SITE_DATA / "ref"


def load_obj(path: Path):
    v, f = [], []
    with open(path) as fh:
        for line in fh:
            if line.startswith("v "):
                v.append(line.split()[1:4])
            elif line.startswith("f "):
                f.append([int(x.split("/")[0]) - 1 for x in line.split()[1:4]])
    return np.array(v, float), np.array(f, np.int64)


def target_faces(n_vox: int, is_root: bool) -> int:
    if is_root:
        return 40000
    return int(np.clip(300 + np.sqrt(n_vox) * 6, 600, 6000))


def main() -> None:
    bbp = C.RAW / "bbp" / C.BBP_PACKAGE
    ann = tifffile.imread(bbp / "annotation.tiff")
    ref = tifffile.imread(bbp / "reference.tiff")
    structs = json.loads((bbp / "structures.json").read_text())
    assert ann.shape == C.GRID

    ids = [s["id"] for s in structs]
    index = {sid: i + 1 for i, sid in enumerate(ids)}           # 0 = outside the brain
    lut = np.zeros(int(ann.max()) + 1, np.uint16)
    for sid, i in index.items():
        if sid <= ann.max():
            lut[sid] = i
    unknown = set(np.unique(ann)) - set(index) - {0}
    assert not unknown, f"annotation ids missing from the ontology: {sorted(unknown)[:10]}"
    ann_idx = lut[ann]

    own = np.bincount(ann_idx.ravel(), minlength=len(ids) + 1)
    total = own.copy()
    parent = {}
    for s in structs:
        p = s["structure_id_path"]
        parent[s["id"]] = p[-2] if len(p) > 1 else None
        for anc in p[:-1]:
            total[index[anc]] += own[index[s["id"]]]

    onto = []
    mesh_dir = bbp / "meshes"
    (OUT / "meshes").mkdir(parents=True, exist_ok=True)
    mesh_bytes = 0
    for s in structs:
        i = index[s["id"]]
        has_mesh = (mesh_dir / f"{s['id']}.obj").exists() and total[i] > 0
        if has_mesh:
            v, f = load_obj(mesh_dir / f"{s['id']}.obj")
            tgt = target_faces(int(total[i]), s["id"] == 997)
            if len(f) > tgt:
                simp = pyfqmr.Simplify()
                simp.setMesh(np.ascontiguousarray(v, np.float32), np.ascontiguousarray(f, np.int32))   # needs f32/i32
                simp.simplify_mesh(target_count=tgt, aggressiveness=7, preserve_border=False, max_iterations=300, verbose=0)
                v, f, _ = simp.getMesh()
            mesh_bytes += write_mesh(OUT / "meshes" / f"{s['id']}.mesh.gz", v, f)
        onto.append({"i": i, "id": s["id"], "acronym": s["acronym"], "name": s["name"],
                     "parent": index.get(parent[s["id"]], 0), "rgb": s["rgb_triplet"],
                     "depth": len(s["structure_id_path"]) - 1, "own_vox": int(own[i]), "vox": int(total[i]),
                     "mesh": bool(has_mesh)})
    write_json(OUT / "ontology.json", onto)
    print(f"ontology: {len(onto)} structures, {sum(o['mesh'] for o in onto)} meshes ({mesh_bytes / 1e6:.1f} MB)")

    man = {"grid": list(C.GRID), "voxel_um": C.VOXEL_UM, "axes": ["AP", "DV", "LR"], "orientation": "ASR",
           "space": "CCFv3BBP (Piluso et al.)", "coords": "um = voxel index x 25 (voxel centres)"}
    man["annotation"] = write_array(OUT / "annotation.u16.gz", ann_idx)
    inside = ann > 0
    lo, hi = np.percentile(ref[inside], [1, 99.5])
    ref8 = np.clip((ref.astype(np.float32) - lo) / (hi - lo) * 255, 0, 255).astype(np.uint8)
    ref8[~inside] = 0
    man["reference"] = write_array(OUT / "reference.u8.gz", ref8)
    write_json(OUT / "manifest.json", man, compact=False)
    print(json.dumps(man, indent=1))


if __name__ == "__main__":
    main()

"""Binary formats shared by the build scripts and the web viewer (docs/assets/js/io.js).

Array files are raw little-endian C-order arrays, gzip-compressed (`.gz`); the browser inflates them with
DecompressionStream('gzip'). Dtype, shape and any value scaling live in the JSON manifests next to them.

Mesh files (`.mesh.gz`): uint32 n_vertices, uint32 n_triangles, uint16 vertex coords (µm, AP DV LR, 3 x n),
padding to 4 bytes, uint32 triangle indices (3 x n).
"""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import numpy as np


def write_array(path: Path, a: np.ndarray, level: int = 6) -> dict:
    path.parent.mkdir(parents=True, exist_ok=True)
    a = np.ascontiguousarray(a)
    raw = a.astype(a.dtype.newbyteorder("<"), copy=False).tobytes()
    path.write_bytes(gzip.compress(raw, compresslevel=level, mtime=0))
    return {"file": path.name, "dtype": a.dtype.name, "shape": list(a.shape), "bytes_gz": path.stat().st_size}


def write_json(path: Path, obj, compact: bool = True) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, separators=(",", ":")) if compact else json.dumps(obj, indent=1))


def write_mesh(path: Path, verts_um: np.ndarray, tris: np.ndarray) -> int:
    v = np.clip(np.round(verts_um), 0, 65535).astype("<u2")
    t = tris.astype("<u4")
    head = np.array([len(v), len(t)], "<u4").tobytes()
    body = v.tobytes()
    body += b"\0" * ((-len(body)) % 4)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(gzip.compress(head + body + t.tobytes(), compresslevel=9, mtime=0))
    return path.stat().st_size


def quantize_u8(x: np.ndarray, vmax: float) -> np.ndarray:
    """0 stays 0; (0, vmax] -> 1..255 (so any detected expression is visible)."""
    q = np.zeros(x.shape, np.uint8)
    pos = x > 0
    q[pos] = np.clip(np.ceil(x[pos] / vmax * 255.0), 1, 255).astype(np.uint8)
    return q

// Loaders for the baked formats (see src/ccfatlas/bake.py) and the original ISH NIfTI volumes.
import { CONFIG } from './config.js?v=20261029210827';

export const GRID = [566, 320, 456];          // CCFv3BBP, ASR (AP, DV, LR), 25 µm
export const VOX = 25;
const DT = { uint8: Uint8Array, uint16: Uint16Array, uint32: Uint32Array, float32: Float32Array, int16: Int16Array };

async function readAll(resp, onProgress) {
  const total = +resp.headers.get('content-length') || 0;
  const reader = resp.body.getReader();
  const chunks = []; let got = 0;
  for (;;) {
    const { done, value } = await reader.read();
    if (done) break;
    chunks.push(value); got += value.length;
    if (onProgress) onProgress(got, total);
  }
  const out = new Uint8Array(got); let o = 0;
  for (const c of chunks) { out.set(c, o); o += c.length; }
  return out;
}

export async function inflate(bytes) {
  // bytes may already be inflated if a server applied Content-Encoding: gzip
  if (!(bytes[0] === 0x1f && bytes[1] === 0x8b)) return bytes.buffer;
  const ds = new Blob([bytes]).stream().pipeThrough(new DecompressionStream('gzip'));
  return await new Response(ds).arrayBuffer();
}

export async function fetchGz(url, onProgress) {
  const r = await fetch(url);
  if (!r.ok) throw new Error(`${r.status} ${url}`);
  return inflate(await readAll(r, onProgress));
}

export async function json(url) {
  const r = await fetch(url);
  if (!r.ok) throw new Error(`${r.status} ${url}`);
  return r.json();
}

export async function array(url, dtype, onProgress) {
  const buf = await fetchGz(url, onProgress);
  if (dtype === 'float16') return f16(new Uint16Array(buf));
  return new DT[dtype](buf);
}

export function f16(u16) {
  const out = new Float32Array(u16.length);
  for (let i = 0; i < u16.length; i++) {
    const h = u16[i], s = h & 0x8000 ? -1 : 1, e = (h >> 10) & 0x1f, m = h & 0x3ff;
    out[i] = e === 0 ? s * m * 5.960464477539063e-8 : e === 31 ? (m ? NaN : s * Infinity) : s * (1 + m / 1024) * Math.pow(2, e - 15);
  }
  return out;
}

export const data = (path) => CONFIG.DATA + path;

export async function mesh(id) {
  const buf = await fetchGz(data(`ref/meshes/${id}.mesh.gz`));
  const head = new Uint32Array(buf, 0, 2), nv = head[0], nt = head[1];
  const v = new Uint16Array(buf, 8, nv * 3);
  const off = 8 + Math.ceil(nv * 6 / 4) * 4;
  return { um: v, index: new Uint32Array(buf, off, nt * 3) };
}

// ISH: original Blue Brain NIfTI (uint8, dims LR x AP x DV, Fortran order) -> ASR uint8 array.
// asr[(ap * 320 + dv) * 456 + lr] = nifti[lr + 456 * ((566 - 1 - ap) + 566 * (320 - 1 - dv))]  (checked against nibabel)
export async function ishVolume(urlTemplate, gene, onProgress) {
  const raw = new Uint8Array(await fetchGz(urlTemplate.replace('{gene}', encodeURIComponent(gene)), onProgress));
  const dv = new DataView(raw.buffer, raw.byteOffset);
  const nx = dv.getInt16(42, true), ny = dv.getInt16(44, true), nz = dv.getInt16(46, true);
  const dtype = dv.getInt16(70, true), off = Math.round(dv.getFloat32(108, true));
  if (nx !== 456 || ny !== 566 || nz !== 320 || dtype !== 2) throw new Error(`unexpected ISH volume layout ${nx}x${ny}x${nz} type ${dtype}`);
  const src = raw.subarray(off);
  const out = new Uint8Array(566 * 320 * 456);
  for (let ap = 0; ap < 566; ap++) {
    for (let d = 0; d < 320; d++) {
      const so = ((565 - ap) + 566 * (319 - d)) * 456, oo = (ap * 320 + d) * 456;
      out.set(src.subarray(so, so + 456), oo);
    }
  }
  return out;
}

// display frame: scene units = mm, centred on the grid; x = LR, y = up (dorsal), z = anterior
export const CENTER = [GRID[0] * VOX / 2, GRID[1] * VOX / 2, GRID[2] * VOX / 2];
export function toScene(ap, dvv, lr) {
  return [(lr - CENTER[2]) / 1000, -(dvv - CENTER[1]) / 1000, -(ap - CENTER[0]) / 1000];
}

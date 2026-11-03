// HGEA-ARA tab: the HGEA domains on the 3D Allen Reference Atlas, with this site's MERFISH cells mapped into it.
import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import * as IO from './io.js?v=20261102182233';
import { CONFIG } from './config.js?v=20261102182233';
import { lut } from './colormaps.js?v=20261102182233';

const $ = (id) => document.getElementById(id);
const EDITION = 'ara';
const DATA = new URL('../../data/hgea_ara/', import.meta.url).href;
const SITE_DATA = new URL('../../data/', import.meta.url).href;
if (new URLSearchParams(location.search).has('embed')) document.body.classList.add('embed');   // inside the atlas page's HGEA-ARA tab

// ───────── data ─────────
async function gz(name) {
  const r = await fetch(DATA + name);
  if (!r.ok) throw new Error(`${name}: HTTP ${r.status}`);
  return new Response(r.body.pipeThrough(new DecompressionStream('gzip'))).arrayBuffer();
}
const json = (name) => fetch(DATA + name).then((r) => r.json());

let META, ONTO, NISSL, ARA, LAB = {}, NA, ND, NL;
const S = { view: 'coronal', layer: 'full', nissl: true, outlines: false, opacity: 0.7, on: new Set(),
  pos: { coronal: 313, sagittal: 140, horizontal: 160 }, zoom: 1, tx: 0, ty: 0, nissl3d: 'volume', density: 0.5, shell: false, svg: false, svgAlpha: 0.8 };
let DOMRGB = new Uint8Array(256 * 3), DOMNAME = {};

async function boot() {
  [META, ONTO] = await Promise.all([json('meta.json'), json('ontology.json')]);
  [NA, ND, NL] = META.grid.shape;
  const [n, a, f, p] = await Promise.all(['nissl', 'ara', 'hgea_full', 'hgea_plates'].map((k) => gz(META.files[k].file)));
  NISSL = new Uint8Array(n); ARA = new Uint16Array(a); LAB.full = new Uint8Array(f); LAB.plates = new Uint8Array(p);
  for (const d of META.domains) {
    const rgb = [1, 3, 5].map((k) => parseInt(d.colour.slice(k, k + 2), 16));
    for (const id of [d.id, d.id + 100]) DOMRGB.set(rgb, id * 3);
    DOMNAME[d.id] = d.name; DOMNAME[d.id + 100] = `${d.name} (left)`;
    S.on.add(d.id);
  }
  S.pos.coronal = META.plates[Math.floor(META.plates.length / 2)].z;
  const hp = LAB.full; let lrs = 0, lrn = 0, dvs = 0;
  for (let i = 0; i < hp.length; i += 97) if (hp[i]) { const r = i % (ND * NL); lrs += Math.min(r % NL, NL - 1 - (r % NL)); dvs += Math.floor(r / NL); lrn++; }
  if (lrn) { S.pos.sagittal = Math.round(lrs / lrn); S.pos.horizontal = Math.round(dvs / lrn); }
  S.opacity = $('opacity').value / 100; ui(); init3D(); setView('coronal');
  await initCells();
  if (META.svg !== false && $('svg').checked) await $('svg').onchange();   // the SVG plates are on by default
  $('loading').style.opacity = 0;
}

// ───────── panel ─────────
function ui() {
  $('domains').innerHTML = META.domains.map((d) => `<div class="dom" data-id="${d.id}"><span class="sw" style="background:${d.colour}"></span>${d.name}</div>`).join('');
  $('domains').querySelectorAll('.dom').forEach((el) => el.addEventListener('click', () => {
    const id = +el.dataset.id; S.on.has(id) ? S.on.delete(id) : S.on.add(id); el.classList.toggle('off', !S.on.has(id)); changed();
  }));
  $('allon').onclick = () => { META.domains.forEach((d) => S.on.add(d.id)); syncDoms(); changed(); };
  $('alloff').onclick = () => { S.on.clear(); syncDoms(); changed(); };
  seg('layer', (v) => { S.layer = v; changed(); });
  seg('view', setView);
  $('nissl').onchange = () => { S.nissl = $('nissl').checked; redraw(); };
  $('outlines').onchange = () => { S.outlines = $('outlines').checked; redraw(); };
  $('opacity').oninput = () => { S.opacity = $('opacity').value / 100; if (hgPts) { hgPts.material.opacity = S.opacity; hgPts.material.transparent = S.opacity < 1; need3D = true; } redraw(); };
  seg('nissl3d', (v) => { S.nissl3d = v; update3DNissl(); });
  $('svg').onchange = async () => {
    S.svg = $('svg').checked; $('svgrow').style.display = S.svg ? 'flex' : 'none';
    if (S.svg) await loadSvgs();
    updateSvg3D(); redraw();
  };
  $('svgop').oninput = () => { S.svgAlpha = $('svgop').value / 100; if (svgPlanes) svgPlanes.forEach((m) => { m.material.opacity = S.svgAlpha; }); need3D = true; redraw(); };
  $('density').oninput = () => { S.density = $('density').value / 100; update3DNissl(); };
  $('shell').onchange = () => { S.shell = $('shell').checked; if (shell) shell.visible = S.shell; need3D = true; };
  $('pos').oninput = () => { S.pos[pl()] = +$('pos').value; planesDirty = true; if (SEL) closeRoi(false); redraw(); };
  $('prev').onclick = () => jump(-1); $('next').onclick = () => jump(1);
  addEventListener('keydown', (e) => { if (e.target.tagName === 'INPUT') return; if (e.key === ']') jump(1); if (e.key === '[') jump(-1); });
  if (META.accuracy_rows) {   // CCFv3BBP edition: rows baked with the data
    $('acc').innerHTML = '<tr><th>measure</th><th class="n">value</th></tr>' + META.accuracy_rows.map(([m, v]) => `<tr><td>${m}</td><td class="n">${v}</td></tr>`).join('');
  } else {
    const A = META.accuracy, L = A.lopo_label_accuracy, D = A.lopo_hgea_dice;
    $('acc').innerHTML = `<tr><th>measure</th><th class="n">value</th></tr>
      <tr><td>placement, ARA lines (median)</td><td class="n">${A.placement_line_median_um} µm</td></tr>
      <tr><td>placement, hippocampus</td><td class="n">${A.placement_hpf_median_um} µm</td></tr>
      <tr><td>fill, ${A.fill_guide}-guided · label accuracy</td><td class="n">${L[A.fill_guide]}</td></tr>
      <tr><td>fill · HGEA overlap (Dice)</td><td class="n">${D[A.fill_guide]}</td></tr>
      <tr><td>plain interpolation, for comparison</td><td class="n">${L.none}</td></tr>`;
  }
  if (META.accuracy_note) $('accnote').textContent = META.accuracy_note;
  if (META.method) $('method').textContent = META.method;
  if (META.title) { const [a, b] = META.title.split(' · '); $('brand').innerHTML = `${a} · <em>${b}</em>`; document.title = META.title; }
  if (META.subtitle) $('subtitle').textContent = META.subtitle;
  if (META.svg === false) $('svgwrap').style.display = 'none';   // the SVG plates are drawn on the 3D ARA slices only
  if (!$('edition')) { $('sources').innerHTML = META.sources.map((x) => `· ${x}`).join('<br>'); return; }
  document.querySelectorAll('#edition button').forEach((b) => {
    b.classList.toggle('active', b.dataset.v === EDITION);
    b.onclick = () => { if (b.dataset.v !== EDITION) location.search = b.dataset.v === 'ara' ? '' : `?edition=${b.dataset.v}`; };
  });
  $('sources').innerHTML = META.sources.map((s) => `· ${s}`).join('<br>');
}
function syncDoms() { $('domains').querySelectorAll('.dom').forEach((el) => el.classList.toggle('off', !S.on.has(+el.dataset.id))); }
function seg(id, fn) {
  document.querySelectorAll(`#${id} button`).forEach((b) => b.addEventListener('click', () => {
    document.querySelectorAll(`#${id} button`).forEach((x) => x.classList.toggle('active', x === b)); fn(b.dataset.v);
  }));
}
function changed() { redraw(); if (S.view === '3d') build3D(); else dirty3D = true; }

// ───────── slices ─────────
const pl = () => (S.view === 'cut3d' ? 'coronal' : S.view);   // slice plane behind the current view
const is3D = () => S.view === '3d' || S.view === 'cut3d';
const PLANES = { coronal: { n: () => NA, W: () => NL, H: () => ND }, sagittal: { n: () => NL, W: () => NA, H: () => ND },
  horizontal: { n: () => ND, W: () => NL, H: () => NA } };
function vox(plane, k, u, v) {   // -> grid index
  if (plane === 'coronal') return (k * ND + v) * NL + u;
  if (plane === 'sagittal') return (u * ND + v) * NL + k;
  return (v * ND + k) * NL + u;
}
function setView(v) {
  S.view = v; if (SEL) closeRoi(false);
  const is3 = is3D();
  $('gl').style.display = is3 ? 'block' : 'none'; $('slice').style.display = is3 ? 'none' : 'block';
  $('slicebar').style.display = v === '3d' ? 'none' : 'flex';
  $('prev').style.display = $('next').style.display = pl() === 'coronal' ? '' : 'none';
  $('hint').textContent = is3 ? 'drag rotate · scroll zoom · right-drag pan' + (v === 'cut3d' ? ' · [ ] plates' : '') : 'scroll zoom · drag pan · [ ] plates';
  if (v !== '3d') { $('pos').max = PLANES[pl()].n() - 1; $('pos').value = S.pos[pl()]; }
  if (!is3) { S.zoom = 1; S.tx = S.ty = 0; }
  else { if (dirty3D || !hgPts) build3D(); if (v === 'cut3d') aimCut(); update3DNissl(); }
  $('sec3d').style.display = is3 ? 'block' : 'none';
  resize(); redraw();
}
function jump(step) {
  if (pl() !== 'coronal') return;
  const z = S.pos.coronal, zs = META.plates.map((p) => p.z);
  const t = step > 0 ? zs.find((p) => p > z) : [...zs].reverse().find((p) => p < z);
  if (t != null) { S.pos.coronal = t; $('pos').value = t; if (SEL) closeRoi(false); redraw(); }
}
function sliceImage(plane = pl(), k = S.pos[pl()], tissueAlpha = 255) {
  const P = PLANES[plane], W = P.W(), H = P.H(), lab = LAB[S.layer];
  const img = new ImageData(W, H), px = img.data, fo = S.focus && G.gene, al = fo ? S.opacity * 0.2 : S.opacity;   // fo: a gene picked in a region popup
  const ish = G.ishAra && $('ishon').checked ? G.ishAra : null, IL = ish ? lut('magma') : null, lo = (fo ? 0.12 : 0.2) * G.p995, span = Math.max(1, G.p995 - lo);
  for (let v = 0; v < H; v++) for (let u = 0; u < W; u++) {
    const i = vox(plane, k, u, v), o = (v * W + u) * 4, a = ARA[i];
    let r = 0, g = 0, b = 0, A = a ? tissueAlpha : 0;
    if (a) { const t = S.nissl ? 25 + NISSL[i] * (fo ? 0.45 : 0.85) : 46; r = g = b = t; }
    if (ish && a && ish[i] > lo) {   // Blue Brain ISH heat under the HGEA
      const f = Math.min(1, (ish[i] - lo) / span), q = Math.round(f * 255) * 3, w = fo ? 0.55 + 0.45 * f : 0.3 + 0.6 * f;
      r = r * (1 - w) + IL[q] * w; g = g * (1 - w) + IL[q + 1] * w; b = b * (1 - w) + IL[q + 2] * w;
    }
    const h = lab[i];
    const edge = S.hgeaLines && h && ((u + 1 < W && lab[vox(plane, k, u + 1, v)] !== h) || (v + 1 < H && lab[vox(plane, k, u, v + 1)] !== h)
      || (u > 0 && lab[vox(plane, k, u - 1, v)] !== h) || (v > 0 && lab[vox(plane, k, u, v - 1)] !== h));
    if (h && S.on.has(h % 100) && (!S.hgeaLines || edge)) { r = r * (1 - al) + DOMRGB[h * 3] * al; g = g * (1 - al) + DOMRGB[h * 3 + 1] * al; b = b * (1 - al) + DOMRGB[h * 3 + 2] * al; A = 255; }
    if (S.outlines && a && ((u + 1 < W && ARA[vox(plane, k, u + 1, v)] !== a) || (v + 1 < H && ARA[vox(plane, k, u, v + 1)] !== a))) { r = 110; g = 170; b = 230; }
    px[o] = r; px[o + 1] = g; px[o + 2] = b; px[o + 3] = A;
  }
  const cv = document.createElement('canvas'); cv.width = W; cv.height = H; cv.getContext('2d').putImageData(img, 0, 0);
  return cv;
}
function fit() {
  const cv = $('slice'), P = PLANES[S.view], W = P.W(), H = P.H();
  const s = Math.min(cv.width / W, cv.height / H) * 0.94 * S.zoom;
  return { s, ox: (cv.width - W * s) / 2 + S.tx, oy: (cv.height - H * s) / 2 + S.ty };
}
function redraw() {
  if (!META) return;
  where();
  if (S.view === 'cut3d') $('posl').textContent = `coronal ${S.pos.coronal} · ${(S.pos.coronal * 0.025).toFixed(2)} mm`;
  if (is3D()) { if (S.view === 'cut3d') { updateCut(); update3DNissl(); } else updateSvg3D(); need3D = true; return; }
  const cv = $('slice'), ctx = cv.getContext('2d'), P = PLANES[S.view];
  ctx.clearRect(0, 0, cv.width, cv.height);
  const { s, ox, oy } = fit();
  ctx.imageSmoothingEnabled = false;
  ctx.drawImage(sliceImage(), ox, oy, P.W() * s, P.H() * s);
  if (S.svg && SVG) drawSvg2D(ctx, s, ox, oy);
  drawCells2D(ctx, s, ox, oy);
  $('posl').textContent = `${S.view} ${S.pos[S.view]} · ${(S.pos[S.view] * 0.025).toFixed(2)} mm`;
}
function where() {
  const el = $('where');
  if (pl() !== 'coronal') { el.innerHTML = `${S.layer === 'full' ? 'filled 3D' : 'plate slices only'} · ${S.view === '3d' ? '3D view' : `${S.view} cut`}`; return; }
  const z = S.pos.coronal, p = META.plates.find((q) => q.z === z), tag = S.view === 'cut3d' ? ' · 3D cut' : '';
  if (p) {
    const fit = EDITION === 'ccfv3bbp' ? `median AP of the plate in CCFv3BBP (its sheet spans ${p.ap_spread_vox} slices)` : `ARA line fit ${p.line_median_um} µm (hippocampus ${p.hpf_median_um} µm)`;
    el.innerHTML = `<b>HGEA level ${p.level}</b> · Bregma ${p.bregma_mm.toFixed(2)} mm<br>${fit}${tag}`; return;
  }
  const a = [...META.plates].reverse().find((q) => q.z < z), b = META.plates.find((q) => q.z > z);
  el.innerHTML = (a && b ? `slice ${z}: <b>filled</b> between levels ${a.level} and ${b.level}` : `slice ${z}: outside the plate series`) + tag;
}

// pan, zoom, hover
(() => {
  const cv = $('slice'); let drag = null;
  cv.addEventListener('wheel', (e) => { e.preventDefault(); const f = e.deltaY < 0 ? 1.15 : 1 / 1.15; S.zoom = Math.max(0.5, Math.min(12, S.zoom * f)); redraw(); }, { passive: false });
  const plateSvg = () => (S.view === 'coronal' && S.svg && SVG ? SVG.find((x) => x.plate.z === S.pos.coronal) : null);
  const at = (e) => { const r = cv.getBoundingClientRect(), { s, ox, oy } = fit(); return [((e.clientX - r.left) * devicePixelRatio - ox) / s - 0.5, ((e.clientY - r.top) * devicePixelRatio - oy) / s - 0.5]; };
  cv.addEventListener('mousedown', (e) => { drag = { x: e.clientX, y: e.clientY, tx: S.tx, ty: S.ty, moved: false }; });
  addEventListener('mouseup', (e) => {
    const d = drag; drag = null;
    if (!d || d.moved || e.target !== cv) return;
    if (SEL) { closeRoi(); return; }   // a click with the popup open: close it, back to the default view
    const sv = plateSvg(); if (!sv) return;
    const roi = roiAt(sv, ...at(e)); if (roi) openRoi(sv, roi);
  });
  cv.addEventListener('mousemove', (e) => {
    if (drag && (drag.moved || Math.hypot(e.clientX - drag.x, e.clientY - drag.y) > 4)) { drag.moved = true; S.tx = drag.tx + (e.clientX - drag.x) * devicePixelRatio; S.ty = drag.ty + (e.clientY - drag.y) * devicePixelRatio; redraw(); return; }
    const r = cv.getBoundingClientRect(), { s, ox, oy } = fit(), P = PLANES[S.view];
    const u = Math.floor(((e.clientX - r.left) * devicePixelRatio - ox) / s), v = Math.floor(((e.clientY - r.top) * devicePixelRatio - oy) / s);
    const tip = $('tip'), sv = SEL ? null : plateSvg(), roi = sv ? roiAt(sv, ...at(e)) : null;
    if ((roi && roi.path) !== (HOV && HOV.roi.path)) { HOV = roi ? { e: sv, roi } : null; cv.style.cursor = roi ? 'pointer' : 'crosshair'; redraw(); }
    if (u < 0 || v < 0 || u >= P.W() || v >= P.H()) { tip.style.display = 'none'; return; }
    const i = vox(S.view, S.pos[S.view], u, v), o = ONTO[ARA[i]];
    const hf = LAB.full[i], hp = LAB.plates[i];
    tip.innerHTML = (roi ? `Bienkowski drawing: <b>${DOMNAME[roi.label] || '—'}</b> · click for its genes and cells<br>` : '') + `<b>${o.acronym}</b> ${o.name}` + (hf ? `<br>HGEA (filled): <b>${DOMNAME[hf]}</b>` : '') + (hp ? `<br>HGEA (plate): <b>${DOMNAME[hp]}</b>` : '')
      + (G.ishAra && G.ishAra[i] ? `<br>ISH ${G.gene}: ${G.ishAra[i]} (${(G.ishAra[i] / G.p995 * 100).toFixed(0)}% of p99.5)` : '')
      + cellTip((e.clientX - r.left) * devicePixelRatio, (e.clientY - r.top) * devicePixelRatio);
    tip.style.display = 'block'; tip.style.left = `${e.clientX - r.left + 14}px`; tip.style.top = `${e.clientY - r.top + 14}px`;
  });
  cv.addEventListener('mouseleave', () => { $('tip').style.display = 'none'; if (HOV) { HOV = null; redraw(); } });
})();

// ───────── 3D ─────────
let cutFace = null, clip = new THREE.Plane(new THREE.Vector3(0, 0, -1), 0);
let renderer, scene, camera, controls, hgPts = null, shell = null, dirty3D = true, need3D = true, vol = null, planes = null, planesDirty = true;
const toScene = (a, d, l) => [(l - NL / 2) * 0.025, -(d - ND / 2) * 0.025, -(a - NA / 2) * 0.025];
function init3D() {
  renderer = new THREE.WebGLRenderer({ canvas: $('gl'), antialias: true, alpha: true });
  renderer.setPixelRatio(Math.min(devicePixelRatio, 2)); renderer.localClippingEnabled = true;
  scene = new THREE.Scene();
  camera = new THREE.PerspectiveCamera(35, 1, 0.01, 200); camera.position.set(9, 6, 8);
  controls = new OrbitControls(camera, $('gl')); controls.enableDamping = true;
  controls.target.set(...toScene(META.plates[16].z, 140, NL / 2)); controls.addEventListener('change', () => { need3D = true; });
  // brain shell: boundary voxels of the annotation, every 3rd voxel
  const pos = [];
  for (let a = 0; a < NA; a += 3) for (let d = 0; d < ND; d += 3) for (let l = 0; l < NL; l += 3) {
    const i = (a * ND + d) * NL + l; if (!ARA[i]) continue;
    if ((l + 3 < NL && !ARA[i + 3]) || (l >= 3 && !ARA[i - 3]) || (d + 3 < ND && !ARA[i + 3 * NL]) || (d >= 3 && !ARA[i - 3 * NL])) pos.push(...toScene(a, d, l));
  }
  const g = new THREE.BufferGeometry(); g.setAttribute('position', new THREE.Float32BufferAttribute(pos, 3));
  shell = new THREE.Points(g, new THREE.PointsMaterial({ size: 0.03, color: 0x8fa3b0, transparent: true, opacity: 0.12, depthWrite: false }));
  shell.visible = S.shell; scene.add(shell);
  scene.add(new THREE.AmbientLight(0xffffff, 0.8));
  const key = new THREE.DirectionalLight(0xffffff, 0.6); key.position.set(4, 8, 6); scene.add(key);
  const fill = new THREE.DirectionalLight(0xffffff, 0.3); fill.position.set(-5, -3, -4); scene.add(fill);
  new ResizeObserver(resize).observe($('stage'));
  (function loop() { requestAnimationFrame(loop); controls.update(); if (need3D && is3D()) { need3D = false; if (vol && vol.visible) vol.material.uniforms.uCam.value.copy(vol.worldToLocal(camera.position.clone())); renderer.render(scene, camera); } })();
}
function build3D() {   // voxel-exact surfaces: every face between a shown domain voxel and a different label (as napari draws labels)
  if (hgPts) { scene.remove(hgPts); hgPts.geometry.dispose(); }
  const lab = LAB[S.layer], show = (h) => h && S.on.has(h % 100);
  const step = [ND * NL, NL, 1], dims = [NA, ND, NL];
  let n = 0;
  const count = (i, c, j, inside) => { if (!inside || lab[j] !== c) n++; };
  for (let a = 0; a < NA; a++) for (let d = 0; d < ND; d++) for (let l = 0; l < NL; l++) {
    const i = (a * ND + d) * NL + l, c = lab[i]; if (!show(c)) continue;
    const q = [a, d, l];
    for (let ax = 0; ax < 3; ax++) { count(i, c, i + step[ax], q[ax] + 1 < dims[ax]); count(i, c, i - step[ax], q[ax] > 0); }
  }
  const pos = new Float32Array(n * 12), nor = new Int8Array(n * 12), col = new Float32Array(n * 12), idx = new Uint32Array(n * 6);
  const lin = (x) => { const c = x / 255; return c <= 0.04045 ? c / 12.92 : ((c + 0.055) / 1.055) ** 2.4; };   // sRGB legend colours -> linear, as three.js expects
  let f = 0;
  const V = 0.025, cx = NL / 2, cy = ND / 2, cz = NA / 2;
  // a face on axis ax at voxel boundary, side s (+1 / -1); corners in grid units (a, d, l)
  const face = (a, d, l, ax, sgn, c) => {
    const o = [a, d, l]; o[ax] += sgn > 0 ? 1 : 0;
    const u = (ax + 1) % 3, w = (ax + 2) % 3, cs = [[0, 0], [1, 0], [1, 1], [0, 1]];
    for (let k = 0; k < 4; k++) {
      const p = [o[0], o[1], o[2]]; p[u] += cs[k][0]; p[w] += cs[k][1];
      const vi = (f * 4 + k) * 3;
      pos[vi] = (p[2] - cx) * V; pos[vi + 1] = -(p[1] - cy) * V; pos[vi + 2] = -(p[0] - cz) * V;
      const nv = [0, 0, 0]; nv[ax] = sgn;   // grid normal -> scene normal (x = LR, y = -DV, z = -AP)
      nor[vi] = nv[2] * 127; nor[vi + 1] = -nv[1] * 127; nor[vi + 2] = -nv[0] * 127;
      col[vi] = lin(DOMRGB[c * 3]); col[vi + 1] = lin(DOMRGB[c * 3 + 1]); col[vi + 2] = lin(DOMRGB[c * 3 + 2]);
    }
    const b0 = f * 4, ii = f * 6;
    // winding so the face is front-facing along its normal (parity of the axis permutation decides)
    if ((sgn > 0) === (ax !== 1)) { idx.set([b0, b0 + 1, b0 + 2, b0, b0 + 2, b0 + 3], ii); } else { idx.set([b0, b0 + 2, b0 + 1, b0, b0 + 3, b0 + 2], ii); }
    f++;
  };
  for (let a = 0; a < NA; a++) for (let d = 0; d < ND; d++) for (let l = 0; l < NL; l++) {
    const i = (a * ND + d) * NL + l, c = lab[i]; if (!show(c)) continue;
    const q = [a, d, l];
    for (let ax = 0; ax < 3; ax++) {
      if (!(q[ax] + 1 < dims[ax]) || lab[i + step[ax]] !== c) face(a, d, l, ax, 1, c);
      if (!(q[ax] > 0) || lab[i - step[ax]] !== c) face(a, d, l, ax, -1, c);
    }
  }
  const g = new THREE.BufferGeometry();
  g.setAttribute('position', new THREE.BufferAttribute(pos, 3));
  g.setAttribute('normal', new THREE.BufferAttribute(nor, 3, true));
  g.setAttribute('color', new THREE.BufferAttribute(col, 3));
  g.setIndex(new THREE.BufferAttribute(idx, 1));
  hgPts = new THREE.Mesh(g, new THREE.MeshLambertMaterial({ vertexColors: true, side: THREE.DoubleSide, transparent: S.opacity < 1, opacity: S.opacity,
    clippingPlanes: S.view === 'cut3d' ? [clip] : [] }));
  scene.add(hgPts); dirty3D = false; need3D = true;
}
function resize() {
  const st = $('stage'), w = st.clientWidth, h = st.clientHeight;
  if (renderer) { renderer.setSize(w, h, false); camera.aspect = w / h; camera.updateProjectionMatrix(); need3D = true; }
  const sc = $('slice'); sc.width = w * devicePixelRatio; sc.height = h * devicePixelRatio; redraw();
}

// 3D Nissl: volume rendering (front-to-back compositing, 50 µm) or three orthogonal slices (25 µm)
const VS = `out vec3 vPos; void main() { vPos = position; gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0); }`;
const FS = `precision highp float; precision highp sampler3D;
in vec3 vPos; out vec4 outColor;
uniform sampler3D uVol; uniform vec3 uCam; uniform float uDensity; uniform float uClipZ;
vec2 hitBox(vec3 o, vec3 d) { vec3 inv = 1.0 / d, t0 = (vec3(-0.5) - o) * inv, t1 = (vec3(0.5) - o) * inv;
  vec3 tmin = min(t0, t1), tmax = max(t0, t1); return vec2(max(max(tmin.x, tmin.y), tmin.z), min(min(tmax.x, tmax.y), tmax.z)); }
void main() {
  vec3 dir = normalize(vPos - uCam); vec2 t = hitBox(uCam, dir); t.x = max(t.x, 0.0);
  if (t.x >= t.y) discard;
  float dt = 1.0 / 320.0; vec4 acc = vec4(0.0);
  for (float s = t.x; s < t.y; s += dt) {
    vec3 p = uCam + dir * s;
    if (p.z > uClipZ) continue;
    float v = texture(uVol, vec3(p.x + 0.5, 0.5 - p.y, 0.5 - p.z)).r;
    float a = smoothstep(0.08, 0.6, v) * uDensity * 0.08;
    vec3 c = vec3(0.35 + 0.9 * v);
    acc.rgb += (1.0 - acc.a) * a * c; acc.a += (1.0 - acc.a) * a;
    if (acc.a > 0.97) break;
  }
  outColor = acc;
}`;
function buildVolume() {
  const W = NL >> 1, H = ND >> 1, D = NA >> 1, v = new Uint8Array(W * H * D);
  for (let z = 0; z < D; z++) for (let y = 0; y < H; y++) for (let x = 0; x < W; x++) {
    let s = 0;
    for (let dz = 0; dz < 2; dz++) for (let dy = 0; dy < 2; dy++) { const b = ((2 * z + dz) * ND + 2 * y + dy) * NL + 2 * x; s += NISSL[b] + NISSL[b + 1]; }
    v[(z * H + y) * W + x] = s >> 3;
  }
  const tex = new THREE.Data3DTexture(v, W, H, D);
  tex.format = THREE.RedFormat; tex.minFilter = tex.magFilter = THREE.LinearFilter; tex.unpackAlignment = 1; tex.needsUpdate = true;
  const mat = new THREE.ShaderMaterial({ glslVersion: THREE.GLSL3, vertexShader: VS, fragmentShader: FS, side: THREE.BackSide,
    transparent: true, depthWrite: false, uniforms: { uVol: { value: tex }, uCam: { value: new THREE.Vector3() }, uDensity: { value: S.density }, uClipZ: { value: 1.0 } } });
  vol = new THREE.Mesh(new THREE.BoxGeometry(1, 1, 1), mat);
  vol.scale.set(NL * 0.025, ND * 0.025, NA * 0.025); vol.renderOrder = -1;
  scene.add(vol);
}
function planeMesh(plane) {
  const P = PLANES[plane], W = P.W(), H = P.H(), k = S.pos[plane], data = new Uint8Array(W * H * 4);
  for (let v = 0; v < H; v++) for (let u = 0; u < W; u++) {
    const i = vox(plane, k, u, v), o = (v * W + u) * 4, g = Math.min(255, 30 + NISSL[i] * 1.1);
    data[o] = data[o + 1] = data[o + 2] = g; data[o + 3] = ARA[i] ? Math.round(255 * Math.min(1, 0.25 + S.density)) : 0;
  }
  const tex = new THREE.DataTexture(data, W, H); tex.flipY = false; tex.magFilter = THREE.LinearFilter; tex.needsUpdate = true;
  const m = new THREE.Mesh(new THREE.PlaneGeometry(W * 0.025, H * 0.025),
    new THREE.MeshBasicMaterial({ map: tex, transparent: true, side: THREE.DoubleSide, depthWrite: false }));
  // plane texture: u right, v down; place each plane at its cut
  if (plane === 'coronal') { m.position.set(0, 0, -(k - NA / 2) * 0.025); m.scale.y = -1; }
  else if (plane === 'sagittal') { m.rotation.y = Math.PI / 2; m.position.set((k - NL / 2) * 0.025, 0, 0); m.scale.y = -1; }
  else { m.rotation.x = -Math.PI / 2; m.position.set(0, -(k - ND / 2) * 0.025, 0); m.scale.y = -1; }
  return m;
}
// coronal 3D: the atlas cut at the coronal position; the cut face shows that slice exactly as in the 2D view
function updateCut() {
  if (!scene) return;
  const k = S.pos.coronal, zc = -(k - NA / 2) * 0.025;
  clip.constant = zc - 0.0125;   // mesh keeps everything posterior to the slice; the cut face shows the slice itself
  if (hgPts) hgPts.material.clippingPlanes = [clip];
  const tex = new THREE.CanvasTexture(sliceImage('coronal', k, 150));   // tissue see-through, domains solid
  tex.magFilter = THREE.NearestFilter; tex.colorSpace = THREE.SRGBColorSpace;
  if (!cutFace) {
    cutFace = new THREE.Mesh(new THREE.PlaneGeometry(NL * 0.025, ND * 0.025), new THREE.MeshBasicMaterial({ map: tex, transparent: true, side: THREE.DoubleSide }));
    scene.add(cutFace);
  } else { cutFace.material.map.dispose(); cutFace.material.map = tex; }
  cutFace.visible = true; cutFace.position.set(0, 0, zc);
  cutCells(k, zc);
  need3D = true;
}
function aimCut() {
  updateCut();
  const zc = -(S.pos.coronal - NA / 2) * 0.025;
  controls.target.set(0, 0.3, zc - 0.8); camera.position.set(8.5, 6, zc + 13); controls.update();
}

function update3DNissl() {
  if (!scene) return;
  updateSvg3D();
  const cut = S.view === 'cut3d', clips = cut ? [clip] : [];
  if (!cut) { if (cutFace) cutFace.visible = false; if (hgPts) hgPts.material.clippingPlanes = []; if (cutPts) { scene.remove(cutPts); cutPts.geometry.dispose(); cutPts = null; } }
  if (S.nissl3d === 'volume' && !vol) buildVolume();
  if (vol) {
    vol.visible = S.nissl3d === 'volume'; vol.material.uniforms.uDensity.value = S.density;
    vol.material.uniforms.uClipZ.value = cut ? (clip.constant / vol.scale.z) : 1.0;   // world z of the cut -> box-local z
  }
  if (planes && (planesDirty || S.nissl3d !== 'planes')) { planes.forEach((m) => { scene.remove(m); m.geometry.dispose(); m.material.map.dispose(); }); planes = null; }
  if (S.nissl3d === 'planes' && !planes) { planes = ['coronal', 'sagittal', 'horizontal'].map(planeMesh); planes.forEach((m) => scene.add(m)); planesDirty = false; }
  if (planes) planes.forEach((m, k) => { m.material.clippingPlanes = clips; m.visible = !(cut && k === 0); });   // the cut face replaces the coronal slice
  if (shell) shell.material.clippingPlanes = clips;
  if (cellPts) cellPts.material.clippingPlanes = clips;
  need3D = true;
}

// ───────── Bienkowski SVG plates (original drawings, aligned by the same transform the build used) ─────────
let SVG = null, svgPlanes = null, svgLoading = null;
const svgClip = new THREE.Plane(new THREE.Vector3(0, 0, -1), 0);
const R = 2;   // raster resolution for 3D and for sagittal / horizontal crossings: R pixels per 25 µm voxel
function loadSvgs() {
  if (svgLoading) return svgLoading;
  $('where').innerHTML += '<br>loading the 33 HGEA SVG plates…';
  svgLoading = Promise.all(META.plates.map((p) => new Promise((ok, bad) => {
    const im = new Image(); im.onload = () => ok(im); im.onerror = () => bad(new Error(p.svg)); im.src = DATA + p.svg;
  }))).then(async (imgs) => {
    const texts = await Promise.all(META.plates.map((p) => fetch(DATA + p.svg).then((r) => r.text()).catch(() => '')));
    SVG = META.plates.map((p, k) => {
      const [A, B, C, D, E, F] = p.svg_affine, cv = document.createElement('canvas');
      cv.width = NL * R; cv.height = ND * R;
      const ctx = cv.getContext('2d'), o = p.brain_outline;
      ctx.beginPath();   // keep the drawing on the brain: drops the legend, inset, axes and labels around the plate
      for (let j = 0; j < o.length; j += 2) ctx[j ? 'lineTo' : 'moveTo'](R * (o[j] + 0.5), R * (o[j + 1] + 0.5));
      ctx.closePath(); ctx.clip();
      ctx.setTransform(R * A, R * D, R * B, R * E, R * (C + 0.5), R * (F + 0.5));   // SVG px -> raster px (voxel centres at +0.5)
      ctx.drawImage(imgs[k], 0, 0);
      return { plate: p, img: imgs[k], raster: cv, data: ctx.getImageData(0, 0, cv.width, cv.height).data, shapes: svgShapes(texts[k]), rois: null };
    });
    where();
  });
  return svgLoading;
}
function drawSvg2D(ctx, s, ox, oy) {
  ctx.save(); ctx.globalAlpha = S.focus && G.gene ? Math.min(0.2, S.svgAlpha) : S.svgAlpha;
  if (S.view === 'coronal') {
    const e = SVG.find((x) => x.plate.z === S.pos.coronal), hv = SEL || HOV;
    if (e) {
      const [A, B, C, D, E, F] = e.plate.svg_affine, on = hv && hv.e === e;
      ctx.imageSmoothingEnabled = true;
      ctx.setTransform(s * A, s * D, s * B, s * E, ox + s * (C + 0.5), oy + s * (F + 0.5));
      ctx.drawImage(e.img, 0, 0);
      if (on && !(S.focus && G.gene)) {   // the hovered region: opaque, lightened and outlined
        ctx.save(); ctx.clip(hv.roi.path, hv.roi.rule); ctx.globalAlpha = 1; ctx.drawImage(e.img, 0, 0);
        ctx.fillStyle = 'rgba(255,255,255,.25)'; ctx.fill(hv.roi.path, hv.roi.rule); ctx.restore();
      }
      if (on) {
        ctx.globalAlpha = 1; ctx.strokeStyle = '#fff'; ctx.lineWidth = 2.2 * devicePixelRatio / (s * Math.hypot(A, D)); ctx.stroke(hv.roi.path);
      }
    }
  } else {   // a 1-voxel strip of each plate's drawing at its slice
    const P = PLANES[S.view], W = P.W(), H = P.H(), k = S.pos[S.view], img = new ImageData(W, H), px = img.data, RW = NL * R;
    for (const e of SVG) {
      const z = e.plate.z;
      for (let t = 0; t < (S.view === 'sagittal' ? ND : NL); t++) {
        const x = S.view === 'sagittal' ? k : t, y = S.view === 'sagittal' ? t : k;
        const si = ((y * R + (R >> 1)) * RW + x * R + (R >> 1)) * 4;
        if (e.data[si + 3] < 20) continue;
        const u = S.view === 'sagittal' ? z : t, v = S.view === 'sagittal' ? t : z, o = (v * W + u) * 4;
        px[o] = e.data[si]; px[o + 1] = e.data[si + 1]; px[o + 2] = e.data[si + 2]; px[o + 3] = 255;
      }
    }
    const cv = document.createElement('canvas'); cv.width = W; cv.height = H; cv.getContext('2d').putImageData(img, 0, 0);
    ctx.imageSmoothingEnabled = false; ctx.drawImage(cv, ox, oy, W * s, H * s);
  }
  ctx.restore();
}
// the drawing's coloured regions (the HGEA domains as Bienkowski drew them), for hover and click in coronal 2D
function svgShapes(text) {
  if (!text) return [];
  const doc = new DOMParser().parseFromString(text, 'image/svg+xml'), out = [];
  for (const el of doc.querySelectorAll('path, polygon')) {
    const fill = (el.getAttribute('fill') || (el.getAttribute('style') || '').match(/fill:\s*(#[0-9a-fA-F]{6})/)?.[1] || '').toLowerCase();
    if (!/^#[0-9a-f]{6}$/.test(fill) || (fill.slice(1, 3) === fill.slice(3, 5) && fill.slice(3, 5) === fill.slice(5, 7))) continue;   // greys: not a domain
    const d = el.tagName === 'polygon' ? `M${el.getAttribute('points').trim().replace(/\s+/g, ' L')}Z` : el.getAttribute('d');
    if (d) out.push({ path: new Path2D(d), rule: el.getAttribute('fill-rule') === 'evenodd' ? 'evenodd' : 'nonzero', fill });
  }
  return out;
}
const svgCtx = document.createElement('canvas').getContext('2d');
function toSvg(e, u, v) {   // voxel-centre coordinates (u = LR, v = DV) on the plate slice -> SVG px
  const [A, B, C, D, E, F] = e.plate.svg_affine, det = A * E - B * D;
  return [(E * (u - C) - B * (v - F)) / det, (-D * (u - C) + A * (v - F)) / det];
}
function roisOf(e) {   // each coloured region -> the plate voxels inside it and its HGEA domain (majority label)
  if (e.rois) return e.rois;
  const z = e.plate.z, lab = LAB.plates, vox = [];
  for (let v = 0; v < ND; v++) for (let u = 0; u < NL; u++) { const i = (z * ND + v) * NL + u; if (lab[i]) vox.push(u, v, i); }
  e.rois = [];
  for (const sh of e.shapes) {
    const inside = [], n = new Map();
    for (let j = 0; j < vox.length; j += 3) {
      const [x, y] = toSvg(e, vox[j], vox[j + 1]);
      if (!svgCtx.isPointInPath(sh.path, x, y, sh.rule)) continue;
      inside.push(vox[j + 2]); const h = lab[vox[j + 2]]; n.set(h, (n.get(h) || 0) + 1);
    }
    if (inside.length < 3) continue;   // legend swatches and slivers
    const [h, c] = [...n].sort((a, b) => b[1] - a[1])[0];
    e.rois.push({ ...sh, vox: inside, label: h, purity: c / inside.length });
  }
  return e.rois;
}
function roiAt(e, u, v) {   // smallest region under a point
  const [x, y] = toSvg(e, u, v);
  return roisOf(e).filter((r) => svgCtx.isPointInPath(r.path, x, y, r.rule)).sort((a, b) => a.vox.length - b.vox.length)[0] || null;
}

// region popup: Table 2 marker genes, the gene's ISH and MERFISH inside the region, the cell types there
let HOV = null, SEL = null, T2 = null;
const t2cols = (name) => name.replace(' (left)', '').split('/');
async function openRoi(e, roi) {
  SEL = { e, roi }; HOV = null; $('tip').style.display = 'none';
  if (!T2) T2 = await fetch(SITE_DATA + 'hgea/hgea.json').then((r) => r.json()).catch(() => ({ columns: [], genes: [] }));
  let u0 = 1e9, u1 = -1, v0 = 1e9, v1 = -1;   // zoom to the region, left of the popup
  for (const i of roi.vox) { const r = i % (ND * NL), u = r % NL, v = Math.floor(r / NL); u0 = Math.min(u0, u); u1 = Math.max(u1, u); v0 = Math.min(v0, v); v1 = Math.max(v1, v); }
  const cv = $('slice'), s0 = Math.min(cv.width / NL, cv.height / ND) * 0.94;
  const s = Math.min(cv.width * 0.45 / (u1 - u0 + 1), cv.height * 0.6 / (v1 - v0 + 1), s0 * 8);
  S.zoom = Math.max(1, s / s0);
  const sz = s0 * S.zoom;
  S.tx = cv.width * 0.36 - ((u0 + u1) / 2 + 0.5) * sz - (cv.width - NL * sz) / 2;
  S.ty = cv.height * 0.5 - ((v0 + v1) / 2 + 0.5) * sz - (cv.height - ND * sz) / 2;
  renderRoi(); redraw();
}
function closeRoi(reset = true) {   // back to the default view: no zoom, no gene
  SEL = null; S.focus = false; $('roipop').style.display = 'none';
  if (!reset) return;
  S.zoom = 1; S.tx = S.ty = 0;
  if (G.gene) { $('genein').value = ''; clearGene(); } else redraw();
}
function renderRoi() {
  if (!SEL) return;
  const { e, roi } = SEL, p = e.plate, name = DOMNAME[roi.label] || '—';
  const cols = t2cols(name), ks = cols.map((c) => T2.columns.indexOf(c)).filter((k) => k >= 0);
  const markers = T2.genes.filter((g) => g.domains.some((k) => ks.includes(k)));
  // cells in the slab whose position falls inside the drawn region
  const [a, b] = apRange(p.z * 25 - S.slab, p.z * 25 + S.slab), here = [], per = new Map();
  let hp = 0, hpDet = 0;
  for (let j = a; j < b; j++) {
    const i = C.order[j], u = C.pos[i * 3 + 2] / 25, v = C.pos[i * 3 + 1] / 25;
    if (C.hgea[i] && G.expr) { hp++; if (G.expr[i]) hpDet++; }
    const [x, y] = toSvg(e, u, v);
    if (!svgCtx.isPointInPath(roi.path, x, y, roi.rule)) continue;
    here.push(i); per.set(C.sub[i], (per.get(C.sub[i]) || 0) + 1);
  }
  const top = [...per].sort((x, y) => y[1] - x[1]).slice(0, 10);
  const hex = (k) => C.tax.subclass[k].color;
  let gene = '';
  if (G.gene) {
    const rows = [];
    if (G.ishAra && $('ishon')) {
      let sr = 0, sa = 0, na = 0;
      for (const i of roi.vox) sr += G.ishAra[i];
      for (let i = p.z * ND * NL, n = (p.z + 1) * ND * NL; i < n; i++) if (LAB.plates[i]) { sa += G.ishAra[i]; na++; }
      const mr = sr / roi.vox.length, ma = sa / Math.max(1, na);
      rows.push(['ISH mean in the region', `${mr.toFixed(1)} (${(mr / G.p995 * 100).toFixed(0)}% of p99.5)`], ['vs the whole HGEA on this plate', `${(mr / Math.max(0.01, ma)).toFixed(2)}×`]);
    } else rows.push(['ISH', 'no Blue Brain ISH volume']);
    if (G.expr) {
      const det = here.filter((i) => G.expr[i]), m = det.reduce((t, i) => t + G.expr[i], 0) / Math.max(1, det.length) / 255 * G.vmax;
      rows.push(['MERFISH cells detecting it', here.length ? `${(det.length / here.length * 100).toFixed(0)}% of ${here.length.toLocaleString()}` : '—'],
        ['  in the whole HGEA (same slab)', hp ? `${(hpDet / hp * 100).toFixed(0)}% of ${hp.toLocaleString()}` : '—'],
        ['mean log2 where detected', det.length ? m.toFixed(2) : '—']);
    } else rows.push(['MERFISH', `not on the ${C.list[S.brain].short} panel`]);
    gene = `<div class="h">${G.gene} in this region</div><table>${rows.map(([k, v]) => `<tr><td>${k}</td><td class="n">${v}</td></tr>`).join('')}</table>`;
  }
  $('roipop').innerHTML = `<button class="x" id="roix" title="close: back to the default view">×</button>
    <h3>${name}</h3>
    <div class="m">HGEA level ${p.level} · Bregma ${p.bregma_mm.toFixed(2)} mm · region in Bienkowski's drawing${roi.purity < 0.8 ? ` · ${(roi.purity * 100).toFixed(0)}% of its voxels carry this label` : ''}</div>
    ${gene}
    <div class="h">MERFISH ${C.list[S.brain].short} · cells within ±${S.slab} µm, inside the region</div>
    ${here.length ? `<table>${top.map(([k, n]) => `<tr><td><span class="sw" style="background:${hex(k)}"></span>${C.tax.subclass[k].name}</td><td class="n">${n.toLocaleString()}</td><td class="n">${(n / here.length * 100).toFixed(0)}%</td></tr>`).join('')}
      <tr><td>all cells</td><td class="n">${here.length.toLocaleString()}</td><td></td></tr></table>` : '<div class="m">no cells of this brain here</div>'}
    <div class="h">Marker genes · Table 2 (${cols.join(' + ')}) · ${markers.length}</div>
    <div class="gs">${markers.map((g) => { const sym = G.symbols.includes(g.symbol) ? g.symbol : G.symbols.includes(g.ish) ? g.ish : null;
      return `<span class="g${sym ? '' : ' dis'}${sym && sym === G.gene ? ' on' : ''}" ${sym ? `data-g="${sym}"` : 'title="not in the ISH atlas or the MERFISH panels"'}>${g.name}</span>`; }).join('') || '<span class="m">none listed</span>'}</div>
    <div class="m" style="margin-top:8px">click a gene to highlight where it is (others fade) · again to turn it off · click the slice or × to close</div>`;
  $('roipop').style.display = 'block';
  $('roix').onclick = () => closeRoi();
  $('roipop').querySelectorAll('.g[data-g]').forEach((el) => { el.onclick = () => {
    if (el.dataset.g === G.gene) { $('genein').value = ''; clearGene(); return; }   // click the shown gene again: off
    S.focus = true; $('genein').value = el.dataset.g; selectGene(el.dataset.g);
  }; });
}

function updateSvg3D() {
  if (!scene) return;
  if (S.svg && SVG && !svgPlanes) {
    svgPlanes = SVG.map((e) => {
      const tex = new THREE.CanvasTexture(e.raster); tex.colorSpace = THREE.SRGBColorSpace;
      const m = new THREE.Mesh(new THREE.PlaneGeometry(NL * 0.025, ND * 0.025),
        new THREE.MeshBasicMaterial({ map: tex, transparent: true, opacity: S.svgAlpha, side: THREE.DoubleSide, depthWrite: false }));
      m.position.set(0, 0, -(e.plate.z - NA / 2) * 0.025 + 0.003);   // just anterior of the slice, so it draws over the cut face
      scene.add(m); return m;
    });
  }
  if (!svgPlanes) return;
  const cut = S.view === 'cut3d';
  svgClip.constant = clip.constant + 0.025 + 0.004;   // in coronal 3D keep the plate at the cut and everything posterior
  svgPlanes.forEach((m) => { m.visible = S.svg; m.material.clippingPlanes = cut ? [svgClip] : []; });
  need3D = true;
}

// ───────── MERFISH cells mapped into the 3D ARA ─────────
// pos: 3D ARA µm (AP, DV, LR) from this site's CCFv3BBP positions, by the same transform as the atlas.
const C = { list: null, tax: null, sub: null, pos: null, sub_of: null, hgea: null, n: 0, order: null, hits: [] };
const HIPP_TYPES = { '037 DG Glut': '#e6194b', '016 CA1-ProS Glut': '#3cb44b', '017 CA3 Glut': '#4363d8', '025 CA2-FC-IG Glut': '#f58231',
  '023 SUB-ProS Glut': '#911eb4', '031 CT SUB Glut': '#42d4f4', '033 NP SUB Glut': '#f032e6' };
S.cells = true; S.slab = 200; S.cellColour = 'all'; S.cellAlpha = 1; S.cells3d = false; S.brain = null; S.hgeaLines = true;
function initModes() {
  document.querySelectorAll('#mode button').forEach((b) => b.addEventListener('click', () => {
    if (b.dataset.v === 'hgeaara') return;
    if (window.parent !== window) window.parent.postMessage({ hgeaMode: b.dataset.v }, '*');
    else location.href = `../atlas/index.html?mode=${b.dataset.v}`;
  }));
}

// ───────── gene lookup: MERFISH expression per cell, Blue Brain ISH carried into the ARA around the hippocampus ─────────
const G = { ish: null, lookup: null, symbols: [], gene: null, expr: null, ishAra: null, p995: 1, panel: {} };
async function initGenes() {
  G.ish = await fetch(SITE_DATA + 'ish/genes.json').then((r) => r.json());
  const ishSet = new Set(G.ish.genes.map((g) => g.gene));
  for (const k of Object.keys(C.list)) {
    const m = await fetch(SITE_DATA + `merfish/${k}/manifest.json`).then((r) => r.json());
    G.panel[k] = new Map(m.genes.filter((g) => !/^Blank-/.test(g.symbol)).map((g) => [g.symbol, g]));
  }
  const all = new Set([...ishSet, ...Object.values(G.panel).flatMap((p) => [...p.keys()])]);
  G.symbols = [...all].sort((a, b) => a.localeCompare(b));
  $('genelist').innerHTML = G.symbols.map((g) => `<option value="${g}">`).join('');
  $('genein').addEventListener('change', () => { const g = $('genein').value.trim(); S.focus = !!SEL; if (g) selectGene(g); else if (G.gene) clearGene(); });
  $('genein').addEventListener('input', () => { if (!$('genein').value.trim() && G.gene) clearGene(); });
  $('geneclear').onclick = () => { $('genein').value = ''; clearGene(); };
  addEventListener('message', (e) => {   // the atlas page's gene, when it changes there
    const g = e.data && e.data.hgeaGene; if (g === undefined || g === G.gene || (!g && !G.gene)) return;
    $('genein').value = g || ''; if (g && G.symbols.includes(g)) selectGene(g); else clearGene();
  });
  $('ishon').onchange = () => { redraw(); };
  $('exprcells').onchange = () => { colourCells(); redraw(); build3DCells(); if (S.view === 'cut3d') updateCut(); };
  const q = new URLSearchParams(location.search).get('gene');
  if (q && all.has(q.split(',')[0])) { $('genein').value = q.split(',')[0]; selectGene(q.split(',')[0]); }
}
async function selectGene(g) {
  if (!G.symbols.includes(g)) { $('genestatus').textContent = `"${g}" is not in the ISH atlas or any MERFISH panel`; return; }
  G.gene = g; G.expr = null; G.ishAra = null;
  const bits = [];
  const ishEntry = G.ish.genes.find((x) => x.gene === g);
  $('genestatus').textContent = `loading ${g}…`;
  if (ishEntry) {
    if (!G.lookup) {
      const [a, c] = await Promise.all(['ara', 'ccf'].map((f) => gz(`ish_lookup.${f}.u32.gz`)));
      G.lookup = { ara: new Uint32Array(a), ccf: new Uint32Array(c) };
    }
    const vol = await IO.ishVolume(G.ish.url_template, g);
    const out = new Uint8Array(NA * ND * NL);
    for (let j = 0; j < G.lookup.ara.length; j++) out[G.lookup.ara[j]] = vol[G.lookup.ccf[j]];
    G.ishAra = out; G.p995 = Math.max(1, ishEntry.p995 || 1);
    bits.push('ISH ✓ (hippocampal region)');
  } else bits.push('no ISH volume');
  await loadExpr();
  bits.push(G.expr ? `MERFISH ${C.list[S.brain].short} ✓` : `not on the ${C.list[S.brain].short} panel`);
  $('genestatus').textContent = `${g}: ${bits.join(' · ')}`;
  colourCells(); redraw(); build3DCells(); if (S.view === 'cut3d') updateCut(); renderRoi();
}
function clearGene() {   // no gene: back to the cell types of the chosen brain, no ISH
  S.focus = false; G.gene = null; G.expr = null; G.ishAra = null;
  $('genestatus').textContent = 'type a gene symbol, e.g. Prox1, Wfs1, Dsp';
  colourCells(); redraw(); build3DCells(); if (S.view === 'cut3d') updateCut(); renderRoi();
}
async function loadExpr() {
  G.expr = null;
  if (!G.gene || !C.idx || !G.panel[S.brain] || !G.panel[S.brain].has(G.gene)) return;
  try {
    const raw = new Uint8Array(await IO.fetchGz(`${CONFIG.EXPR}merfish/${S.brain}/expr/${encodeURIComponent(G.gene)}.u8.gz`));
    G.expr = new Uint8Array(C.n);
    for (let i = 0; i < C.n; i++) G.expr[i] = raw[C.idx[i]];
    G.vmax = G.panel[S.brain].get(G.gene).vmax;
  } catch (e) { console.warn(e); }
}

async function initCells() {
  initModes();
  [C.list, C.tax] = await Promise.all([json('cells.json'), fetch(SITE_DATA + 'merfish/taxonomy.json').then((r) => r.json())]);
  C.sub_of = Uint16Array.from(C.tax.cluster.map((c) => C.tax.supertype[c.parent].parent));
  C.subRGB = C.tax.subclass.map((x) => [1, 3, 5].map((k) => parseInt(x.color.slice(k, k + 2), 16)));
  C.hippRGB = C.tax.subclass.map((x) => (HIPP_TYPES[x.name] ? [1, 3, 5].map((k) => parseInt(HIPP_TYPES[x.name].slice(k, k + 2), 16)) : null));
  const sel = $('brain');
  sel.innerHTML = Object.entries(C.list).map(([k, d]) => `<option value="${k}">${d.short} · ${d.cells.toLocaleString()} cells</option>`).join('');
  sel.onchange = () => loadBrain(sel.value);
  const restyle = () => { if (S.view === 'cut3d') updateCut(); };
  $('cells').onchange = () => { S.cells = $('cells').checked; redraw(); build3DCells(); restyle(); };
  $('slab').onchange = () => { S.slab = +$('slab').value; redraw(); restyle(); };
  const legend = () => { $('celllegend').style.display = S.cellColour === 'hipp' ? '' : 'none'; };
  $('cellcolour').onchange = () => { S.cellColour = $('cellcolour').value; legend(); colourCells(); redraw(); build3DCells(); restyle(); };
  $('cells3d').onchange = () => { S.cells3d = $('cells3d').checked; build3DCells(); };
  $('cellop').oninput = () => {
    S.cellAlpha = $('cellop').value / 100;
    for (const m of [cellPts, cutPts]) if (m) { m.material.opacity = S.cellAlpha; m.material.transparent = S.cellAlpha < 1; }
    need3D = true; redraw();
  };
  $('hgealines').onchange = () => { S.hgeaLines = $('hgealines').checked; redraw(); };
  $('celllegend').innerHTML = Object.entries(HIPP_TYPES).map(([n, c]) => `<div class="dom" style="cursor:default"><span class="sw" style="background:${c}"></span>${n}</div>`).join('');
  $('cellcolour').value = S.cellColour; legend();
  const acc = await json('accuracy.json');
  $('cellacc').innerHTML = '<tr><th>brain</th><th class="n">fields agree</th><th class="n">DG Glut in DG</th><th class="n">CA1-ProS in CA1</th></tr>' +
    Object.entries(acc).map(([k, a]) => `<tr><td>${C.list[k].short}</td><td class="n">${(a.region_agreement * 100).toFixed(0)}%</td>` +
      `<td class="n">${(a.cell_types['037 DG Glut'].in_ara_field * 100).toFixed(0)}%</td><td class="n">${(a.cell_types['016 CA1-ProS Glut'].in_ara_field * 100).toFixed(0)}%</td></tr>`).join('');
  const want = new URLSearchParams(location.search).get('brain');
  await loadBrain(C.list[want] ? want : C.list.zhuang1 ? 'zhuang1' : Object.keys(C.list)[0]);   // Zhuang 1 by default
  await initGenes();
}
async function loadBrain(k) {
  S.brain = k; $('brain').value = k; $('cellstatus').textContent = `loading ${C.list[k].short}…`;
  const [p, c, h, ix] = await Promise.all([['pos', 16], ['cluster', 16], ['hgea', 8], ['idx', 32]].map(([f, b]) => gz(`merfish/${k}/${f}.u${b}.gz`)));
  C.pos = new Uint16Array(p); C.sub = new Uint16Array(c).map((x) => C.sub_of[x]); C.hgea = new Uint8Array(h); C.idx = new Uint32Array(ix); C.n = C.sub.length;
  await loadExpr();
  C.order = Uint32Array.from({ length: C.n }, (_, i) => i).sort((a, b) => C.pos[a * 3] - C.pos[b * 3]);   // by AP, for slab lookups
  colourCells();
  $('cellstatus').textContent = `${C.list[k].label}: ${C.n.toLocaleString()} cells in the hippocampal region, mapped into the 3D ARA`;
  redraw(); build3DCells(); if (S.view === 'cut3d') updateCut();
}
function colourCells() {
  C.rgb = new Uint8Array(C.n * 3); C.show = new Uint8Array(C.n);
  if (G.expr && $('exprcells').checked) {   // gene expression (magma), cells without detection dim
    const L = lut('magma'), hi = 0.6 * 255;
    for (let i = 0; i < C.n; i++) {
      const v = G.expr[i];
      if (!v) { C.rgb.set([60, 66, 72], i * 3); C.show[i] = 1; continue; }
      const k = Math.max(64, Math.min(255, Math.round(v / hi * 255))) * 3;
      C.rgb.set([L[k], L[k + 1], L[k + 2]], i * 3); C.show[i] = 1;
    }
    return;
  }
  for (let i = 0; i < C.n; i++) {
    const c = S.cellColour === 'hipp' ? C.hippRGB[C.sub[i]] : C.subRGB[C.sub[i]];
    if (c) { C.rgb.set(c, i * 3); C.show[i] = 1; }
  }
}
function apRange(lo, hi) {   // indices into C.order with AP µm in [lo, hi]
  const lb = (x) => { let a = 0, b = C.n; while (a < b) { const m = (a + b) >> 1; if (C.pos[C.order[m] * 3] < x) a = m + 1; else b = m; } return a; };
  return [lb(lo), lb(hi + 1)];
}
function drawCells2D(ctx, s, ox, oy) {
  C.hits = [];
  if (!S.cells || !C.pos) return;
  const k = S.pos[S.view] * 25, slab = S.slab, ax = S.view === 'coronal' ? 0 : S.view === 'sagittal' ? 2 : 1;
  const fo = S.focus && G.gene, r0 = Math.max(1.2 * devicePixelRatio, s * 0.16);
  let r = r0, pass = 0;   // focus: pass 0 the other cells, faint; pass 1 the cells detecting the gene, larger
  const each = (i) => {
    if (!C.show[i] || Math.abs(C.pos[i * 3 + ax] - k) > slab) return;
    if (fo && (pass === 1) !== !!(G.expr && G.expr[i])) return;
    const ap = C.pos[i * 3] / 25, dv = C.pos[i * 3 + 1] / 25, lr = C.pos[i * 3 + 2] / 25;
    const [u, v] = S.view === 'coronal' ? [lr, dv] : S.view === 'sagittal' ? [ap, dv] : [lr, ap];
    const x = ox + (u + 0.5) * s, y = oy + (v + 0.5) * s;
    ctx.fillStyle = `rgb(${C.rgb[i * 3]},${C.rgb[i * 3 + 1]},${C.rgb[i * 3 + 2]})`;
    ctx.fillRect(x - r, y - r, 2 * r, 2 * r);
    if (C.hits.length < 400000) C.hits.push(i, x, y);
  };
  const run = () => { if (ax === 0) { const [a, b] = apRange(k - slab, k + slab); for (let j = a; j < b; j++) each(C.order[j]); } else for (let i = 0; i < C.n; i++) each(i); };
  ctx.save(); ctx.globalAlpha = S.cellAlpha;
  if (!fo) run();
  else { ctx.globalAlpha = 0.12 * S.cellAlpha; run(); if (G.expr) { ctx.globalAlpha = S.cellAlpha; pass = 1; r = r0 * 1.5; run(); } }
  ctx.restore();
}
function cellTip(x, y) {
  let best = -1, bd = 64 * devicePixelRatio * devicePixelRatio;
  for (let j = 0; j < C.hits.length; j += 3) { const d = (C.hits[j + 1] - x) ** 2 + (C.hits[j + 2] - y) ** 2; if (d < bd) { bd = d; best = C.hits[j]; } }
  if (best < 0) return '';
  const h = C.hgea[best];
  const ex = G.expr ? `<br>&nbsp;${G.gene} log2: ${G.expr[best] ? (G.expr[best] / 255 * G.vmax).toFixed(2) : 'not detected'}` : '';
  return `<br><b>cell</b> ${C.tax.subclass[C.sub[best]].name}${ex}<br>&nbsp;HGEA domain at the cell: ${h ? DOMNAME[h] : '—'}`
    + `<br>&nbsp;AP ${(C.pos[best * 3] / 1000).toFixed(2)} mm (${Math.round(C.pos[best * 3] - S.pos.coronal * 25)} µm from this slice)`;
}
let cellPts = null, cutPts = null;
function cutCells(k, zc) {   // the slab's cells drawn on the coronal 3D cut face, over the HGEA drawing
  if (cutPts) { scene.remove(cutPts); cutPts.geometry.dispose(); cutPts = null; }
  if (!S.cells || !C.pos || S.view !== 'cut3d') return;
  const [a, b] = apRange(k * 25 - S.slab, k * 25 + S.slab), pos = [], col = [];
  for (let j = a; j < b; j++) {
    const i = C.order[j]; if (!C.show[i]) continue;
    const p = toScene(k, C.pos[i * 3 + 1] / 25, C.pos[i * 3 + 2] / 25);
    pos.push(p[0], p[1], zc + 0.006);
    col.push(C.rgb[i * 3] / 255, C.rgb[i * 3 + 1] / 255, C.rgb[i * 3 + 2] / 255);
  }
  const g = new THREE.BufferGeometry();
  g.setAttribute('position', new THREE.Float32BufferAttribute(pos, 3)); g.setAttribute('color', new THREE.Float32BufferAttribute(col, 3));
  cutPts = new THREE.Points(g, new THREE.PointsMaterial({ size: 0.032, vertexColors: true, opacity: S.cellAlpha, transparent: S.cellAlpha < 1 }));
  cutPts.renderOrder = 2;
  scene.add(cutPts);
}
function build3DCells() {
  if (!scene) return;
  if (cellPts) { scene.remove(cellPts); cellPts.geometry.dispose(); cellPts = null; }
  if (!S.cells || !S.cells3d || !C.pos) { need3D = true; return; }
  const pos = [], col = [];
  for (let i = 0; i < C.n; i++) {
    if (!C.show[i]) continue;
    pos.push(...toScene(C.pos[i * 3] / 25, C.pos[i * 3 + 1] / 25, C.pos[i * 3 + 2] / 25));
    col.push(C.rgb[i * 3] / 255, C.rgb[i * 3 + 1] / 255, C.rgb[i * 3 + 2] / 255);
  }
  const g = new THREE.BufferGeometry();
  g.setAttribute('position', new THREE.Float32BufferAttribute(pos, 3)); g.setAttribute('color', new THREE.Float32BufferAttribute(col, 3));
  cellPts = new THREE.Points(g, new THREE.PointsMaterial({ size: 0.02, vertexColors: true, opacity: S.cellAlpha, transparent: S.cellAlpha < 1, clippingPlanes: S.view === 'cut3d' ? [clip] : [] }));
  scene.add(cellPts); need3D = true;
}

boot().catch((e) => { $('loading').textContent = `could not load: ${e.message}`; console.error(e); });

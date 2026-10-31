// CCF Gene Atlas viewer: Blue Brain ISH (CCFv3BBP) x MERFISH (Allen Brain Cell Atlas), one brain at a time.
import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { CONFIG } from './config.js?v=20261031092245';
import * as IO from './io.js?v=20261031092245';
import { lut, NAMES as CMAPS, hexToRgb, gradientCss } from './colormaps.js?v=20261031092245';
import { initHgea, selectDomain, loadSignature, hgeaBox } from './hgea.js?v=20261031092245';

const $ = (id) => document.getElementById(id);
const [NA, ND, NL] = IO.GRID;
const LEVELS = ['class', 'subclass', 'supertype', 'cluster'];
const ishOnly = () => S.mode === 'ish' || S.mode === 'hgea';   // modes without MERFISH cells

// ───────────────────────── state (mirrored in the URL) ─────────────────────────
const S = {
  mode: 'ish', brain: null, labelsrc: 'ccf', quality: 0, psize: 14,
  gene: null, cmap: 'magma', colorby: 'type', level: 'subclass', typeq: '',
  region: 0, onlyRegion: false, showMesh: true, view: '3d', slab: 200, shell: true,
  pos: { coronal: 283, sagittal: 228, horizontal: 160 },
  ishlo: 20, ishhi: 100, ishthr: 55, merlo: 0, merhi: 60, cellshow: 'all', hpfOnly: true,
};
const URLKEYS = ['mode', 'brain', 'labelsrc', 'gene', 'cmap', 'colorby', 'level', 'view'];
function readURL() {
  const q = new URLSearchParams(location.search);
  for (const k of URLKEYS) if (q.has(k)) S[k] = q.get(k);
  if (q.has('region')) S.regionAcr = q.get('region');
}
function writeURL() {
  const q = new URLSearchParams(location.search);
  for (const k of URLKEYS) { if (S[k] != null && S[k] !== '') q.set(k, S[k]); else q.delete(k); }
  if (S.region) q.set('region', ONTO[S.region].acronym); else q.delete('region');
  history.replaceState(null, '', `${location.pathname}?${q}`);
}

// ───────────────────────── data ─────────────────────────
let REF, ONTO, TAX, ISH, DS, ANN = null, REFVOL = null, annPromise = null, refPromise = null;
const B = {};                      // brain key -> loaded arrays
const ishCache = new Map(), merCache = new Map();
let GENES = new Map();             // symbol -> {ish: index|-1, mer: {key: index}}
let PARENT, DEPTH, CL2 = {};       // ontology parent/depth arrays; cluster -> level index
let LVCOL = {};                    // level -> Uint8Array rgb

async function boot() {
  readURL();
  [REF, TAX, ISH, DS] = await Promise.all([IO.json(IO.data('ref/manifest.json')), IO.json(IO.data('merfish/taxonomy.json')),
    IO.json(IO.data('ish/genes.json')), IO.json(IO.data('datasets.json'))]);
  const onto = await IO.json(IO.data('ref/ontology.json'));
  ONTO = []; onto.forEach((o) => { ONTO[o.i] = o; });
  ONTO[0] = { i: 0, id: 0, acronym: '—', name: 'outside the brain', parent: 0, depth: -1, rgb: [0, 0, 0], vox: 0 };
  PARENT = new Int32Array(ONTO.length); DEPTH = new Int32Array(ONTO.length);
  ONTO.forEach((o) => { PARENT[o.i] = o.parent; DEPTH[o.i] = o.depth; });
  buildTaxonomy();
  const mans = await Promise.all(DS.merfish.map((d) => IO.json(IO.data(`merfish/${d.key}/manifest.json`))));
  DS.merfish.forEach((d, k) => { B[d.key] = { man: mans[k] }; });
  buildGeneIndex();
  if (!S.brain || !B[S.brain]) S.brain = DS.merfish[0].key;
  if (S.regionAcr) { const r = ONTO.find((o) => o && o.acronym === S.regionAcr); if (r) S.region = r.i; }
  initUI();
  init3D();
  annPromise = IO.array(IO.data('ref/annotation.u16.gz'), 'uint16').then((a) => { ANN = a; requestDraw(); });
  await loadBrain(S.brain);
  if (!S.gene && !new URLSearchParams(location.search).has('nogene') && S.mode !== 'hgea') S.gene = 'Pvalb';   // open with a gene in both atlases
  if (S.mode === 'hgea') await enterHgea();
  const urlGenes = (S.gene || '').split(',').filter(Boolean);
  for (const sym of urlGenes) {
    const m = sym.match(/^(sig|con):(.+)$/);
    if (m) await showSignature(m[2], m[1] === 'sig' ? 'signature' : 'contrast');
  }
  await addGenes(urlGenes.filter((sym) => !/^(sig|con):/.test(sym)));
  if (S.region) selectRegion(S.region, false);
  setView(S.view);
  $('loading').style.opacity = 0;
  requestDraw();
}

function buildTaxonomy() {
  const nC = TAX.cluster.length;
  CL2.cluster = new Uint16Array(nC).map((_, i) => i);
  CL2.supertype = new Uint16Array(nC).map((_, i) => TAX.cluster[i].parent);
  CL2.subclass = CL2.supertype.map((st) => TAX.supertype[st].parent);
  CL2.class = CL2.subclass.map((sc) => TAX.subclass[sc].parent);
  for (const lv of LEVELS) {
    const c = new Uint8Array(TAX[lv].length * 3);
    TAX[lv].forEach((e, i) => { const [r, g, b] = hexToRgb(e.color || '#888888'); c[i * 3] = r; c[i * 3 + 1] = g; c[i * 3 + 2] = b; });
    LVCOL[lv] = c;
  }
}

const isControl = (sym) => /^Blank-/.test(sym);   // Allen 638850 panel: 50 blank control probes, not genes

function buildGeneIndex() {
  GENES = new Map();
  ISH.genes.forEach((g, i) => GENES.set(g.gene, { ish: i, mer: {} }));
  for (const d of DS.merfish) {
    (B[d.key].man.genes || []).forEach((g, i) => {
      if (isControl(g.symbol)) return;
      if (!GENES.has(g.symbol)) GENES.set(g.symbol, { ish: -1, mer: {} });
      GENES.get(g.symbol).mer[d.key] = i;
    });
  }
}

async function loadBrain(key) {
  const b = B[key];
  if (!b.pos) {
    status(`loading ${b.man.label}…`);
    const f = b.man.files, base = IO.data(`merfish/${key}/`);
    const [um, cluster, rccf, rbbp, quality, hq] = await Promise.all([
      IO.array(base + f.pos.file, 'uint16'), IO.array(base + f.cluster.file, 'uint16'),
      IO.array(base + f.region_ccf.file, 'uint16'), IO.array(base + f.region_bbp.file, 'uint16'),
      IO.array(base + f.quality.file, 'uint8'), f.high_quality ? IO.array(base + f.high_quality.file, 'uint8') : null]);
    const n = cluster.length, pos = new Float32Array(n * 3);
    for (let i = 0; i < n; i++) {
      const s = IO.toScene(um[i * 3], um[i * 3 + 1], um[i * 3 + 2]);
      pos[i * 3] = s[0]; pos[i * 3 + 1] = s[1]; pos[i * 3 + 2] = s[2];
    }
    let side = 0;
    for (let i = 0; i < n; i++) if (um[i * 3 + 2] < IO.CENTER[2]) side++;
    b.oneSide = Math.max(side, n - side) / n;
    Object.assign(b, { um, pos, cluster, region_ccf: rccf, region_bbp: rbbp, quality, hq, n, colors: new Uint8Array(n * 3) });
    b.sections = await sectionPlanes(b, base + f.section.file);
  }
  S.brain = key;
  $('lblR').textContent = `MERFISH · ${b.man.dataset}`;
  setupCells();
  for (const g of SEL) await loadMer(g);
  recolorCells();
  refreshProfile();
  renderGeneRegions();
  status('');
}

async function sectionPlanes(b, url) {
  // each original section is a plane: coronal brains -> fixed AP, sagittal brains -> fixed |LR - midline|
  const sec = await IO.array(url, 'uint16'), ns = b.man.sections.length, n = b.n, um = b.um;
  const sum = new Float64Array(ns * 2), sq = new Float64Array(ns * 2), cnt = new Float64Array(ns), mid = IO.CENTER[2];
  for (let i = 0; i < n; i++) {
    const s = sec[i], ap = um[i * 3], lr = Math.abs(um[i * 3 + 2] - mid);
    sum[s * 2] += ap; sum[s * 2 + 1] += lr; sq[s * 2] += ap * ap; sq[s * 2 + 1] += lr * lr; cnt[s]++;
  }
  const sd = (s, k) => Math.sqrt(Math.max(0, sq[s * 2 + k] / cnt[s] - (sum[s * 2 + k] / cnt[s]) ** 2));
  const sdAP = [], sdLR = [];
  for (let s = 0; s < ns; s++) if (cnt[s]) { sdAP.push(sd(s, 0)); sdLR.push(sd(s, 1)); }
  const med = (a) => a.sort((x, y) => x - y)[a.length >> 1];
  const orient = med(sdAP) < med(sdLR) ? 'coronal' : 'sagittal';
  const pos = [];
  for (let s = 0; s < ns; s++) if (cnt[s]) pos.push(sum[s * 2 + (orient === 'coronal' ? 0 : 1)] / cnt[s]);
  pos.sort((x, y) => x - y);
  const gaps = pos.slice(1).map((x, k) => x - pos[k]);
  return { orient, pos, spacing: gaps.length ? med([...gaps]) : 0, maxGap: gaps.length ? Math.max(...gaps) : 0 };
}

function sectionTicks() {
  const el = $('ticks'); el.innerHTML = '';
  const b = B[S.brain]; if (!b || !b.sections || S.view === '3d') return;
  const sp = b.sections, P = PLANES[S.view];
  if (sp.orient !== S.view) return;          // ticks only where the brain's sections are planes of this view
  for (const x of sp.pos) {
    const k = S.view === 'coronal' ? x / 25 : (IO.CENTER[2] + x) / 25;
    const i = document.createElement('i'); i.style.left = `${k / (P.n - 1) * 100}%`; el.appendChild(i);
  }
}

function jumpSection(dir) {
  const b = B[S.brain]; if (!b || !b.sections || S.view === '3d') return;
  const sp = b.sections;
  if (sp.orient !== S.view) { status(`${b.man.dataset} was cut ${sp.orient}ly; switch to the ${sp.orient} view to step through its sections`); return; }
  const cur = S.pos[S.view] * 25, toUm = (x) => (S.view === 'coronal' ? x : IO.CENTER[2] + x);
  const list = sp.pos.map(toUm);
  const next = dir > 0 ? list.find((x) => x > cur + 12.5) : [...list].reverse().find((x) => x < cur - 12.5);
  if (next == null) return;
  S.pos[S.view] = Math.round(next / 25); $('slicepos').value = S.pos[S.view]; requestDraw();
}

async function loadRegionCluster(key, src) {
  const b = B[key], k = `rc_${src}`;
  if (b[k]) return b[k];
  const f = b.man.files[`region_cluster_${src}`], base = IO.data(`merfish/${key}/`);
  const [region, cluster, count] = await Promise.all([IO.array(base + f.region.file, 'uint16'),
    IO.array(base + f.cluster.file, 'uint16'), IO.array(base + f.count.file, 'uint32')]);
  return (b[k] = { region, cluster, count });
}

async function loadRegionGene(key) {
  const b = B[key];
  if (b.rmean) return b;
  const f = b.man.files, base = IO.data(`merfish/${key}/`);
  if (!f.region_gene_mean) return b;
  [b.rmean, b.rfrac] = await Promise.all([IO.array(base + f.region_gene_mean.file, 'float16'), IO.array(base + f.region_gene_frac.file, 'float16')]);
  return b;
}

let ISHR = null;
async function loadIshRegions() {
  if (ISHR) return ISHR;
  const [mean, frac] = await Promise.all([IO.array(IO.data('ish/' + ISH.files.region_mean.file), 'float16'),
    IO.array(IO.data('ish/' + ISH.files.region_frac.file), 'float16')]);
  return (ISHR = { mean, frac, G: ISH.genes.length });
}

// ───────────────────────── gene selection: one or several genes, each a colour channel ─────────────────────────
const CHANNELS = ['#ff3dcf', '#45f06a', '#3db8ff', '#ffc93d'];   // magenta, green, cyan, amber (additive)
const ADDITIVE_MAX = 4;   // up to here colours add up; beyond it they would saturate, so each spot takes its strongest gene's colour
const MAXG = 50;
const LOAD_POOL = 4;      // volumes fetched at once when adding a list of genes
const KEEP_REMOVED_ISH = 8, KEEP_REMOVED_MER = 50;   // removed genes kept in memory for quick re-adding
const PALETTE = [...CHANNELS];
for (let k = 0; PALETTE.length < MAXG; k++) {        // golden-angle hues, three lightness steps
  const h = (205 + k * 137.508) % 360, l = [0.62, 0.72, 0.52][k % 3], a = 0.85 * Math.min(l, 1 - l);
  const f = (n) => { const m = (n + h / 30) % 12; return Math.round(255 * (l - a * Math.max(-1, Math.min(m - 3, 9 - m, 1)))).toString(16).padStart(2, '0'); };
  PALETTE.push(`#${f(0)}${f(8)}${f(4)}`);
}
const nextColor = () => { const used = new Set(SEL.map((x) => x.color)); return PALETTE.find((c) => !used.has(c)); };
let SEL = [];          // selected genes: {symbol, ishIdx, ish, mer, merIdx, color, rgb}
let GENE = null;       // the active gene (availability chips, single-gene colour map)
const multi = () => SEL.length > 1;
const shown = () => SEL.filter((g) => !g.hidden);     // genes currently toggled on
let BATCH = false;     // adding a list of genes: redraw once, at the end

// ISH volumes are held packed: only the voxels inside the brain (32.7 M of 82.6 M), in grid order, one shared map.
// Nothing outside the brain is ever displayed, so this loses nothing visible: 33 MB per gene instead of 82 MB.
let BM = null;
function brainMap() {
  if (BM) return BM;
  const N = ANN.length, words = (N + 31) >>> 5, bits = new Uint32Array(words), cum = new Uint32Array(words), starts = [], lens = [];
  let n = 0, run = -1;
  for (let i = 0; i < N; i++) {
    if (!(i & 31)) cum[i >>> 5] = n;
    if (ANN[i]) { bits[i >>> 5] |= 1 << (i & 31); n++; if (run < 0) run = i; }
    else if (run >= 0) { starts.push(run); lens.push(i - run); run = -1; }
  }
  if (run >= 0) { starts.push(run); lens.push(N - run); }
  return (BM = { n, bits, cum, starts: Uint32Array.from(starts), lens: Uint32Array.from(lens) });
}
function packVol(full, out = null) {
  const { n, starts, lens } = brainMap();
  if (!out) out = new Uint8Array(n);
  for (let r = 0, p = 0; r < starts.length; r++) { out.set(full.subarray(starts[r], starts[r] + lens[r]), p); p += lens[r]; }
  return out;
}
function packedIndex(i) {   // grid index -> position in a packed volume (-1 outside the brain)
  const w = i >>> 5, b = BM.bits[w], m = 1 << (i & 31);
  if (!(b & m)) return -1;
  let x = b & (m - 1);
  x -= (x >>> 1) & 0x55555555; x = (x & 0x33333333) + ((x >>> 2) & 0x33333333);
  return BM.cum[w] + ((((x + (x >>> 4)) & 0x0f0f0f0f) * 0x01010101) >>> 24);
}
function pruneCaches() {    // selected genes stay; only the most recently used removed genes are kept
  const trim = (cache, sel, keep) => { const old = [...cache.keys()].filter((k) => !sel.has(k)); old.slice(0, Math.max(0, old.length - keep)).forEach((k) => cache.delete(k)); };
  trim(ishCache, new Set(SEL.map((g) => g.symbol)), KEEP_REMOVED_ISH);
  trim(merCache, new Set(SEL.map((g) => `${S.brain}|${g.symbol}`)), KEEP_REMOVED_MER);
}

async function addGenes(symbols) {   // several genes in one go: limited parallel downloads, one redraw
  const list = [...new Set(symbols)];
  if (!list.length) return;
  const unknown = list.filter((x) => !GENES.has(x)), known = list.filter((x) => GENES.has(x)).slice(0, MAXG);
  BATCH = true;
  let next = 0, done = 0;
  const worker = async () => { while (next < known.length) { const sym = known[next++]; await selectGene(sym); done++; if (known.length > 1) status(`loaded ${done} of ${known.length} genes…`); } };
  try { await Promise.all(Array.from({ length: LOAD_POOL }, worker)); } finally { BATCH = false; }
  refreshGenes();
  if (await showOnPanelBrain()) return;
  if (known.length === 1 && !unknown.length) return;   // selectGene already reported it
  status(`${known.length} gene${known.length === 1 ? '' : 's'} added${unknown.length ? ` · not found: ${unknown.join(', ')}` : ''}${list.length - unknown.length > MAXG ? ` · only the first ${MAXG} kept` : ''}`);
}

async function selectGene(sym, replace = false) {
  const e = GENES.get(sym);
  if (!e) { status(`"${sym}" is not in the ISH atlas or any MERFISH panel`); return; }
  if (replace) SEL = [];
  let g = SEL.find((x) => x.symbol === sym);
  if (!g) {
    makeRoom();
    const color = nextColor();
    g = { symbol: sym, ishIdx: e.ish, ish: null, mer: null, merIdx: null, color, rgb: hexToRgb(color) };
    SEL.push(g);
  }
  GENE = g; $('gene').value = ''; syncGeneURL(); renderGeneChips(); renderAvail();
  const jobs = [];
  if (e.ish >= 0 && !g.ish) {
    jobs.push((async () => {
      if (!ishCache.has(sym)) {
        status(`streaming ${sym} ISH volume from EBRAINS (${(ISH.genes[e.ish].bytes_gz / 1e6).toFixed(1)} MB)…`);
        const vol = await IO.ishVolume(ISH.url_template, sym, progress);
        await annPromise;
        ishCache.set(sym, packVol(vol));
      }
      g.ish = ishCache.get(sym);
    })());
  }
  jobs.push(loadMer(g));
  try { await Promise.all(jobs); } catch (err) { status(`could not load ${sym}: ${err.message}`); console.error(err); }
  progress(0, 0);
  if (g.mer && S.colorby === 'type') setSeg('colorby', 'gene');   // also in ISH mode, so the MERFISH tab opens on this gene
  refreshGenes();
  const bits = [];
  if (g.ish) bits.push('ISH ✓');
  if (e.mer[S.brain] != null) bits.push(g.mer ? `MERFISH ${S.brain} ✓` : `MERFISH ${S.brain}: expression host unreachable`);
  else bits.push(`not on the ${B[S.brain].man.dataset} panel`);
  status(`${sym}: ${bits.join(' · ')}${multi() ? ` · ${SEL.length} genes selected` : ''}`);
  if (!BATCH) await showOnPanelBrain();
}

// MERFISH tab: panels differ (Allen 550 genes, Zhuang 1,122), so go to a brain that measured the active gene.
async function showOnPanelBrain() {
  if (S.mode !== 'merfish' || !GENE || GENE.pseudo) return false;
  const sym = GENE.symbol, e = GENES.get(sym);
  if (!e || e.mer[S.brain] != null) return false;
  const opts = DS.merfish.filter((d) => e.mer[d.key] != null);
  if (!opts.length) { status(`${sym}: not on any MERFISH panel · ISH only`); return false; }
  const cover = (d) => SEL.filter((g) => !g.pseudo && GENES.get(g.symbol).mer[d.key] != null).length;
  const to = opts.reduce((a, d) => (cover(d) > cover(a) ? d : a));
  const from = B[S.brain].man.dataset;
  if (S.colorby === 'type') setSeg('colorby', 'gene');
  $('brain').value = to.key;
  await loadBrain(to.key);
  writeURL(); renderAvail(); buildIshPoints(); sectionTicks(); SL.key = ''; requestDraw();
  status(`${sym} is not on the ${from} panel · switched to ${B[to.key].man.dataset}`);
  return true;
}

function makeRoom() {   // at the cap, the oldest gene/layer gives way instead of refusing silently
  while (SEL.length >= MAXG) {
    const old = SEL.shift();
    if (GENE === old) GENE = null;
    status(`up to ${MAXG} at once — removed ${old.label || old.symbol}`);
  }
}

function removeGene(sym) {
  SEL = SEL.filter((g) => g.symbol !== sym);
  pruneCaches();
  if (!GENE || GENE.symbol === sym) GENE = SEL[SEL.length - 1] || null;
  syncGeneURL(); renderGeneChips(); renderAvail(); refreshGenes();
  status(SEL.length ? `${SEL.length} gene${SEL.length > 1 ? 's' : ''} selected` : 'type a gene symbol, e.g. Prox1, Gad2, Slc17a7');
}

function toggleGene(g) {
  if (!g) return;
  g.hidden = !g.hidden;
  if (g.hidden && !shown().length) { g.hidden = false; status('at least one gene stays visible'); return; }
  renderGeneChips(); SL.key = ''; buildIshPoints(); recolorCells(); renderLegend(); requestDraw();
}

// ───────────────────────── HGEA-ISH mode ─────────────────────────
let HG = null;
const HPF_MAIN = ['CA1', 'CA2', 'CA3', 'DG', 'SUB', 'ProS', 'FC', 'IG'];      // the HGEA territory in CCF
const HPF_CONTEXT = ['ENTl', 'ENTm', 'PRE', 'POST', 'PAR', 'HATA', 'APr'];   // retrohippocampal neighbours
const hpfMesh = new Map();     // acronym -> {obj, on}
let HPFMASK = null;

function hpfRestrict() {       // Uint8Array over ontology indices (HGEA territory), or null when not restricting
  if (S.mode !== 'hgea' || !S.hpfOnly) return null;
  if (!HPFMASK) {
    HPFMASK = new Uint8Array(ONTO.length);
    for (const acr of HPF_MAIN) { const o = ONTO.find((x) => x && x.acronym === acr); if (o) { const m = regionMask(o.i); for (let i = 0; i < m.length; i++) if (m[i]) HPFMASK[i] = 1; } }
  }
  return HPFMASK;
}

async function setHpfMesh(acr, on) {
  let e = hpfMesh.get(acr);
  if (!e) { e = { obj: null, on: false }; hpfMesh.set(acr, e); }
  e.on = on;
  if (on && !e.obj) {
    const o = ONTO.find((x) => x && x.acronym === acr);
    const m = await IO.mesh(o.id);
    const [r, g, b] = o.rgb;
    e.obj = meshObject(m, new THREE.Color(r / 255, g / 255, b / 255), 0.24);
    scene.add(e.obj);
  }
  if (e.obj) e.obj.visible = on && S.mode === 'hgea';
  document.querySelectorAll(`.hmesh[data-acr="${acr}"]`).forEach((b) => b.classList.toggle('on', on));
  requestDraw();
}

function showHpfMeshes(visible) { hpfMesh.forEach((e) => { if (e.obj) e.obj.visible = visible && e.on; }); requestDraw(); }

function renderHpfControls() {
  const chip = (acr) => { const o = ONTO.find((x) => x && x.acronym === acr);
    return `<button class="hmesh" data-acr="${acr}" title="${o.name} (CCFv3BBP mesh)"><span class="sw" style="background:rgb(${o.rgb})"></span>${acr}</button>`; };
  $('hgeahpf').innerHTML = `<div class="row" style="margin-top:0"><label class="mini" style="display:flex;gap:5px;align-items:center"><input type="checkbox" id="hpfonly" ${S.hpfOnly ? 'checked' : ''}> only hippocampus (CA1–3, DG, SUB, ProS, FC, IG)</label></div>
    <div class="mini" style="margin-top:6px">CCF meshes · click to toggle</div><div class="hrow">${HPF_MAIN.map(chip).join('')}</div>
    <div class="mini" style="margin-top:4px">retrohippocampal context</div><div class="hrow">${HPF_CONTEXT.map(chip).join('')}</div>
    <div class="row"><button class="btn sm" id="hpfall">all hippocampus meshes</button><button class="btn sm" id="hpfctx">+ context</button><button class="btn sm" id="hpfnone">none</button></div>`;
  $('hgeahpf').querySelectorAll('.hmesh').forEach((b) => b.addEventListener('click', () => setHpfMesh(b.dataset.acr, !(hpfMesh.get(b.dataset.acr) || {}).on)));
  $('hpfonly').addEventListener('change', () => { S.hpfOnly = $('hpfonly').checked; SL.key = ''; buildIshPoints(); requestDraw(); });
  $('hpfall').addEventListener('click', () => HPF_MAIN.forEach((a) => setHpfMesh(a, true)));
  $('hpfctx').addEventListener('click', () => HPF_CONTEXT.forEach((a) => setHpfMesh(a, true)));
  $('hpfnone').addEventListener('click', () => [...HPF_MAIN, ...HPF_CONTEXT].forEach((a) => setHpfMesh(a, false)));
}

async function enterHgea() {
  $('hgeasec').style.display = 'block';
  if (!$('hgeahpf').innerHTML) renderHpfControls();
  showHpfMeshes(true);
  if (!HG) {
    HG = await initHgea({ selectGene: (sym) => selectGene(sym), showSignature, status, upsertLayer, removeLayer: (key) => { const g = SEL.find((x) => x.symbol === key); if (g) removeGene(key); },
      setDomain: (c) => { S.hdom = c; const q = new URLSearchParams(location.search); q.set('hdom', c); history.replaceState(null, '', `${location.pathname}?${q}`); } });
    const hd = new URLSearchParams(location.search).get('hdom');
    selectDomain(hd && HG.columns.includes(hd) ? hd : 'DGd');
  }
  focusHippocampus();
}

function focusHippocampus() {
  const box = hgeaBox(); if (!box) return;
  const c = [(box[0][0] + box[1][0]) / 2, (box[0][1] + box[1][1]) / 2, (box[0][2] + box[1][2]) / 2];
  const t = IO.toScene(c[0] * 25, c[1] * 25, c[2] * 25);
  controls.target.set(t[0], t[1], t[2]);
  camera.position.set(t[0] + 7.5, t[1] + 5, t[2] + 7.5);
  S.pos.coronal = Math.round(c[0]); S.pos.sagittal = Math.round(box[0][2] + (box[1][2] - box[0][2]) * 0.22); S.pos.horizontal = Math.round(c[1]);
  SL.focus = c; SL.zoom = 1.3;
  if (S.view !== '3d') { $('slicepos').value = S.pos[S.view]; SL.key = ''; }
  requestDraw();
}

async function upsertLayer({ key, label, vol, p995, tip, legend, lo = null, thr3d = null }) {   // a derived ISH-space layer (e.g. all HGEA genes)
  await annPromise;
  let g = SEL.find((x) => x.symbol === key);
  if (!g) {
    makeRoom();
    const color = nextColor();
    g = { symbol: key, pseudo: true, kind: 'layer', mer: null, color, rgb: hexToRgb(color) };
    SEL.push(g);
  }
  Object.assign(g, { label, ish: packVol(vol, g.ish), p995, tip, legend, lo, thr3d });
  GENE = g; syncGeneURL(); renderGeneChips(); renderAvail(); refreshGenes();
  status(label);
}

async function showSignature(col, kind) {
  const key = `${kind === 'signature' ? 'sig' : 'con'}:${col}`;
  let g = SEL.find((x) => x.symbol === key);
  if (!g) {
    makeRoom();
    status(`loading the ${col} ${kind}…`);
    const sg = await loadSignature(col, kind);
    await annPromise;
    const color = nextColor();
    g = { symbol: key, label: `${col} ${kind}`, pseudo: true, kind, col, ish: packVol(sg.vol), p995: sg.p995, mer: null, color, rgb: hexToRgb(color) };
    SEL.push(g);
  }
  GENE = g; syncGeneURL(); renderGeneChips(); renderAvail(); refreshGenes();
  status(`${g.label}: Blue Brain ISH of the genes annotated in ${col} (Bienkowski et al. 2018, Supplementary Table 2)`);
}

function refreshGenes() { if (BATCH) { renderGeneChips(); return; } SL.key = ''; buildIshPoints(); recolorCells(); refreshProfile(); renderLegend(); renderGeneRegions(); requestDraw(); }
function syncGeneURL() { S.gene = SEL.filter((g) => g.symbol !== 'hgea:all').map((g) => g.symbol).join(','); writeURL(); }

function renderGeneChips() {
  const el = $('genechips'); el.innerHTML = '';
  for (const g of SEL) {
    const c = document.createElement('span');
    c.className = 'gchip' + (g === GENE ? ' active' : '') + (g.hidden ? ' hidden' : '');
    c.innerHTML = `<span class="dot" style="background:${multi() ? g.color : 'var(--accent)'}"></span>${g.label || g.symbol}<span class="x" title="remove">×</span>`;
    c.addEventListener('click', (ev) => {
      if (ev.target.classList.contains('x')) { removeGene(g.symbol); return; }
      GENE = g; renderGeneChips(); renderAvail(); renderLegend(); renderGeneRegions(); if (!multi()) { recolorCells(); SL.key = ''; buildIshPoints(); }
    });
    el.appendChild(c);
  }
  if (SEL.length > 1) {
    const clr = document.createElement('a'); clr.href = '#'; clr.className = 'mini'; clr.textContent = 'clear all';
    clr.addEventListener('click', (ev) => { ev.preventDefault(); [...SEL].forEach((g) => removeGene(g.symbol)); });
    el.appendChild(clr);
  }
  $('gene').placeholder = SEL.length ? `Add a gene, or paste a list (up to ${MAXG})…` : 'Search a gene, or paste a list (ISH ∪ MERFISH panels)…';
}

async function loadMer(g) {
  if (g.pseudo) return;
  const e = GENES.get(g.symbol);
  const gi = e.mer[S.brain];
  g.mer = null; g.merIdx = gi;
  if (gi == null) return;
  const brain = S.brain, key = `${brain}|${g.symbol}`;
  if (!merCache.has(key)) {
    const url = `${CONFIG.EXPR}merfish/${brain}/expr/${encodeURIComponent(g.symbol)}.u8.gz`;
    try { merCache.set(key, await IO.array(url, 'uint8', progress)); } catch (err) { console.warn(err); return; }
  }
  if (S.brain !== brain) return;   // the brain changed during the download: this array belongs to other cells
  g.mer = merCache.get(key);
}

const LOW_ISH = 10;        // ISH p99.5 at or below this (of 255) is flagged: such genes agree less with MERFISH, see the accuracy page
const ishP995 = (g) => (g && !g.pseudo && g.ishIdx >= 0 ? ISH.genes[g.ishIdx].p995 : null);
const lowIsh = (g) => { const p = ishP995(g); return p != null && p <= LOW_ISH; };
const lowIshNote = (g) => `Low ISH signal: p99.5 = ${ishP995(g)} of 255. Colours are scaled to this gene's own range, so faint signal can look strong.`;

function renderAvail() {
  const el = $('avail'); el.innerHTML = '';
  if (GENE && GENE.pseudo) { el.innerHTML = `<span class="chip ish on">HGEA ${GENE.kind}</span><span class="mini">derived from the paper's marker genes (ISH)</span>`; return; }
  const e = GENE && GENES.get(GENE.symbol); if (!e) return;
  const chip = (cls, on, txt, title) => { const c = document.createElement('span'); c.className = `chip ${cls} ${on ? 'on' : 'off'}`; c.textContent = txt; c.title = title; el.appendChild(c); };
  chip('ish', e.ish >= 0, 'ISH', e.ish >= 0 ? 'Blue Brain ISH volume (25 µm)' : 'no ISH volume');
  for (const d of DS.merfish) chip('mer', e.mer[d.key] != null, d.short, `${d.label}: ${e.mer[d.key] != null ? 'on the panel' : 'not on the panel'}`);
  if (lowIsh(GENE)) { const m = document.createElement('div'); m.className = 'mini'; m.style.cssText = 'flex-basis:100%;margin-top:4px;color:#a8591a'; m.textContent = lowIshNote(GENE); el.appendChild(m); }
}

// ───────────────────────── 3D ─────────────────────────
let renderer, scene, camera, controls, shellMesh = null, regionMesh = null, cellPts = null, bgPts = null, ishPts = null, needDraw = true;
function init3D() {
  const canvas = $('gl');
  renderer = new THREE.WebGLRenderer({ canvas, antialias: true, alpha: true });
  renderer.setPixelRatio(Math.min(devicePixelRatio, 2));
  scene = new THREE.Scene();
  camera = new THREE.PerspectiveCamera(35, 1, 0.01, 200);
  camera.position.set(14, 7, 12);
  controls = new OrbitControls(camera, canvas);
  controls.enableDamping = true; controls.dampingFactor = 0.12;
  controls.addEventListener('change', requestDraw);
  camera.layers.enableAll();
  scene.add(new THREE.AmbientLight(0xffffff, 0.75));
  const dl = new THREE.DirectionalLight(0xffffff, 0.8); dl.position.set(5, 10, 7); scene.add(dl);
  IO.mesh(997).then((m) => { shellMesh = meshObject(m, 0x9fb7c9, 0.08); shellMesh.visible = S.shell; scene.add(shellMesh); requestDraw(); });
  new ResizeObserver(resize).observe($('stage'));
  resize();
  (function loop() { requestAnimationFrame(loop); controls.update(); if (needDraw) { needDraw = false; draw(); } })();
}
function requestDraw() { needDraw = true; }
function resize() {
  const st = $('stage'), w = st.clientWidth, h = st.clientHeight;
  renderer.setSize(w, h, false);
  const sc = $('slice'); sc.width = w * devicePixelRatio; sc.height = h * devicePixelRatio;
  requestDraw();
}

function meshObject(m, color, opacity) {
  const n = m.um.length / 3, pos = new Float32Array(m.um.length);
  for (let i = 0; i < n; i++) { const s = IO.toScene(m.um[i * 3], m.um[i * 3 + 1], m.um[i * 3 + 2]); pos.set(s, i * 3); }
  const g = new THREE.BufferGeometry();
  g.setAttribute('position', new THREE.BufferAttribute(pos, 3));
  g.setIndex(new THREE.BufferAttribute(m.index, 1));
  g.computeVertexNormals();
  const mat = new THREE.MeshLambertMaterial({ color, transparent: true, opacity, depthWrite: false, side: THREE.DoubleSide });
  const o = new THREE.Mesh(g, mat);
  o.layers.enableAll();
  return o;
}

function setupCells() {
  const b = B[S.brain];
  if (cellPts) { scene.remove(cellPts); cellPts.geometry.dispose(); }
  if (bgPts) { scene.remove(bgPts); bgPts.geometry.dispose(); }
  const g = new THREE.BufferGeometry();
  g.setAttribute('position', new THREE.BufferAttribute(b.pos, 3));
  g.setAttribute('color', new THREE.BufferAttribute(b.colors, 3, true));
  g.setIndex(new THREE.BufferAttribute(new Uint32Array(b.n), 1));
  cellPts = new THREE.Points(g, new THREE.PointsMaterial({ size: S.psize / 1000, vertexColors: true, sizeAttenuation: true }));
  cellPts.layers.set(2);
  scene.add(cellPts);
  // MERFISH tab, cells by gene: cells without signal as a faint haze, so they cannot hide a sparse gene's cells
  const bg = new THREE.BufferGeometry();
  bg.setAttribute('position', g.attributes.position);
  bgPts = new THREE.Points(bg, new THREE.PointsMaterial({ size: S.psize / 1000, color: 0x8fa3b0, sizeAttenuation: true, transparent: true, opacity: 0.035, depthWrite: false }));
  bgPts.layers.set(2); bgPts.visible = false;
  scene.add(bgPts);
}

function regionMask(r) {           // Uint8Array over ontology indices: 1 if index is r or a descendant
  const m = new Uint8Array(ONTO.length);
  for (let i = 1; i < ONTO.length; i++) { let a = i; while (a) { if (a === r) { m[i] = 1; break; } a = PARENT[a]; } }
  return m;
}

function typeMatch() {             // Uint8Array over clusters, or null (no filter)
  const q = S.typeq.trim().toLowerCase(); if (!q) return null;
  const lvIdx = CL2[S.level], names = TAX[S.level], m = new Uint8Array(TAX.cluster.length);
  const hit = names.map((e) => e.name.toLowerCase().includes(q));
  for (let c = 0; c < m.length; c++) m[c] = hit[lvIdx[c]] ? 1 : 0;
  return m;
}

function exprFilter() {              // cells shown by signal: all | signal for any selected gene | signal for all of them
  if (S.cellshow === 'all' || S.colorby !== 'gene') return null;
  const arrs = shown().map((g) => g.mer).filter(Boolean);
  if (!arrs.length) return null;
  if (S.cellshow === 'any') return (i) => { for (const a of arrs) if (a[i]) return true; return false; };
  return (i) => { for (const a of arrs) if (!a[i]) return false; return true; };
}

function cellVisible() {           // predicate over cell index, shared by 3D and slices
  const b = B[S.brain], reg = S.labelsrc === 'ccf' ? b.region_ccf : b.region_bbp;
  const rm = S.onlyRegion && S.region ? regionMask(S.region) : null;
  const tm = typeMatch(), qmin = S.quality * 255, ef = exprFilter();
  return (i) => (!rm || rm[reg[i]]) && (!tm || tm[b.cluster[i]]) && b.quality[i] >= qmin && (!ef || ef(i));
}

const LOW_VISIBLE = 64;    // MERFISH tab: lowest colour-map step for a cell with signal
function recolorCells() {
  const b = B[S.brain]; if (!b || !b.pos) return;
  const c = b.colors, n = b.n, gs = shown().filter((g) => g.mer);
  const lo = S.merlo / 100 * 255, hi = Math.max(lo + 1, S.merhi / 100 * 255);
  if (S.colorby === 'gene' && gs.length && !multi()) {
    const L = lut(S.cmap), v = gs[0].mer, floor = S.mode === 'merfish' ? LOW_VISIBLE : 0;
    for (let i = 0; i < n; i++) {
      if (v[i] === 0) { c[i * 3] = 38; c[i * 3 + 1] = 44; c[i * 3 + 2] = 50; continue; }
      const k = v[i] <= lo ? 0 : Math.max(floor, Math.min(255, Math.round((v[i] - lo) / (hi - lo) * 255))) * 3;   // any detection above the range stays visible
      c[i * 3] = L[k]; c[i * 3 + 1] = L[k + 1]; c[i * 3 + 2] = L[k + 2];
    }
  } else if (S.colorby === 'gene' && gs.length > ADDITIVE_MAX) {
    // many genes: each cell takes the colour of its strongest gene (each gene relative to its own maximum)
    const best = new Uint8Array(n), who = new Uint8Array(n);
    gs.forEach((g, j) => { const v = g.mer; for (let i = 0; i < n; i++) if (v[i] > best[i]) { best[i] = v[i]; who[i] = j; } });
    for (let i = 0; i < n; i++) {
      if (!best[i]) { c[i * 3] = 38; c[i * 3 + 1] = 44; c[i * 3 + 2] = 50; continue; }
      const f = Math.max(0.25, Math.min(1, (best[i] - lo) / (hi - lo))), rgb = gs[who[i]].rgb;
      c[i * 3] = rgb[0] * f; c[i * 3 + 1] = rgb[1] * f; c[i * 3 + 2] = rgb[2] * f;
    }
  } else if (S.colorby === 'gene' && gs.length) {
    // additive channels: each gene adds its colour scaled by its expression; overlap mixes (magenta + green -> white)
    for (let i = 0; i < n; i++) {
      let r = 0, gg = 0, bb = 0, any = false;
      for (const g of gs) {
        const v = g.mer[i]; if (!v) continue;
        any = true; const f = Math.max(0.12, Math.min(1, (v - lo) / (hi - lo)));
        r += g.rgb[0] * f; gg += g.rgb[1] * f; bb += g.rgb[2] * f;
      }
      if (!any) { c[i * 3] = 38; c[i * 3 + 1] = 44; c[i * 3 + 2] = 50; continue; }
      c[i * 3] = Math.min(255, r); c[i * 3 + 1] = Math.min(255, gg); c[i * 3 + 2] = Math.min(255, bb);
    }
  } else {
    const lv = CL2[S.level], col = LVCOL[S.level];
    for (let i = 0; i < n; i++) { const k = lv[b.cluster[i]] * 3; c[i * 3] = col[k]; c[i * 3 + 1] = col[k + 1]; c[i * 3 + 2] = col[k + 2]; }
  }
  const vis = cellVisible(), idx = new Uint32Array(n), split = S.mode === 'merfish' && S.colorby === 'gene' && gs.length > 0;
  let m = 0, nb = 0, rest = null;
  if (split) {           // cells with signal solid, the others behind them as a haze
    const arrs = gs.map((g) => g.mer), expr = new Uint8Array(n);
    for (const a of arrs) for (let i = 0; i < n; i++) if (a[i] > lo) expr[i] = 1;
    rest = new Uint32Array(n);
    for (let i = 0; i < n; i++) if (vis(i)) { if (expr[i]) idx[m++] = i; else rest[nb++] = i; }
  } else for (let i = 0; i < n; i++) if (vis(i)) idx[m++] = i;
  if (cellPts) {
    cellPts.geometry.attributes.color.needsUpdate = true;
    cellPts.geometry.setIndex(new THREE.BufferAttribute(idx.subarray(0, m), 1));
    cellPts.material.size = S.psize / 1000;
    bgPts.visible = split; bgPts.material.size = S.psize / 1000;
    if (split) bgPts.geometry.setIndex(new THREE.BufferAttribute(rest.subarray(0, nb), 1));
  }
  b.nVisible = m + nb; b.nDetected = split ? m : null;
  hud(); renderLegend(); requestDraw();
}

function ishRange(g = GENE) {
  const top = Math.max(1, g.p995 != null ? g.p995 : ISH.genes[g.ishIdx].p995);
  if (g.lo != null) return [g.lo, Math.max(g.lo + 1, top), top];   // derived layer with its own floor (e.g. gene counts)
  return [S.ishlo / 100 * top, Math.max(S.ishlo / 100 * top + 1, S.ishhi / 100 * top), top];
}

function ishShade(p, gs, rngs, base) {   // colour of packed voxel p for the selected ISH genes, blended over base [r,g,b]; null = below range
  if (gs.length === 1 && !multi()) {
    const x = gs[0].ish[p], [lo, hi] = rngs[0]; if (x <= lo) return null;
    const f = Math.min(1, (x - lo) / (hi - lo)), k = Math.round(f * 255) * 3, L = lut(S.cmap), al = 0.25 + 0.75 * f;
    return [base[0] * (1 - al) + L[k] * al, base[1] * (1 - al) + L[k + 1] * al, base[2] * (1 - al) + L[k + 2] * al];
  }
  if (gs.length > ADDITIVE_MAX) {       // the strongest gene here, each relative to its own range
    let best = 0, bj = -1;
    for (let j = 0; j < gs.length; j++) {
      const x = gs[j].ish[p], lo = rngs[j][0]; if (x <= lo) continue;
      const f = Math.min(1, (x - lo) / (rngs[j][1] - lo)); if (f > best) { best = f; bj = j; }
    }
    if (bj < 0) return null;
    const al = 0.25 + 0.75 * best, c = gs[bj].rgb;
    return [base[0] * (1 - al) + c[0] * al, base[1] * (1 - al) + c[1] * al, base[2] * (1 - al) + c[2] * al];
  }
  let r = base[0], g = base[1], b = base[2], any = false;
  for (let j = 0; j < gs.length; j++) {
    const x = gs[j].ish[p], [lo, hi] = rngs[j]; if (x <= lo) continue;
    any = true; const f = Math.min(1, (x - lo) / (hi - lo));
    r += gs[j].rgb[0] * f; g += gs[j].rgb[1] * f; b += gs[j].rgb[2] * f;
  }
  return any ? [Math.min(255, r), Math.min(255, g), Math.min(255, b)] : null;
}

// Per-voxel summary of all shown genes, one tight pass per gene; adding a gene to the end only merges that gene.
let HITS = null, BEST = null, WHO = null, COMP = null;
function mergeIsh(gs, rngs, thr, dominant) {
  const n = brainMap().n, params = `${dominant}|${S.ishlo}|${S.ishhi}|${S.ishthr}`, syms = gs.map((g) => g.symbol);
  if (!HITS) HITS = new Uint8Array(n);
  if (dominant && !BEST) { BEST = new Uint8Array(n); WHO = new Uint8Array(n); }
  const extend = COMP && COMP.params === params && !gs.some((g) => g.pseudo) && COMP.syms.length <= syms.length && COMP.syms.every((x, k) => x === syms[k]);
  const from = extend ? COMP.syms.length : 0;
  if (!from) { HITS.fill(0); if (dominant) { BEST.fill(0); WHO.fill(0); } }
  const f = new Uint8Array(256);
  for (let j = from; j < gs.length; j++) {
    const v = gs[j].ish, t = thr[j];
    if (!dominant) { for (let p = 0; p < n; p++) if (v[p] >= t) HITS[p]++; continue; }
    const [lo, hi] = rngs[j];
    for (let x = 0; x < 256; x++) f[x] = x <= lo ? 0 : 1 + Math.round(Math.min(1, (x - lo) / (hi - lo)) * 254);
    for (let p = 0; p < n; p++) { const x = v[p]; if (x >= t) { HITS[p]++; if (f[x] > BEST[p]) { BEST[p] = f[x]; WHO[p] = j; } } }
  }
  COMP = { params, syms };
}

function buildIshPoints() {
  if (ishPts) { scene.remove(ishPts); ishPts.geometry.dispose(); ishPts = null; }
  const gs = shown().filter((g) => g.ish);
  if (!gs.length || !ANN) { hud(); return; }
  const rngs = gs.map((g) => ishRange(g)), thr = gs.map((g, j) => (g.thr3d != null ? g.thr3d : S.ishthr / 100 * rngs[j][2]));
  const want = S.cellshow === 'every' && gs.length > 1 ? gs.length : 1, dominant = gs.length > ADDITIVE_MAX;
  const rm = S.onlyRegion && S.region ? regionMask(S.region) : hpfRestrict();
  const { starts, lens } = brainMap();
  mergeIsh(gs, rngs, thr, dominant);
  let n = 0;
  for (let r = 0, p = 0; r < starts.length; r++) for (let i = starts[r], e = i + lens[r]; i < e; i++, p++) if (HITS[p] >= want && (!rm || rm[ANN[i]])) n++;
  const cap = 3_000_000, stride = n > cap ? Math.ceil(n / cap) : 1;
  const pos = new Float32Array(Math.ceil(n / stride) * 3), col = new Uint8Array(pos.length);
  let m = 0, seen = 0;
  const base = [0, 0, 0], grey = [40, 40, 40], sliceN = ND * NL;
  for (let r = 0, p = 0; r < starts.length; r++) for (let i = starts[r], e = i + lens[r]; i < e; i++, p++) {
    if (HITS[p] < want || (rm && !rm[ANN[i]])) continue;
    if ((seen++ % stride) !== 0) continue;
    const a = Math.floor(i / sliceN), rest = i - a * sliceN, d = Math.floor(rest / NL), s = IO.toScene(a * 25, d * 25, (rest - d * NL) * 25);
    pos[m * 3] = s[0]; pos[m * 3 + 1] = s[1]; pos[m * 3 + 2] = s[2];
    if (dominant) {
      const c = BEST[p] ? gs[WHO[p]].rgb : grey, al = BEST[p] ? 0.25 + 0.75 * (BEST[p] - 1) / 254 : 1;
      col[m * 3] = c[0] * al; col[m * 3 + 1] = c[1] * al; col[m * 3 + 2] = c[2] * al;
    } else {
      const c = ishShade(p, gs, rngs, base) || grey;
      col[m * 3] = c[0]; col[m * 3 + 1] = c[1]; col[m * 3 + 2] = c[2];
    }
    m++;
  }
  const g = new THREE.BufferGeometry();
  g.setAttribute('position', new THREE.BufferAttribute(pos.subarray(0, m * 3), 3));
  g.setAttribute('color', new THREE.BufferAttribute(col.subarray(0, m * 3), 3, true));
  ishPts = new THREE.Points(g, new THREE.PointsMaterial({ size: 0.032, vertexColors: true, sizeAttenuation: true, transparent: true, opacity: 0.85, depthWrite: false }));
  ishPts.layers.set(1);
  ishPts.userData = { n: m, stride };
  scene.add(ishPts);
  hud(); requestDraw();
}

async function updateRegionMesh() {
  if (regionMesh) { scene.remove(regionMesh); regionMesh.geometry.dispose(); regionMesh = null; }
  if (!S.region || !S.showMesh || !ONTO[S.region].mesh) { requestDraw(); return; }
  const r = S.region, m = await IO.mesh(ONTO[r].id);
  if (r !== S.region) return;
  const [cr, cg, cb] = ONTO[r].rgb;
  regionMesh = meshObject(m, new THREE.Color(cr / 255, cg / 255, cb / 255), 0.28);
  scene.add(regionMesh);
  requestDraw();
}

function draw() {
  if (S.view === '3d') draw3D(); else drawSlice();
}

function draw3D() {
  const w = renderer.domElement.clientWidth, h = renderer.domElement.clientHeight;
  const split = S.mode === 'split';
  $('lblL').style.display = $('lblR').style.display = split ? 'block' : 'none';
  if (shellMesh) shellMesh.visible = S.shell;
  renderer.setScissorTest(split);
  const pass = (x, width, layers) => {
    renderer.setViewport(x, 0, width, h); renderer.setScissor(x, 0, width, h);
    camera.aspect = width / h; camera.updateProjectionMatrix();
    camera.layers.set(0); layers.forEach((l) => camera.layers.enable(l));
    renderer.render(scene, camera);
  };
  renderer.clear();
  if (split) { pass(0, w / 2, [1]); pass(w / 2, w / 2, [2]); }
  else pass(0, w, ishOnly() ? [1] : S.mode === 'merfish' ? [2] : [1, 2]);
  camera.layers.enableAll();
}

// ───────────────────────── slices ─────────────────────────
const SL = { zoom: 1, tx: 0, ty: 0, img: null, key: '', hits: null, layout: null };
const PLANES = {
  coronal: { W: NL, H: ND, n: NA, lbl: 'AP', fix: 0 },
  sagittal: { W: NA, H: ND, n: NL, lbl: 'LR', fix: 2 },
  horizontal: { W: NL, H: NA, n: ND, lbl: 'DV', fix: 1 },
};
function voxAt(plane, k, u, v) {   // plane coords -> [ap, dv, lr]
  if (plane === 'coronal') return [k, v, u];
  if (plane === 'sagittal') return [u, v, k];
  return [v, k, u];
}

function sliceImage(plane, k, withIsh) {
  const P = PLANES[plane], img = new ImageData(P.W, P.H), px = img.data;
  const gs = withIsh ? shown().filter((g) => g.ish) : [], rngs = gs.map((g) => ishRange(g)), base = [0, 0, 0], hm = hpfRestrict();
  const lab = new Uint16Array(P.W * P.H);
  if (gs.length) brainMap();
  for (let v = 0; v < P.H; v++) for (let u = 0; u < P.W; u++) {
    const [a, d, l] = voxAt(plane, k, u, v), i = (a * ND + d) * NL + l, o = (v * P.W + u) * 4;
    const an = ANN ? ANN[i] : 0; lab[v * P.W + u] = an;
    let r = 0, g = 0, b = 0;
    if (an) { const t = REFVOL ? REFVOL[i] : 70; r = g = b = 18 + t * 0.42; }
    if (gs.length && an && (!hm || hm[an])) {
      base[0] = r; base[1] = g; base[2] = b;
      const c = ishShade(packedIndex(i), gs, rngs, base); if (c) { r = c[0]; g = c[1]; b = c[2]; }
    }
    px[o] = r; px[o + 1] = g; px[o + 2] = b; px[o + 3] = an ? 255 : 0;
  }
  const sel = S.region ? regionMask(S.region) : null;
  for (let v = 0; v < P.H; v++) for (let u = 0; u < P.W; u++) {
    const c = lab[v * P.W + u], rr = u + 1 < P.W ? lab[v * P.W + u + 1] : c, dd = v + 1 < P.H ? lab[(v + 1) * P.W + u] : c;
    const o = (v * P.W + u) * 4;
    if (sel && sel[c] && (!sel[rr] || !sel[dd])) { px[o] = 127; px[o + 1] = 230; px[o + 2] = 211; px[o + 3] = 255; }
    else if (c !== rr || c !== dd) { px[o] = Math.min(255, px[o] + 70); px[o + 1] = Math.min(255, px[o + 1] + 70); px[o + 2] = Math.min(255, px[o + 2] + 70); px[o + 3] = 255; }
  }
  const cv = document.createElement('canvas'); cv.width = P.W; cv.height = P.H; cv.getContext('2d').putImageData(img, 0, 0);
  return { canvas: cv, lab };
}

async function ensureRef() {
  if (REFVOL || refPromise) return refPromise;
  status('loading the CCFv3BBP template for slices…');
  refPromise = IO.array(IO.data('ref/' + REF.reference.file), 'uint8', progress).then((a) => { REFVOL = a; SL.key = ''; progress(0, 0); status(''); requestDraw(); });
  return refPromise;
}

function drawSlice() {
  const cv = $('slice'), ctx = cv.getContext('2d'), dpr = devicePixelRatio, W = cv.width, H = cv.height;
  ctx.setTransform(1, 0, 0, 1, 0, 0); ctx.clearRect(0, 0, W, H);
  if (!ANN) { ctx.fillStyle = '#9fb0aa'; ctx.font = `${13 * dpr}px JetBrains Mono`; ctx.fillText('loading annotation…', 20 * dpr, 40 * dpr); return; }
  const plane = S.view, P = PLANES[plane], k = S.pos[plane];
  const panels = S.mode === 'split' ? [{ ish: true, cells: false }, { ish: false, cells: true }]
    : [{ ish: S.mode !== 'merfish', cells: !ishOnly() }];
  const pw = W / panels.length;
  const fit = Math.min(pw / P.W, H / P.H) * 0.94 * SL.zoom;
  if (SL.focus) {   // keep a focus point (voxel AP, DV, LR) in the middle of each panel
    const [fa, fd, fl] = SL.focus, fu = plane === 'sagittal' ? fa : fl, fv = plane === 'horizontal' ? fa : fd;
    SL.tx = fit * (P.W / 2 - fu) / dpr; SL.ty = fit * (P.H / 2 - fv) / dpr;
  }
  const key = `${plane}|${k}|${S.gene}|${S.cmap}|${S.ishlo}|${S.ishhi}|${S.region}|${!!REFVOL}|${shown().map((g) => g.symbol).join('+')}|${!!hpfRestrict()}`;
  if (SL.key !== key) {
    SL.key = key;
    SL.imgIsh = sliceImage(plane, k, true);
    SL.imgPlain = S.mode === 'split' || S.mode === 'merfish' ? sliceImage(plane, k, false) : SL.imgIsh;
  }
  SL.layout = { plane, k, fit, pw, panels };
  SL.hits = [];
  panels.forEach((p, j) => {
    const ox = j * pw + (pw - P.W * fit) / 2 + SL.tx * dpr, oy = (H - P.H * fit) / 2 + SL.ty * dpr;
    p.ox = ox; p.oy = oy;
    ctx.save(); ctx.beginPath(); ctx.rect(j * pw, 0, pw, H); ctx.clip();
    ctx.imageSmoothingEnabled = false;
    ctx.drawImage((p.ish ? SL.imgIsh : SL.imgPlain).canvas, ox, oy, P.W * fit, P.H * fit);
    if (p.cells) drawCellsOnSlice(ctx, plane, k, ox, oy, fit, j);
    ctx.restore();
  });
  if (panels.length > 1) { ctx.fillStyle = 'rgba(255,255,255,.12)'; ctx.fillRect(pw - dpr, 0, 2 * dpr, H); }
  $('slicev').textContent = `${P.lbl} ${(k * 25 / 1000).toFixed(3)} mm (#${k})`;
  gapMessage(plane, k);
}

function gapMessage(plane, k) {
  const el = $('gapmsg'), b = B[S.brain];
  if (ishOnly() || !b || !b.sections) { el.style.display = 'none'; return; }
  const sp = b.sections, drawn = SL.hits && SL.hits.some((h) => h && h.grid.size);
  if (drawn) { el.style.display = 'none'; return; }
  let msg = `no ${b.man.dataset} cells within ±${S.slab} µm of this ${plane} plane`;
  if (sp.orient === plane) {
    const cur = k * 25, list = sp.pos.map((x) => (plane === 'coronal' ? x : IO.CENTER[2] + x));
    const near = list.reduce((a, x) => (Math.abs(x - cur) < Math.abs(a - cur) ? x : a), list[0]);
    msg += ` — this brain was imaged in ${sp.pos.length} ${sp.orient} sections ~${Math.round(sp.spacing)} µm apart; nearest section ${(near / 1000).toFixed(2)} mm (◀ ▶ to step)`;
  } else msg += ` — this brain was cut ${sp.orient}ly (${sp.pos.length} sections); use the ${sp.orient} view or a thicker slab`;
  el.textContent = msg; el.style.display = 'block';
}

function drawCellsOnSlice(ctx, plane, k, ox, oy, fit, panel) {
  const b = B[S.brain]; if (!b || !b.um) return;
  const um = b.um, col = b.colors, vis = cellVisible(), slab = S.slab, c0 = k * 25;
  const ax = PLANES[plane].fix;
  let i0 = 0, i1 = b.n;
  if (plane === 'coronal') { i0 = lowerBound(um, c0 - slab); i1 = lowerBound(um, c0 + slab + 1); }
  const r = Math.max(0.9 * devicePixelRatio, S.psize / 25 * fit / 2);
  const grid = new Map(), cell = 10 * devicePixelRatio;
  for (let i = i0; i < i1; i++) {
    if (Math.abs(um[i * 3 + ax] - c0) > slab || !vis(i)) continue;
    let u, v;
    if (plane === 'coronal') { u = um[i * 3 + 2]; v = um[i * 3 + 1]; }
    else if (plane === 'sagittal') { u = um[i * 3]; v = um[i * 3 + 1]; }
    else { u = um[i * 3 + 2]; v = um[i * 3]; }
    const x = ox + u / 25 * fit, y = oy + v / 25 * fit;
    ctx.fillStyle = `rgb(${col[i * 3]},${col[i * 3 + 1]},${col[i * 3 + 2]})`;
    ctx.fillRect(x - r, y - r, 2 * r, 2 * r);
    const gk = `${panel}:${Math.floor(x / cell)}:${Math.floor(y / cell)}`;
    if (!grid.has(gk)) grid.set(gk, []);
    grid.get(gk).push(i, x, y);
  }
  SL.hits[panel] = { grid, cell };
}
function lowerBound(um, ap) { let lo = 0, hi = um.length / 3; while (lo < hi) { const m = (lo + hi) >> 1; if (um[m * 3] < ap) lo = m + 1; else hi = m; } return lo; }

function sliceHit(ev) {
  if (!SL.layout || !ANN) return null;
  const dpr = devicePixelRatio, rect = $('slice').getBoundingClientRect();
  const x = (ev.clientX - rect.left) * dpr, y = (ev.clientY - rect.top) * dpr;
  const { plane, k, fit, pw, panels } = SL.layout, j = Math.min(panels.length - 1, Math.floor(x / pw)), p = panels[j];
  const u = Math.floor((x - p.ox) / fit), v = Math.floor((y - p.oy) / fit), P = PLANES[plane];
  if (u < 0 || v < 0 || u >= P.W || v >= P.H) return { x, y, j };
  const [a, d, l] = voxAt(plane, k, u, v), i = (a * ND + d) * NL + l;
  let cellHit = null;
  const h = SL.hits && SL.hits[j];
  if (h) {
    let best = (6 * dpr) ** 2;
    for (let gx = -1; gx <= 1; gx++) for (let gy = -1; gy <= 1; gy++) {
      const arr = h.grid.get(`${j}:${Math.floor(x / h.cell) + gx}:${Math.floor(y / h.cell) + gy}`); if (!arr) continue;
      for (let q = 0; q < arr.length; q += 3) { const dd = (arr[q + 1] - x) ** 2 + (arr[q + 2] - y) ** 2; if (dd < best) { best = dd; cellHit = arr[q]; } }
    }
  }
  return { x, y, j, a, d, l, i, region: ANN[i], cell: cellHit };
}

const TIP_GENES = 6;       // hover lists at most this many genes, strongest first
function tipFor(hit) {
  const lines = [];
  if (hit.region != null) {
    const o = ONTO[hit.region];
    lines.push(`<b>${o.acronym}</b> ${o.name}`);
    lines.push(`AP ${(hit.a * 25 / 1000).toFixed(2)} · DV ${(hit.d * 25 / 1000).toFixed(2)} · LR ${(hit.l * 25 / 1000).toFixed(2)} mm`);
    if (hit.region) {
      const p = packedIndex(hit.i), vals = [];
      for (const g of SEL.filter((x) => x.ish)) {
        if (g.tip) { const t = g.tip(hit.i); if (t) lines.push(t); continue; }
        const x = g.ish[p], top = g.p995 != null ? g.p995 : ISH.genes[g.ishIdx].p995;
        vals.push({ g, x, rel: x / Math.max(1, top) });
      }
      const many = vals.length > TIP_GENES, list = many ? vals.filter((v) => v.x).sort((a, b) => b.rel - a.rel) : vals;
      for (const { g, x, rel } of list.slice(0, TIP_GENES)) lines.push(`${g.pseudo ? '' : 'ISH '}${g.label || g.symbol}: ${x} (${(rel * 100).toFixed(0)}% of p99.5)`);
      if (many) lines.push(`&nbsp;ISH detected in ${list.length} of ${vals.length} genes here${list.length > TIP_GENES ? ` · strongest ${TIP_GENES} shown` : ''}`);
    }
  }
  if (hit.cell != null) {
    const b = B[S.brain], c = b.cluster[hit.cell], q = b.quality[hit.cell] / 255;
    lines.push(`<b>cell</b> ${b.man.dataset} #${hit.cell}`);
    for (const lv of ['class', 'subclass', 'supertype', 'cluster']) lines.push(`&nbsp;${lv}: ${TAX[lv][CL2[lv][c]].name}`);
    lines.push(`&nbsp;own CCF: ${ONTO[b.region_ccf[hit.cell]].acronym} · CCFv3BBP: ${ONTO[b.region_bbp[hit.cell]].acronym}`);
    lines.push(`&nbsp;${b.man.quality_field}: ${q.toFixed(2)}`);
    const mers = SEL.filter((x) => x.mer), manyM = mers.length > TIP_GENES;
    const mlist = manyM ? mers.filter((x) => x.mer[hit.cell]).sort((x, y) => y.mer[hit.cell] - x.mer[hit.cell]) : mers;
    for (const sg of mlist.slice(0, TIP_GENES)) {
      const g = b.man.genes[sg.merIdx], q = sg.mer[hit.cell];
      lines.push(`&nbsp;${sg.symbol} log2: ${q ? (q / 255 * g.vmax).toFixed(2) : 'not detected'}`);
    }
    if (manyM) lines.push(`&nbsp;detected: ${mlist.length} of ${mers.length} genes${mlist.length > TIP_GENES ? ` · strongest ${TIP_GENES} shown` : ''}`);
  }
  return lines.join('<br>');
}

// ───────────────────────── region profile ─────────────────────────
function regionCounts(rc, r) {
  let lo = 0, hi = rc.region.length;
  while (lo < hi) { const m = (lo + hi) >> 1; if (rc.region[m] < r) lo = m + 1; else hi = m; }
  const start = lo; hi = rc.region.length;
  while (lo < hi) { const m = (lo + hi) >> 1; if (rc.region[m] <= r) lo = m + 1; else hi = m; }
  return [start, lo];
}

async function refreshProfile() {
  const el = $('profile');
  if (!S.region) { el.innerHTML = '<div class="status">pick a region (search, or click a slice) to see each brain\'s cell-type profile and the gene\'s values there</div>'; return; }
  const r = S.region, o = ONTO[r], lv = S.level, toLv = CL2[lv];
  let html = `<h5><span class="sw" style="background:rgb(${o.rgb})"></span>${o.acronym} · ${o.name}</h5>
    <div class="mini">${(o.vox * 25 ** 3 / 1e9).toFixed(3)} mm³ in CCFv3BBP · cells by ${S.labelsrc === 'ccf' ? 'own Allen CCF label' : 'CCFv3BBP label'} · ${lv}</div>`;
  const rcs = await Promise.all(DS.merfish.map((d) => loadRegionCluster(d.key, S.labelsrc)));
  if (S.region !== r) return;
  DS.merfish.forEach((d, j) => {
    const rc = rcs[j], [s, e] = regionCounts(rc, r), agg = new Map(); let tot = 0;
    for (let q = s; q < e; q++) { const t = toLv[rc.cluster[q]]; agg.set(t, (agg.get(t) || 0) + rc.count[q]); tot += rc.count[q]; }
    const top = [...agg.entries()].sort((x, y) => y[1] - x[1]);
    const cur = d.key === S.brain ? ' style="font-weight:600;color:var(--ink)"' : '';
    html += `<div class="brainrow"><div class="head"><span${cur}>${d.short} · ${d.dataset}</span><span>${tot.toLocaleString()} cells · ${agg.size} ${lv === 'class' ? 'classes' : lv === 'subclass' ? 'subclasses' : lv + 's'}</span></div>`;
    if (tot) {
      html += '<div class="stackbar">' + top.slice(0, 14).map(([t, n]) => `<div title="${TAX[lv][t].name}: ${n} (${(n / tot * 100).toFixed(1)}%)" style="width:${n / tot * 100}%;background:${TAX[lv][t].color}"></div>`).join('') + '<div style="flex:1;background:#d9ddda"></div></div>';
      html += '<div class="types">' + top.slice(0, d.key === S.brain ? 10 : 4).map(([t, n]) => `<div><span class="sw" style="background:${TAX[lv][t].color}"></span>${TAX[lv][t].name}<span class="pct">${(n / tot * 100).toFixed(1)}% · ${n}</span></div>`).join('') + '</div>';
    }
    html += '</div>';
  });
  html += compareBrainsHtml(rcs, r, lv, toLv);
  for (const g of SEL) html += await geneInRegionHtml(r, g);
  html += await topGenesHtml(r);
  if (S.region === r) el.innerHTML = html;
  el.querySelectorAll('[data-gene]').forEach((a) => a.addEventListener('click', (ev) => { ev.preventDefault(); selectGene(a.dataset.gene); }));
}

function compareBrainsHtml(rcs, r, lv, toLv) {
  // each brain's composition is computed on its own; the table only places them side by side
  const per = rcs.map((rc) => {
    const [s, e] = regionCounts(rc, r), agg = new Map(); let tot = 0;
    for (let q = s; q < e; q++) { const t = toLv[rc.cluster[q]]; agg.set(t, (agg.get(t) || 0) + rc.count[q]); tot += rc.count[q]; }
    return { agg, tot };
  });
  const score = new Map();
  per.forEach(({ agg, tot }) => { if (tot) agg.forEach((n, t) => score.set(t, Math.max(score.get(t) || 0, n / tot))); });
  const top = [...score.entries()].sort((x, y) => y[1] - x[1]).slice(0, 12).map(([t]) => t);
  if (!top.length) return '';
  const head = DS.merfish.map((d) => `<th class="n" title="${d.label}">${d.short.replace('Allen ', 'A').replace('Zhuang ', 'Z')}</th>`).join('');
  const rows = top.map((t) => `<tr><td><span class="sw" style="background:${TAX[lv][t].color}"></span>${TAX[lv][t].name}</td>` +
    per.map(({ agg, tot }) => { const f = tot ? (agg.get(t) || 0) / tot : NaN; return `<td class="n" style="background:rgba(20,125,112,${isNaN(f) ? 0 : Math.min(0.85, f * 1.6)});color:${f > 0.3 ? '#fff' : 'inherit'}">${isNaN(f) ? '—' : (f * 100).toFixed(1)}</td>`; }).join('') + '</tr>').join('');
  return `<h5>${lv} composition, brain by brain (% of each brain's cells in ${ONTO[r].acronym})</h5><table class="t"><tr><th>${lv}</th>${head}</tr>${rows}</table>`;
}

let CONC = null;
async function renderGeneRegions() {
  const el = $('generegions');
  if (!SEL.length) { el.innerHTML = ''; return; }
  const want = S.gene;
  let html = '';
  const tables = SEL.length > TABLE_GENES ? SEL.filter((g) => g === GENE) : SEL;   // many genes: the clicked chip's table only
  if (tables.length < SEL.length) html += `<div class="mini" style="margin:6px 0">${SEL.length} genes selected · click a gene above to see its top regions</div>`;
  for (const g of tables) html += await geneRegionsHtml(g);
  if (want !== S.gene) return;
  el.innerHTML = html;
  el.querySelectorAll('[data-region]').forEach((a) => a.addEventListener('click', (ev) => { ev.preventDefault(); selectRegion(+a.dataset.region); }));
}

const TABLE_GENES = 6;
async function geneRegionsHtml(G0) {
  if (G0.pseudo && G0.kind === 'layer') return `<h5><span class="sw" style="background:${G0.color}"></span>${G0.label}</h5><div class="mini">combined from the HGEA-ISH “All genes” tab (Blue Brain ISH, 50 µm, hippocampus)</div>`;
  if (G0.pseudo) return `<h5><span class="sw" style="background:${G0.color}"></span>${G0.label}</h5><div class="mini">${G0.kind === 'signature' ? 'mean' : 'contrast of the'} Blue Brain ISH of the genes Bienkowski et al. annotate in ${G0.col} (Supplementary Table 2)</div>`;
  const e = GENES.get(G0.symbol), sym = G0.symbol;
  if (!CONC) CONC = await IO.json(IO.data('qc/concordance.json')).catch(() => ({ ish_vs_merfish: {} }));
  const okReg = (i) => ONTO[i] && DEPTH[i] >= 5 && DEPTH[i] <= 8 && ONTO[i].vox >= 500;
  let html = '';
  const cols = [];
  if (e.ish >= 0) {
    const R = await loadIshRegions(), G = R.G, v = [];
    for (let i = 1; i < ONTO.length; i++) if (okReg(i)) v.push([i, R.mean[i * G + e.ish]]);
    v.sort((x, y) => y[1] - x[1]); cols.push({ name: 'ISH', top: v.slice(0, 10), fmt: (x) => x.toFixed(1) });
  }
  const gi = e.mer[S.brain];
  if (gi != null) {
    const b = await loadRegionGene(S.brain);
    if (b.rmean) {
      const G = b.man.genes.length, v = [];
      for (let i = 1; i < ONTO.length; i++) if (okReg(i) && b.man.region_cells[i] >= 50) v.push([i, b.rmean[i * G + gi]]);
      v.sort((x, y) => y[1] - x[1]); cols.push({ name: `MERFISH ${DS.merfish.find((d) => d.key === S.brain).short}`, top: v.slice(0, 10), fmt: (x) => x.toFixed(2) });
    }
  }
  if (cols.length) {
    html += `<h5>${multi() ? `<span class="sw" style="background:${G0.color}"></span>` : ''}Where is ${sym}? top regions (depth 5–8)</h5><table class="t"><tr>${cols.map((c) => `<th>${c.name}</th>`).join('')}</tr>`;
    for (let k = 0; k < 10; k++) html += '<tr>' + cols.map((c) => { const x = c.top[k]; return x ? `<td><a href="#" data-region="${x[0]}">${ONTO[x[0]].acronym}</a> <span class="mini">${c.fmt(x[1])}</span></td>` : '<td></td>'; }).join('') + '</tr>';
    html += '</table>';
  }
  const rho = DS.merfish.map((d) => { const r = (CONC.ish_vs_merfish[d.key] || {}).genes; if (!r || !r[sym]) return null; const x = r[sym];
    return `<tr><td>${d.short}</td><td class="n">${x[0].toFixed(2)}</td><td class="n">${x[2].toFixed(2)}</td><td class="n">${Math.min(x[3], x[4]).toFixed(2)}</td></tr>`; }).filter(Boolean);
  if (rho.length) html += `<h5>ISH vs MERFISH agreement for ${sym}</h5><table class="t"><tr><th>brain</th><th>Spearman</th><th>Pearson</th><th>peak</th></tr>${rho.join('')}</table>
    <div class="mini" style="margin-top:3px">across CCF regions (depth 5–8); peak = percentile of each atlas' top region in the other (1.00 = same top region). <a href="../qc/index.html">details</a></div>`;
  return html;
}

async function geneInRegionHtml(r, G0) {
  if (G0.pseudo) return '';
  const e = GENES.get(G0.symbol); let rows = '';
  if (e.ish >= 0) {
    const R = await loadIshRegions(), G = R.G, v = R.mean[r * G + e.ish], root = R.mean[ONTO.findIndex((o) => o && o.id === 997) * G + e.ish];
    rows += `<tr><td>ISH (Blue Brain)</td><td class="n">${v.toFixed(2)}</td><td class="n">${(R.frac[r * G + e.ish] * 100).toFixed(0)}% vox</td><td class="n">${(v / Math.max(root, 1e-6)).toFixed(2)}×</td></tr>`;
  }
  for (const d of DS.merfish) {
    const gi = e.mer[d.key]; if (gi == null) continue;
    const b = await loadRegionGene(d.key); if (!b.rmean) continue;
    const G = b.man.genes.length, n = b.man.region_cells[r], v = b.rmean[r * G + gi], root = b.rmean[ONTO.findIndex((o) => o && o.id === 997) * G + gi];
    rows += `<tr><td>${d.short}</td><td class="n">${n ? v.toFixed(2) : '—'}</td><td class="n">${n ? (b.rfrac[r * G + gi] * 100).toFixed(0) + '% cells' : '0 cells'}</td><td class="n">${n ? (v / Math.max(root, 1e-6)).toFixed(2) + '×' : ''}</td></tr>`;
  }
  return `<h5>${multi() ? `<span class="sw" style="background:${G0.color}"></span>` : ''}${G0.symbol} in ${ONTO[r].acronym}</h5><table class="t"><tr><th>dataset</th><th>mean</th><th>detected</th><th>vs brain</th></tr>${rows}</table>
    <div class="mini" style="margin-top:3px">ISH: mean voxel value (0–255); MERFISH: mean log2 expression per cell (own CCF labels). Not comparable in absolute units — compare the "vs brain" enrichment.</div>`;
}

async function topGenesHtml(r) {
  const R = await loadIshRegions(), G = R.G, root = ONTO.findIndex((o) => o && o.id === 997);
  const enr = [];
  for (let g = 0; g < G; g++) { const m0 = R.mean[root * G + g], m = R.mean[r * G + g]; if (m0 > 0.5 && m > 2) enr.push([g, m / m0]); }
  enr.sort((x, y) => y[1] - x[1]);
  let html = `<h5>Most enriched genes in ${ONTO[r].acronym}</h5><div class="mini">ISH: ` +
    enr.slice(0, 14).map(([g, x]) => `<a href="#" data-gene="${ISH.genes[g].gene}">${ISH.genes[g].gene}</a> ${x.toFixed(1)}×`).join(' · ') + '</div>';
  const b = await loadRegionGene(S.brain);
  if (b.rmean && b.man.region_cells[r] >= 30) {
    const Gb = b.man.genes.length, me = [];
    for (let g = 0; g < Gb; g++) { if (isControl(b.man.genes[g].symbol)) continue; const m0 = b.rmean[root * Gb + g], m = b.rmean[r * Gb + g]; if (m0 > 0.02 && b.rfrac[r * Gb + g] > 0.1) me.push([g, m / m0]); }
    me.sort((x, y) => y[1] - x[1]);
    html += `<div class="mini" style="margin-top:4px">MERFISH ${DS.merfish.find((d) => d.key === S.brain).short}: ` +
      me.slice(0, 14).map(([g, x]) => `<a href="#" data-gene="${b.man.genes[g].symbol}">${b.man.genes[g].symbol}</a> ${x.toFixed(1)}×`).join(' · ') + '</div>';
  }
  return html;
}

function selectRegion(r, fly = true) {
  S.region = r; $('region').value = r ? `${ONTO[r].acronym} · ${ONTO[r].name}` : ''; writeURL();
  updateRegionMesh(); if (S.onlyRegion) { recolorCells(); buildIshPoints(); }
  SL.key = ''; refreshProfile(); requestDraw();
}

// ───────────────────────── legend, HUD, status ─────────────────────────
function renderLegend() {
  const el = $('legend'); let html = '';
  if (multi()) {
    const dom = shown().length > ADDITIVE_MAX;
    html += `<h4>${SEL.length} genes (${dom ? 'strongest gene per spot' : 'additive colour'}) · click to toggle</h4><div class="glegend">` + SEL.map((g, j) => `<div class="item gtoggle${g.hidden ? ' off' : ''}" data-g="${j}" title="click to ${g.hidden ? 'show' : 'hide'} ${g.symbol}"><span class="sw" style="background:${g.hidden ? 'transparent' : g.color};border-color:${g.color}"></span>${g.label || g.symbol}<span style="margin-left:auto;font-family:var(--mono);font-size:10px;color:#9fb0aa">${g.ish ? (lowIsh(g) ? 'ISH low' : 'ISH') : ''}${g.ish && g.mer ? '·' : ''}${g.mer ? 'MER' : ''}</span></div>`).join('') + '</div>' +
      `<div class="mini" style="color:#9fb0aa;margin-top:4px">${dom ? `with more than ${ADDITIVE_MAX} genes shown, each spot takes the colour of its strongest gene` : `overlap mixes colours${SEL.length === 2 ? ' (both ≈ white)' : ''}`}; brightness = level within each gene's range</div>`;
  } else if (GENE && (GENE.ish || GENE.mer)) {
    html += `<h4>${GENE.label || GENE.symbol}</h4>`;
    if (GENE.legend) { const L = GENE.legend; html += `<div class="mini" style="color:#cdd8d3">${L.unit}</div><div class="cb" style="background:${gradientCss(S.cmap)}"></div><div class="ticks"><span>${L.lo}</span><span>${L.hi}</span></div>`; }
    else if (GENE.ish && S.mode !== 'merfish') { const [lo, hi] = ishRange(); html += `<div class="mini" style="color:#cdd8d3">${GENE.pseudo ? `${GENE.kind} (×255)` : 'ISH value'}</div><div class="cb" style="background:${gradientCss(S.cmap)}"></div><div class="ticks"><span>${lo.toFixed(0)}</span><span>${hi.toFixed(0)}</span></div>`; if (lowIsh(GENE)) html += `<div class="mini" style="color:#f0b27a;margin-top:4px">${lowIshNote(GENE)}</div>`; }
    if (GENE.mer && !ishOnly() && S.colorby === 'gene') {
      const g = B[S.brain].man.genes[GENE.merIdx];
      html += `<div class="mini" style="color:#cdd8d3;margin-top:6px">MERFISH log2 (${S.brain})</div><div class="cb" style="background:${gradientCss(S.cmap)}"></div><div class="ticks"><span>${(S.merlo / 100 * g.vmax).toFixed(1)}</span><span>${(S.merhi / 100 * g.vmax).toFixed(1)}</span></div>`;
      if (S.mode === 'merfish') html += '<div class="mini" style="color:#9fb0aa;margin-top:4px">weakest colours are raised so every cell with signal stays visible; signal is not proof of expression, and each brain has its own scale</div>';
    }
  }
  if (S.colorby === 'type' && !ishOnly() && B[S.brain] && B[S.brain].cluster) {
    const b = B[S.brain], lv = CL2[S.level], cnt = new Map();
    const idx = cellPts ? cellPts.geometry.index.array : [];
    for (let q = 0; q < idx.length; q += Math.max(1, Math.floor(idx.length / 200000))) { const t = lv[b.cluster[idx[q]]]; cnt.set(t, (cnt.get(t) || 0) + 1); }
    const top = [...cnt.entries()].sort((x, y) => y[1] - x[1]).slice(0, 18);
    html += `<h4 style="margin-top:8px">${S.level} (top in view)</h4>` + top.map(([t]) => `<div class="item"><span class="sw" style="background:${TAX[S.level][t].color}"></span>${TAX[S.level][t].name}</div>`).join('');
  }
  el.innerHTML = html || '<h4>legend</h4><div class="mini" style="color:#9fb0aa">select a gene or colour cells by type</div>';
  el.querySelectorAll('.gtoggle').forEach((it) => it.addEventListener('click', () => toggleGene(SEL[+it.dataset.g])));
}
const ONE_SIDE = 0.9;      // share of a brain's cells on one side above which it counts as one hemisphere
function hud() {
  const b = B[S.brain]; const bits = [];
  if (b && b.n && !ishOnly()) bits.push(`<b>${(b.nVisible || 0).toLocaleString()}</b> / ${b.n.toLocaleString()} cells · ${b.man.dataset}` + (b.sections ? ` · ${b.sections.pos.length} ${b.sections.orient} sections ~${Math.round(b.sections.spacing)} µm apart` : ''));
  if (S.mode !== 'merfish' && !SEL.some((g) => g.ish)) bits.push(SEL.length ? `no ISH volume for ${SEL.map((g) => g.symbol).join(', ')}` : 'pick a gene to see ISH');
  if (ishPts && S.mode !== 'merfish') bits.push(`<b>${ishPts.userData.n.toLocaleString()}</b> ISH voxels${ishPts.userData.stride > 1 ? ` (1/${ishPts.userData.stride})` : ''}`);
  if (b && b.n && !ishOnly() && b.oneSide > ONE_SIDE) bits.push(`one hemisphere sampled (${Math.round(b.oneSide * 100)}% of cells): the other side is unsampled, not empty`);
  if (S.mode === 'merfish' && b && b.nDetected != null) bits.push(`<b>${b.nDetected.toLocaleString()}</b> with signal`);
  if (S.mode === 'merfish' && b) {      // say which selected genes this brain cannot show
    const real = SEL.filter((g) => !g.pseudo), absent = real.filter((g) => GENES.get(g.symbol).mer[S.brain] == null);
    const none = absent.filter((g) => !Object.keys(GENES.get(g.symbol).mer).length), other = absent.filter((g) => !none.includes(g));
    if (none.length) bits.push(`${none.map((g) => g.symbol).join(', ')}: on no MERFISH panel (ISH only)`);
    if (other.length) bits.push(`${other.map((g) => g.symbol).join(', ')}: not on this brain's panel`);
    if (absent.length === real.length && real.length) bits.push('cells coloured by type');
  }
  $('hud').innerHTML = bits.join(' · ') || 'no data shown';
}
function status(t) { $('genestatus').textContent = t; const h = $('hgeastatus'); if (h) h.textContent = t; }
function progress(got, total) {
  const p = $('prog'); if (!total) { p.style.display = 'none'; return; }
  p.style.display = 'block'; p.firstElementChild.style.width = `${Math.min(100, got / total * 100)}%`;
}

// ───────────────────────── UI ─────────────────────────
function setSeg(id, v) {
  document.querySelectorAll(`#${id} button`).forEach((b) => b.classList.toggle('active', b.dataset.v === v));
  const map = { mode: 'mode', labelsrc: 'labelsrc', colorby: 'colorby', level: 'level', view: 'view' };
  S[map[id]] = v;
}
function seg(id, fn) {
  document.querySelectorAll(`#${id} button`).forEach((b) => b.addEventListener('click', () => { setSeg(id, b.dataset.v); writeURL(); fn(b.dataset.v); }));
  setSeg(id, S[id]);
}
function setView(v) {
  setSeg('view', v);
  const is3 = v === '3d';
  $('gl').style.display = is3 ? 'block' : 'none'; $('slice').style.display = is3 ? 'none' : 'block';
  $('slicebar').style.display = is3 ? 'none' : 'flex'; $('hint').textContent = is3 ? 'drag rotate · scroll zoom · right-drag pan' : 'scroll zoom · drag pan · click a region to select · hover for details';
  if (is3) $('gapmsg').style.display = 'none';
  if (!is3) {
    const P = PLANES[v]; $('slicepos').max = P.n - 1; $('slicepos').value = S.pos[v]; $('slicelbl').textContent = P.lbl;
    sectionTicks();
    ensureRef(); SL.key = '';
  }
  writeURL(); requestDraw();
}

function suggest(inputId, listId, items, onPick) {
  const inp = $(inputId), list = $(listId); let hi = -1, cur = [];
  const show = () => {
    const q = inp.value.trim().toLowerCase(); list.innerHTML = '';
    if (!q) { list.style.display = 'none'; return; }
    cur = items(q).slice(0, 40); hi = -1;
    cur.forEach((it, k) => { const d = document.createElement('div'); d.innerHTML = it.html; d.addEventListener('mousedown', (e) => { e.preventDefault(); pick(k); }); list.appendChild(d); });
    list.style.display = cur.length ? 'block' : 'none';
  };
  const pick = (k) => { if (cur[k]) { list.style.display = 'none'; onPick(cur[k].value); inp.blur(); } };
  inp.addEventListener('input', show);
  inp.addEventListener('focus', show);
  inp.addEventListener('blur', () => setTimeout(() => { list.style.display = 'none'; }, 120));
  inp.addEventListener('keydown', (e) => {
    if (e.key === 'ArrowDown' || e.key === 'ArrowUp') { hi = (hi + (e.key === 'ArrowDown' ? 1 : -1) + cur.length) % cur.length; [...list.children].forEach((c, k) => c.classList.toggle('hi', k === hi)); e.preventDefault(); }
    if (e.key === 'Enter') { pick(hi >= 0 ? hi : 0); e.preventDefault(); }
  });
}

function initUI() {
  const sel = $('brain');
  DS.merfish.forEach((d) => { const o = document.createElement('option'); o.value = d.key; o.textContent = `${d.short} · ${d.label} · ${d.cells.toLocaleString()} cells`; sel.appendChild(o); });
  sel.value = S.brain; sel.addEventListener('change', async () => { await loadBrain(sel.value); writeURL(); renderAvail(); buildIshPoints(); sectionTicks(); SL.key = ''; requestDraw(); });
  CMAPS.forEach((c) => { const o = document.createElement('option'); o.value = c; o.textContent = c; $('cmap').appendChild(o); });
  $('cmap').value = S.cmap; $('cmap').addEventListener('change', () => { S.cmap = $('cmap').value; writeURL(); recolorCells(); buildIshPoints(); SL.key = ''; renderLegend(); });
  seg('mode', async (v) => { if (v === 'merfish') await showOnPanelBrain(); recolorCells(); if (v === 'hgea') await enterHgea(); else { $('hgeasec').style.display = 'none'; showHpfMeshes(false); SL.focus = null; SL.zoom = 1; SL.tx = SL.ty = 0; } buildIshPoints(); hud(); renderLegend(); SL.key = ''; requestDraw(); });
  seg('labelsrc', () => { recolorCells(); refreshProfile(); });
  seg('colorby', () => { recolorCells(); });
  seg('level', () => { recolorCells(); refreshProfile(); });
  document.querySelectorAll('#view button').forEach((b) => b.addEventListener('click', () => setView(b.dataset.v)));
  const slider = (id, key, fn, fmt) => { const el = $(id); el.value = S[key]; el.addEventListener('input', () => { S[key] = +el.value; if (fmt) fmt(); fn(); }); };
  let ishTimer = null;
  const ishLater = () => { clearTimeout(ishTimer); ishTimer = setTimeout(() => { buildIshPoints(); renderLegend(); }, 180); };
  slider('quality', 'quality', recolorCells, () => { $('qualityv').textContent = S.quality.toFixed(2); });
  slider('psize', 'psize', () => { recolorCells(); }, () => { $('psizev').textContent = `${S.psize} µm`; });
  slider('ishlo', 'ishlo', () => { SL.key = ''; ishLater(); requestDraw(); });
  slider('ishhi', 'ishhi', () => { SL.key = ''; ishLater(); requestDraw(); });
  slider('ishthr', 'ishthr', ishLater, () => { $('ishthrv').textContent = `${S.ishthr}%`; });
  slider('merlo', 'merlo', recolorCells); slider('merhi', 'merhi', recolorCells);
  $('cellshow').addEventListener('change', () => { S.cellshow = $('cellshow').value; recolorCells(); buildIshPoints(); SL.key = ''; });
  $('shell').addEventListener('change', () => { S.shell = $('shell').checked; requestDraw(); });
  $('onlyregion').addEventListener('change', () => { S.onlyRegion = $('onlyregion').checked; recolorCells(); buildIshPoints(); });
  $('showmesh').addEventListener('change', () => { S.showMesh = $('showmesh').checked; updateRegionMesh(); });
  $('clearregion').addEventListener('click', () => selectRegion(0));
  let tq = null;
  $('typeq').addEventListener('input', () => { clearTimeout(tq); tq = setTimeout(() => { S.typeq = $('typeq').value; recolorCells(); $('typestatus').textContent = S.typeq ? `${(B[S.brain].nVisible || 0).toLocaleString()} cells match "${S.typeq}" at ${S.level} level` : ''; }, 220); });
  $('slicepos').addEventListener('input', () => { S.pos[S.view] = +$('slicepos').value; requestDraw(); });
  $('slab').addEventListener('change', () => { S.slab = +$('slab').value; requestDraw(); });
  $('secprev').addEventListener('click', () => jumpSection(-1)); $('secnext').addEventListener('click', () => jumpSection(1));

  const symbols = [...GENES.keys()].sort((a, b) => a.localeCompare(b));
  const geneList = (text) => text.split(/[\s,;]+/).filter(Boolean);
  $('gene').addEventListener('paste', (ev) => {
    const list = geneList(ev.clipboardData.getData('text'));
    if (list.length > 1) { ev.preventDefault(); addGenes(list); }
  });
  $('gene').addEventListener('keydown', (ev) => {
    const list = geneList($('gene').value);
    if (ev.key === 'Enter' && list.length > 1) { ev.preventDefault(); ev.stopImmediatePropagation(); $('gene').value = ''; addGenes(list); }
  });
  suggest('gene', 'genelist', (q) => {
    const starts = symbols.filter((g) => g.toLowerCase().startsWith(q)), inc = symbols.filter((g) => !g.toLowerCase().startsWith(q) && g.toLowerCase().includes(q));
    return [...starts, ...inc].map((g) => { const e = GENES.get(g), nm = Object.keys(e.mer).length;
      return { value: g, html: `<span class="g">${g}</span><span class="n">${e.ish >= 0 ? 'ISH' : ''}${e.ish >= 0 && nm ? ' · ' : ''}${nm ? `MERFISH ×${nm}` : ''}</span>` }; });
  }, (g) => selectGene(g));
  renderGeneChips();
  suggest('region', 'regionlist', (q) => ONTO.filter((o) => o && o.i && (o.acronym.toLowerCase().includes(q) || o.name.toLowerCase().includes(q)))
    .sort((a, b) => (b.acronym.toLowerCase() === q) - (a.acronym.toLowerCase() === q) || a.depth - b.depth)
    .map((o) => ({ value: o.i, html: `<span class="sw" style="background:rgb(${o.rgb})"></span><span class="g">${o.acronym}</span> ${o.name}` })), (r) => selectRegion(r));

  const sc = $('slice'); let drag = null;
  sc.addEventListener('wheel', (e) => { e.preventDefault(); SL.zoom = Math.max(0.5, Math.min(30, SL.zoom * (e.deltaY < 0 ? 1.15 : 1 / 1.15))); requestDraw(); }, { passive: false });
  sc.addEventListener('mousedown', (e) => { SL.focus = null; drag = { x: e.clientX, y: e.clientY, tx: SL.tx, ty: SL.ty, moved: false }; });
  window.addEventListener('mouseup', (e) => {
    if (drag && !drag.moved && e.target === sc) { const h = sliceHit(e); if (h && h.region) selectRegion(h.region); }
    drag = null;
  });
  sc.addEventListener('mousemove', (e) => {
    if (drag) { const dx = e.clientX - drag.x, dy = e.clientY - drag.y; if (Math.abs(dx) + Math.abs(dy) > 3) drag.moved = true; SL.tx = drag.tx + dx; SL.ty = drag.ty + dy; requestDraw(); return; }
    const h = sliceHit(e), tip = $('tip');
    if (!h || (h.region == null && h.cell == null)) { tip.style.display = 'none'; return; }
    tip.innerHTML = tipFor(h); tip.style.display = 'block';
    const rect = $('stage').getBoundingClientRect();
    tip.style.left = `${Math.min(e.clientX - rect.left + 14, rect.width - 350)}px`; tip.style.top = `${Math.min(e.clientY - rect.top + 14, rect.height - 200)}px`;
  });
  sc.addEventListener('mouseleave', () => { $('tip').style.display = 'none'; });
  sc.addEventListener('dblclick', () => { SL.focus = null; SL.zoom = 1; SL.tx = SL.ty = 0; requestDraw(); });
  renderAvail(); renderLegend();
}

boot().catch((e) => { console.error(e); $('loading').textContent = `failed to load: ${e.message}`; });

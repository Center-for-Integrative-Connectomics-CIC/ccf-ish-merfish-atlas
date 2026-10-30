// HGEA-ISH: the Hippocampus Gene Expression Atlas gene set (Bienkowski et al. 2018, Supplementary Tables 2 and 7)
// on the Blue Brain ISH volumes. Gene names, domain names, presence marks and Table 7 are shown as published;
// signatures and the all-genes layer are derived from the ISH (documented in data/hgea/hgea.json).
import * as IO from './io.js?v=20261029212214';

const $ = (id) => document.getElementById(id);
let H = null, api = null, current = null;
const esc = (s) => String(s).replace(/[&<>"]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
const short = (c) => (c === 'Putative HPF Interneurons' ? 'Interneurons' : c);

export async function initHgea(callbacks) {
  api = callbacks;
  H = await IO.json(IO.data('hgea/hgea.json'));
  document.querySelectorAll('#htabs button').forEach((b) => b.addEventListener('click', () => openTab(b.dataset.t)));
  renderDomains();
  renderSource();
  return H;
}

export function openTab(t) {
  document.querySelectorAll('#htabs button').forEach((b) => b.classList.toggle('active', b.dataset.t === t));
  document.querySelectorAll('.hpane').forEach((p) => { p.style.display = p.dataset.t === t ? 'block' : 'none'; });
  if (t === 'all') ensureAll().then(() => { if (!ALL.opened) { ALL.opened = true; changed(); } });   // opening the tab shows the layer
}

const colIndex = (col) => H.columns.indexOf(col);
const genesOf = (col) => { const k = colIndex(col); return H.genes.filter((g) => g.domains.includes(k)); };
function t7For(col) {
  if (col.startsWith('SUB_')) { const L = +col.slice(4); return H.table7.filter((d) => d.layers && d.layers.includes(L)); }
  return H.table7.filter((d) => (d.table2_column || d.name) === col);
}
const domainGrid = (cls, extra = () => '') => H.groups.map((grp) => `<div class="hgrp"><div class="mini">${grp.group}</div><div class="hrow">` +
  grp.columns.map((c) => `<button class="${cls}" data-col="${esc(c)}" title="${genesOf(c).length} genes annotated in ${esc(c)} (Supplementary Table 2)">${esc(short(c))}<span>${genesOf(c).length}</span>${extra(c)}</button>`).join('') + '</div></div>').join('');

// ───────────── Domains tab ─────────────
function renderDomains() {
  $('hgeadomains').innerHTML = '<div class="mini" style="margin-bottom:2px"><b>1 ·</b> pick a domain (number = genes annotated there)</div>' + domainGrid('hdom');
  $('hgeadomains').querySelectorAll('.hdom').forEach((b) => b.addEventListener('click', () => selectDomain(b.dataset.col)));
}

export function selectDomain(col) {
  current = col;
  document.querySelectorAll('#hgeadomains .hdom').forEach((b) => b.classList.toggle('on', b.dataset.col === col));
  api.setDomain(col);
  const genes = genesOf(col), t7 = t7For(col), sig = H.signatures.files[col];
  let html = `<div class="hcard"><h5>${esc(col)} · ${genes.length} genes</h5>`;
  if (sig) html += `<div class="mini"><b>2 ·</b> see the whole domain at once</div><div class="row" style="margin-top:4px"><button class="btn sm" data-sig="signature"><i class="fa-solid fa-wave-square"></i> signature</button><button class="btn sm" data-sig="contrast"><i class="fa-solid fa-circle-half-stroke"></i> contrast</button><button class="btn sm" id="dom2all"><i class="fa-solid fa-dna"></i> these genes in All genes</button></div>
    <div class="mini" style="margin-top:3px">signature = mean ISH of the ${sig.markers} marker genes · contrast = minus the ${sig.absent_reference} genes the paper marks in sibling domains only</div>`;
  html += '</div>';
  html += `<div class="hcard"><div class="mini"><b>3 ·</b> or look at single genes (click to add, up to 4) · ★ named in a paper figure · * as published</div>
    <input class="search" id="hgfilter" placeholder="filter ${genes.length} genes…" style="margin-top:6px"><div id="hglist" class="glist"></div></div>`;
  if (t7.length) {
    html += '<div class="hcard"><div class="mini" style="margin-bottom:3px">connectivity · Supplementary Table 7</div>';
    for (const d of t7) {
      const gl = Object.entries(d.genes).map(([h, v]) => `${esc(h.replace(' Genes', ''))} ${v.length}`).join(', ');
      html += `<details class="t7"><summary><b>${esc(d.name)}</b> <span class="mini">${d.inputs.length} inputs · ${d.outputs.length} outputs · genes: ${gl}${d.shares_genes_with ? ` (list shared with ${esc(d.shares_genes_with)})` : ''}</span></summary>
        <div class="mini"><b>inputs</b>: ${d.inputs.map(esc).join(', ')}</div><div class="mini"><b>outputs</b>: ${d.outputs.map(esc).join(', ')}</div></details>`;
    }
    html += '</div>';
  }
  $('hgeadomain').innerHTML = html;
  $('hgeadomain').querySelectorAll('[data-sig]').forEach((b) => b.addEventListener('click', () => api.showSignature(col, b.dataset.sig)));
  const d2a = $('dom2all'); if (d2a) d2a.addEventListener('click', async () => { openTab('all'); await ensureAll(); setOnly(genesOf(col).filter((g) => g.ish).map((g) => g.ish)); });
  const list = () => {
    const q = ($('hgfilter').value || '').toLowerCase();
    $('hglist').innerHTML = genes.filter((g) => !q || g.name.toLowerCase().includes(q) || (g.ish || '').toLowerCase().includes(q))
      .map((g) => `<label class="${g.ish ? '' : 'dis'}" data-ish="${esc(g.ish || '')}" title="annotated in: ${esc(g.domains.map((j) => H.columns[j]).join(', '))}">
        ${g.figures.length ? '★' : ''}<b>${esc(g.name)}</b><span class="tags">${g.ish ? (g.match === 'exact' ? '' : esc(g.match.startsWith('alias') ? `ISH: ${g.ish}` : g.match)) : 'no coronal ISH'} · ${g.domains.length} dom.</span></label>`).join('');
    $('hglist').querySelectorAll('label').forEach((r) => r.addEventListener('click', () => r.dataset.ish ? api.selectGene(r.dataset.ish)
      : api.status(`${r.querySelector('b').textContent}: Allen imaged this gene sagittally only, so there is no Blue Brain coronal volume`)));
  };
  $('hgfilter').addEventListener('input', list);
  list();
}

// ───────────── Paper tab ─────────────
function renderSource() {   // Table 2 and gene coverage: always visible above the tabs
  $('hgeamatrix').addEventListener('click', openMatrix);
  const al = H.genes.filter((g) => g.match.startsWith('alias'));
  const no = H.genes.filter((g) => !g.ish);
  $('hgeanotes').innerHTML = `<summary><b>${H.genes.filter((g) => g.ish).length} of ${H.genes.length}</b> genes have Blue Brain ISH · missing: ${no.map((g) => esc(g.name)).join(', ')}</summary>
    ${H.genes.filter((g) => g.match === 'exact').length} exact, ${H.genes.filter((g) => g.match === 'asterisk removed').length} with the asterisk removed, ${al.length} under a current name:
    ${al.map((g) => `${esc(g.name)} → ${esc(g.ish)}`).join(', ')}. No coronal volume for ${no.map((g) => esc(g.name)).join(', ')} (Allen imaged them sagittally only).
    ★ = named in the paper's figure legends.<br>${esc(H.asterisk_note)}<br>Supplementary Table 7 lists genes, inputs and outputs for ${H.table7.length} domains; its gene lists match Table 2 exactly.`;
}

function openMatrix() {
  const ov = $('hmatrix');
  if (ov.style.display === 'flex') { ov.style.display = 'none'; return; }
  const cols = H.columns;
  const first = new Set(H.groups.map((g) => g.columns[0]));
  const gb = (c) => (first.has(c) ? ' gb' : '');
  const colgroup = `<colgroup><col style="width:150px">${cols.map(() => '<col>').join('')}<col style="width:120px"></colgroup>`;
  const head = `<tr><th class="grp gname" rowspan="2" style="text-align:left;vertical-align:bottom">gene (as published)</th>` +
    H.groups.map((g) => (g.columns.length === 1    // single-column group (CA2, interneurons): one header cell spanning both rows
      ? `<th class="vc gb solo" rowspan="2" data-col="${esc(g.columns[0])}" title="${esc(g.columns[0])} · open this domain">${esc(g.columns[0] === 'Putative HPF Interneurons' ? 'Int.' : g.columns[0])}</th>`
      : `<th class="grp gb" colspan="${g.columns.length}">${esc(g.group)}</th>`)).join('') + `<th class="grp gb" rowspan="2" style="vertical-align:bottom">ISH volume</th></tr>` +
    `<tr>${cols.filter((c) => !H.groups.some((g) => g.columns.length === 1 && g.columns[0] === c)).map((c) => `<th class="vc${gb(c)}" data-col="${esc(c)}" title="open the ${esc(c)} domain">${esc(short(c))}</th>`).join('')}</tr>`;
  const rows = H.genes.map((g) => `<tr data-ish="${esc(g.ish || '')}" data-name="${esc(g.name.toLowerCase())}"><td class="gname">${g.figures.length ? '★ ' : ''}${esc(g.name)}</td>` +
    cols.map((c, k) => `<td class="${(g.domains.includes(k) ? 'x' : '') + gb(c)}"></td>`).join('') + `<td class="ishc gb" title="${esc(g.ish || 'no coronal ISH volume')}">${g.ish ? (g.match === 'exact' ? '✓' : esc(g.ish)) : '— none'}</td></tr>`).join('');
  ov.innerHTML = `<div class="mh"><b>Supplementary Table 2</b> · ${H.genes.length} genes × ${cols.length} columns, exactly as published (filled = annotated). Click a gene to view it, a column to open that domain.
    <input class="search" id="mfilter" placeholder="filter genes…" style="width:200px;margin-left:10px"><button class="btn sm" id="mclose" style="float:right">close</button></div>
    <div class="mwrap"><table class="mtab">${colgroup}${head}${rows}</table></div>`;
  ov.style.display = 'flex';
  $('mclose').addEventListener('click', () => { ov.style.display = 'none'; });
  $('mfilter').addEventListener('input', () => { const q = $('mfilter').value.toLowerCase(); ov.querySelectorAll('tr[data-name]').forEach((r) => { r.style.display = !q || r.dataset.name.includes(q) ? '' : 'none'; }); });
  ov.querySelectorAll('tr[data-ish]').forEach((r) => r.addEventListener('click', () => { if (r.dataset.ish) api.selectGene(r.dataset.ish); }));
  ov.querySelectorAll('th[data-col]').forEach((h) => h.addEventListener('click', () => { openTab('domains'); selectDomain(h.dataset.col); ov.style.display = 'none'; }));
}

// ───────────── All genes tab: every HGEA gene together, each toggleable ─────────────
const ALL = { stack: null, on: null, mode: 'count', thr: 30, mingenes: 1, stats: null, shown: false, vol: null, lookup: null, timer: null };

function ensureAll() { return (ALL.loading ||= loadAll()); }   // load the stack once, however often the tab is opened

async function loadAll() {
  const s = H.stack50;
  const [vox, val] = await Promise.all([IO.array(IO.data(`hgea/${s.voxels.file}`), 'uint16'), IO.array(IO.data(`hgea/${s.values.file}`), 'uint8')]);
  ALL.stack = { vox, val, n: s.n_voxels, genes: s.genes, origin: s.origin_vox25, shape: s.shape50 };
  ALL.on = new Uint8Array(s.genes.length).fill(1);
  ALL.lookup = new Int32Array(s.shape50[0] * s.shape50[1] * s.shape50[2]).fill(-1);
  for (let j = 0; j < ALL.stack.n; j++) ALL.lookup[(vox[j * 3] * s.shape50[1] + vox[j * 3 + 1]) * s.shape50[2] + vox[j * 3 + 2]] = j;
  renderAll();
}

const byIsh = () => { const m = new Map(); H.genes.forEach((g) => g.ish && m.set(g.ish, g)); return m; };

function renderAll() {
  const G = ALL.stack.genes, gi = byIsh();
  $('hgeaall').innerHTML = `
    <div class="hcard"><div class="mini"><b>1 ·</b> combine the genes that are switched on below</div>
      <div class="row"><label class="mini" style="display:flex;gap:6px;align-items:center"><input type="checkbox" id="allshow"> show all-genes layer</label></div>
      <div class="row"><label class="mini">combine</label><select class="sel" id="allmode" style="flex:1">
        <option value="count">how many genes are expressed here</option><option value="max">strongest gene (max)</option><option value="mean">average of the genes (mean)</option></select></div>
      <div class="row" id="allthrrow"><label class="mini">expressed if ≥</label><input type="range" id="allthr" min="5" max="90" value="${ALL.thr}"><span class="mini" id="allthrv">${ALL.thr}% of p99.5</span></div>
      <div class="row" id="allminrow"><label class="mini">show where ≥</label><input type="range" id="allmin" min="1" max="60" value="${ALL.mingenes}"><span class="mini" id="allminv">${ALL.mingenes} gene${ALL.mingenes > 1 ? 's' : ''}</span></div>
      <div class="mini" id="allstats" style="margin-top:4px;color:#0f6e5c"></div>
      <div class="mini" style="margin-top:2px">50 µm, inside the hippocampus; hover a slice to see the strongest genes at a spot</div></div>
    <div class="hcard"><div class="mini"><b>2 ·</b> switch genes on / off — click a domain to show only its genes, click more domains to add them, click again to remove; or tick genes one by one</div>
      <div class="row"><button class="btn sm" id="allon">all on</button><button class="btn sm" id="alloff">all off</button><button class="btn sm" id="allinv">invert</button><span class="mini" id="allcount" style="margin-left:auto"></span></div>
      <div id="alldoms" style="margin-top:4px">${domainGrid('hdom adom')}</div>
      <input class="search" id="allfilter" placeholder="filter genes…" style="margin-top:6px"><div class="glist" id="allgenes"></div></div>`;
  $('allgenes').innerHTML = H.genes.map((g) => {
    const k = g.ish ? G.indexOf(g.ish) : -1;
    return `<label class="${k < 0 ? 'dis' : ''}" data-k="${k}" data-name="${esc(g.name.toLowerCase())}"><input type="checkbox" ${k < 0 ? 'disabled' : 'checked'}>${g.figures.length ? '★' : ''}${esc(g.name)}<span class="tags">${esc(g.domains.map((j) => short(H.columns[j])).slice(0, 4).join(' ') + (g.domains.length > 4 ? ` +${g.domains.length - 4}` : ''))}${k < 0 ? ' · no coronal ISH' : ''}</span></label>`;
  }).join('');
  $('allgenes').querySelectorAll('input[type=checkbox]').forEach((cb) => cb.addEventListener('change', () => { const k = +cb.parentElement.dataset.k; ALL.on[k] = cb.checked ? 1 : 0; changed(); }));
  $('allfilter').addEventListener('input', () => { const q = $('allfilter').value.toLowerCase(); $('allgenes').querySelectorAll('label').forEach((l) => { l.style.display = !q || l.dataset.name.includes(q) ? '' : 'none'; }); });
  $('allon').addEventListener('click', () => { ALL.on.fill(1); changed(); });
  $('alloff').addEventListener('click', () => { ALL.on.fill(0); changed(); });
  $('allinv').addEventListener('click', () => { ALL.on.forEach((v, k) => { ALL.on[k] = v ? 0 : 1; }); changed(); });
  $('alldoms').querySelectorAll('.adom').forEach((b) => b.addEventListener('click', () => {
    const ks = genesOf(b.dataset.col).filter((g) => g.ish).map((g) => G.indexOf(g.ish));
    if (ALL.on.every((v) => v)) ALL.on.fill(0);   // from "everything on", a domain click means "only this domain"
    const allOn = ks.every((k) => ALL.on[k]);
    ks.forEach((k) => { ALL.on[k] = allOn ? 0 : 1; });
    changed();
  }));
  $('allshow').addEventListener('change', () => { ALL.shown = $('allshow').checked; ALL.shown ? update() : api.removeLayer('hgea:all'); });
  $('allmode').addEventListener('change', () => { ALL.mode = $('allmode').value; $('allthrrow').style.display = $('allminrow').style.display = ALL.mode === 'count' ? 'flex' : 'none'; changed(); });
  $('allmin').addEventListener('input', () => { ALL.mingenes = +$('allmin').value; $('allminv').textContent = `${ALL.mingenes} gene${ALL.mingenes > 1 ? 's' : ''}`; changed(); });
  $('allthr').addEventListener('input', () => { ALL.thr = +$('allthr').value; $('allthrv').textContent = `${ALL.thr}% of p99.5`; changed(); });
  void gi;
  syncAll();
}

function setOnly(ishNames) {
  const G = ALL.stack.genes; ALL.on.fill(0);
  ishNames.forEach((n) => { const k = G.indexOf(n); if (k >= 0) ALL.on[k] = 1; });
  changed();
}

function syncAll() {
  const G = ALL.stack.genes, n = ALL.on.reduce((a, b) => a + b, 0);
  $('allcount').textContent = `${n} / ${G.length} on`;
  $('allgenes').querySelectorAll('label').forEach((l) => { const k = +l.dataset.k; if (k >= 0) l.querySelector('input').checked = !!ALL.on[k]; });
  $('alldoms').querySelectorAll('.adom').forEach((b) => {
    const ks = genesOf(b.dataset.col).filter((g) => g.ish).map((g) => G.indexOf(g.ish));
    const c = ks.filter((k) => ALL.on[k]).length;
    b.classList.toggle('on', c === ks.length && c > 0); b.classList.toggle('part', c > 0 && c < ks.length);
  });
}

function showAll() { ALL.shown = true; $('allshow').checked = true; }

function changed() { showAll(); syncAll(); clearTimeout(ALL.timer); ALL.timer = setTimeout(update, 200); }   // any toggle shows the layer

function composite() {
  const { val, n } = ALL.stack, sel = [];
  ALL.on.forEach((v, k) => v && sel.push(k));
  const out = new Float32Array(n), thr = ALL.thr / 100 * 255;
  if (!sel.length) return { out, sel };
  for (const k of sel) {
    const row = val.subarray(k * n, (k + 1) * n);
    if (ALL.mode === 'count') { for (let j = 0; j < n; j++) if (row[j] >= thr) out[j] += 1; }
    else if (ALL.mode === 'max') { for (let j = 0; j < n; j++) if (row[j] > out[j]) out[j] = row[j]; }
    else for (let j = 0; j < n; j++) out[j] += row[j];
  }
  if (ALL.mode === 'mean') for (let j = 0; j < n; j++) out[j] /= sel.length;
  return { out, sel };
}

function update() {
  const { out, sel } = composite(), { vox, origin, n, val } = ALL.stack, [NA, ND, NL] = IO.GRID;
  if (!ALL.vol) ALL.vol = new Uint8Array(NA * ND * NL); else ALL.vol.fill(0);
  const count = ALL.mode === 'count';
  const scale = count ? Math.min(1, 255 / Math.max(1, sel.length)) : 1;   // counts are stored as-is (≤255 genes on) so 1 gene = 1
  let top = 1, covered = 0;
  for (let j = 0; j < n; j++) { if (out[j] > top) top = out[j]; if (out[j] >= 1) covered++; }
  for (let j = 0; j < n; j++) {
    const v = Math.min(255, Math.max(out[j] > 0 ? 1 : 0, Math.round(out[j] * scale))); if (!v) continue;
    const a0 = origin[0] + vox[j * 3] * 2, d0 = origin[1] + vox[j * 3 + 1] * 2, l0 = origin[2] + vox[j * 3 + 2] * 2;
    for (let da = 0; da < 2; da++) for (let dd = 0; dd < 2; dd++) {
      const base = ((a0 + da) * ND + d0 + dd) * NL + l0; ALL.vol[base] = v; ALL.vol[base + 1] = v;
    }
  }
  if (count) {   // every switched-on gene must reach the layer: report how many are expressed anywhere
    const thr = ALL.thr / 100 * 255; let present = 0;
    for (const k of sel) { const row = val.subarray(k * n, (k + 1) * n); for (let j = 0; j < n; j++) if (row[j] >= thr) { present++; break; } }
    $('allstats').textContent = sel.length ? `${present} of ${sel.length} switched-on genes are expressed (≥ ${ALL.thr}%) somewhere in the hippocampus · ` +
      `${(100 * covered / n).toFixed(0)}% of hippocampal voxels have ≥ 1 gene · up to ${Math.round(top)} gene${Math.round(top) > 1 ? 's' : ''} at one spot` : 'no genes switched on';
  } else $('allstats').textContent = `${sel.length} genes combined (${ALL.mode})`;
  const topc = Math.max(ALL.mingenes + 1, Math.round(top * scale));
  const label = count ? `HGEA genes expressed (${sel.length} on)` : `HGEA ${ALL.mode} of ${sel.length} genes`;
  const legend = count ? { lo: ALL.mingenes, hi: Math.round(topc / scale), unit: 'genes ≥ ' + ALL.thr + '%' } : { lo: 0, hi: 1, unit: 'ISH / p99.5' };
  const extra = count ? { lo: Math.max(0.5, (ALL.mingenes - 1) * scale + 0.5), thr3d: Math.max(1, Math.round(ALL.mingenes * scale)), p995: topc } : { p995: Math.max(1, top) };
  api.upsertLayer({ key: 'hgea:all', label, vol: ALL.vol, tip: tipAt, legend, ...extra });
}

function tipAt(i25) {        // strongest switched-on genes at a 25 µm voxel
  const { origin, shape, val, n, genes } = ALL.stack, [, ND, NL] = IO.GRID;
  const a = Math.floor(i25 / (ND * NL)), d = Math.floor((i25 % (ND * NL)) / NL), l = i25 % NL;
  const a2 = (a - origin[0]) >> 1, d2 = (d - origin[1]) >> 1, l2 = (l - origin[2]) >> 1;
  if (a2 < 0 || d2 < 0 || l2 < 0 || a2 >= shape[0] || d2 >= shape[1] || l2 >= shape[2]) return null;
  const j = ALL.lookup[(a2 * shape[1] + d2) * shape[2] + l2]; if (j < 0) return null;
  const vals = []; const thr = ALL.thr / 100 * 255; let cnt = 0, on = 0;
  ALL.on.forEach((v, k) => { if (!v) return; on++; const x = val[k * n + j]; if (x >= thr) cnt++; vals.push([x, k]); });
  vals.sort((x, y) => y[0] - x[0]);
  const gi = byIsh();
  return `HGEA genes here: ${cnt} of ${on} ≥ ${ALL.thr}%<br>&nbsp;top: ` + vals.slice(0, 6).map(([x, k]) => `${esc((gi.get(genes[k]) || {}).name || genes[k])} ${Math.round(x / 2.55)}%`).join(', ');
}

export async function loadSignature(col, kind) {
  const f = H.signatures.files[col][kind], lo = H.signatures.box_lo, hi = H.signatures.box_hi;
  const crop = await IO.array(IO.data(`hgea/signatures/${f.file}`), 'uint8');
  const [NA, ND, NL] = IO.GRID, full = new Uint8Array(NA * ND * NL), w = hi[2] - lo[2], h = hi[1] - lo[1];
  for (let a = lo[0]; a < hi[0]; a++) for (let d = lo[1]; d < hi[1]; d++) {
    const src = ((a - lo[0]) * h + (d - lo[1])) * w;
    full.set(crop.subarray(src, src + w), (a * ND + d) * NL + lo[2]);
  }
  const p = H.signatures.files[col][kind === 'signature' ? 'p995_sig' : 'p995_con'] * 255;
  return { vol: full, p995: Math.max(1, p) };
}

export function hgeaBox() { return H ? [H.signatures.box_lo, H.signatures.box_hi] : null; }

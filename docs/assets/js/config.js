// Where the data lives.
// - DATA: small baked files served with the site (GitHub Pages): reference space, cells, cell types, region tables.
// - EXPR: per-cell MERFISH expression (too large for Pages) on a Hugging Face dataset; when the site runs
//   locally (python -m http.server from the repo root) the staged copy in /hf_dataset/ is used instead.
// - ISH volumes stream straight from the original EBRAINS bucket (see data/ish/genes.json url_template).
// Override the expression host with ?expr=<base-url> (remembered in localStorage).
const HF_DATASET = 'ybb2297/ccf-ish-merfish-data';     // set to your Hugging Face dataset repo

const here = new URL('.', import.meta.url);
const siteRoot = new URL('../../', here);                       // docs/
const isLocal = ['localhost', '127.0.0.1', ''].includes(location.hostname);

function exprBase() {
  const q = new URLSearchParams(location.search).get('expr');
  try { if (q) localStorage.setItem('ccfatlas.expr', q); } catch (e) { /* storage blocked */ }
  let saved = null;
  try { saved = localStorage.getItem('ccfatlas.expr'); } catch (e) { /* storage blocked */ }
  if (saved) return saved.endsWith('/') ? saved : saved + '/';
  if (isLocal) return new URL('../hf_dataset/', siteRoot).href;   // repo-root dev server
  return `https://huggingface.co/datasets/${HF_DATASET}/resolve/main/`;
}

export const CONFIG = {
  DATA: new URL('data/', siteRoot).href,
  EXPR: exprBase(),
  HF_DATASET,
};

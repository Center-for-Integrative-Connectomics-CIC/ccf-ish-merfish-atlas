// Where the data lives.
// - DATA: small baked files served with the site (GitHub Pages): reference space, cells, cell types, region tables.
// - EXPR: per-cell MERFISH expression (too large for Pages) on a Hugging Face dataset.
// - ISH volumes stream straight from the original EBRAINS bucket (see data/ish/genes.json url_template).
// Override the expression host with ?expr=<base-url>, ?expr=local (a locally built hf_dataset/ next to docs/)
// or ?expr=hf (back to the default); the choice is remembered in localStorage.
const HF_DATASET = 'ybb2297/ccf-ish-merfish-data';     // Hugging Face dataset repo

const here = new URL('.', import.meta.url);
const siteRoot = new URL('../../', here);                       // docs/
const HF_BASE = `https://huggingface.co/datasets/${HF_DATASET}/resolve/main/`;

function exprBase() {
  const q = new URLSearchParams(location.search).get('expr');
  try {
    if (q === 'hf') localStorage.removeItem('ccfatlas.expr');
    else if (q) localStorage.setItem('ccfatlas.expr', q);
  } catch (e) { /* storage blocked */ }
  let saved = null;
  try { saved = localStorage.getItem('ccfatlas.expr'); } catch (e) { /* storage blocked */ }
  if (saved === 'local') return new URL('../hf_dataset/', siteRoot).href;
  if (saved) return saved.endsWith('/') ? saved : saved + '/';
  return HF_BASE;
}

export const CONFIG = {
  DATA: new URL('data/', siteRoot).href,
  EXPR: exprBase(),
  HF_DATASET,
};

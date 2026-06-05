# E156-PROTOCOL — 786-MIII (Multilevel & Multivariate Synthesis)

- **Project:** Mulitlevelhtmlfinal (GitHub repo `Mulitlevelhtmlfinal`, user `mahmood726-cyber`)
- **Revived:** 2026-06-05 (from a single-file `786MIIImulti.html` "Gold Master" dump)
- **Type:** single-file offline browser tool + Node-testable engine
- **Dashboard:** GitHub Pages (`index.html`)

## What changed in the revival

- Made **fully offline**: vendored Plotly 2.24.1 locally (`plotly.min.js`,
  was a `cdn.plot.ly` `<script>`); removed the Google Fonts `<link>` and the two
  `preconnect` hints (system fonts fall back). No `https?://` `src`/`href` or
  remote `@import` remains.
- Extracted the statistical core (`Matrix` / `Metric` / `Engine`) **verbatim**
  into a pure `engine.js` (single source of truth; the inline duplicates were
  removed and the page now loads `engine.js`).
- Added `tests.js` (66 hand-derived assertions, all passing).
- No correctness bug was found; the math (log-scale pooling, IGLS τ²
  partition, RVE sandwich, Cheung I²) is preserved unchanged and locked by tests.
- Added Pages scaffold (`.nojekyll`, README, `.gitignore`); renamed
  `786MIIImulti.html` → `index.html`.

## Body (E156 draft — CURRENT BODY)

When effect sizes are nested within studies, does ignoring that dependence
distort a meta-analytic summary, and how is the heterogeneity actually split? This
offline dashboard ingests binary, continuous, proportion, or survival data and
forms log-scale effects with appropriate variances and small-sample corrections. It fits a three-level random-effects model by iterative generalised
least squares over a block-diagonal covariance, adds a robust variance sandwich,
and decomposes heterogeneity into between- and within-study I² following Cheung.
Across the bundled examples the multilevel fit recovers the inverse-variance
result when dependence is absent yet widens intervals once effects cluster within
studies. A 66-assertion Node suite, every value
hand-derived, confirms log-scale pooling, matrix inversion, the I² partition, the
k=1 fixed-effect degrade, and the single-cluster robust-variance limit. The honest
read is that dependence changes the uncertainty more than the point estimate, so
reporting both τ² components matters. The tool is a transparent synthesis aid for
that decomposition, not a confirmatory replacement for REML software.

SUBMITTED: [ ]

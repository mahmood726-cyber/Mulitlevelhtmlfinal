# 786-MIII — Multilevel & Multivariate Synthesis

> **⚠️ Archived duplicate.** This repository is a byte-identical copy of the
> canonical repo **[`Metamvhtml`](https://github.com/mahmood726-cyber/Metamvhtml)**
> (live app: <https://mahmood726-cyber.github.io/Metamvhtml/>). It is kept
> read-only for history; please use the canonical repo for any updates.

A single-file, **fully offline** meta-analysis dashboard that fits a
**three-level (multilevel) model** to dependent effect sizes — multiple effects
nested within studies — using an **iterative GLS (IGLS)** solver with a
**robust variance estimate (RVE) sandwich**, **Cheung's I² decomposition**
(between- vs within-cluster), and a **multilevel Egger's** small-study test. It
also provides forest/funnel/cumulative/leave-one-out plots, a power panel, and a
fragility search.

**Live app:** open `index.html` (or the GitHub Pages link). No build step, no
network, no external CDN — Plotly is vendored locally.

## Layout

```
index.html      single-file UI (loads engine.js + the local plotly.min.js)
engine.js       pure statistical core — runs in Node and the browser
tests.js        Node test harness, 66 assertions (all hand-derived)
plotly.min.js   vendored Plotly 2.24.1 (offline)
LICENSE         Apache-2.0
```

## Statistical core (`engine.js`)

| Symbol | What it does |
|---|---|
| `Matrix` | dense matrix algebra (zeros / identity / transpose / dot / Gauss-Jordan `inv` / scalarMult / add / sub) |
| `Metric.isLogRatio(m)` | which metrics are pooled on the **log scale** (`RR`, `OR`, `HR`) |
| `Engine.calcEffects(rows,type,met)` | per-row effect size + sampling variance; log scale for ratios and logit for proportions; Hedges' *g* small-sample correction; conditional 0.5 continuity correction only when a cell is zero/full |
| `Engine.pool(eff,mod,…)` | multilevel **IGLS** solver over a block-diagonal Σ (`τ²_between` shared within a study, `τ²_within` on the diagonal), an **RVE sandwich** SE with Tipton's `√(m/(m−1))` small-sample correction, **Cheung's I²** partition, and back-transformed model + robust CIs |
| `Engine.multilevelEgger(eff)` | precision-weighted OLS Egger's intercept test (guarded for k<3) |
| `Engine.getP(z)` | two-sided normal-tail p (Zelen-Severo approximation) |
| `Engine.getT(df)` | small-k t-multiplier `1.96 + 2.37/df` (1.96 floor for df≤0) |

## Fixes applied during revival (2026-06-05)

- **Made fully offline.** Plotly was loaded from `https://cdn.plot.ly/...`; it is
  now vendored as a local `plotly.min.js` (2.24.1). The Google Fonts `<link>` and
  the two `preconnect` hints were removed (system fonts fall back via the CSS
  `serif` / `sans-serif` declarations). `grep` for any `https?://` `src`/`href`
  or remote `@import` now returns nothing.
- **Single source of truth for the math.** The `Matrix` / `Metric` / `Engine`
  objects were extracted **verbatim** from the inline script into `engine.js`,
  loaded before the app script; the inline duplicates were deleted so the browser
  and the Node tests run identical code.
- **Added a 66-assertion Node test suite** (`tests.js`), every expected value
  derived by hand independently of the engine.
- Renamed `786MIIImulti.html` → `index.html` and added the Pages scaffold
  (`.nojekyll`, `.gitignore`, this README).

No correctness bug was found in the statistical core during this revival: the
pooling is on the log scale and back-transformed, `getP` is a two-sided p (≈1.0
at z=0, ≈0.05 at z=1.96), `Matrix.inv·A = I` holds, the Cheung I² components sum
to the total, k=1 degrades to a finite fixed-effect estimate, and identical
inputs give τ²=0 / I²=0. The math is preserved unchanged and locked by the tests.

## Tests

```
node tests.js
# 66 passed, 0 failed
```

Checks include a fully hand-worked fixed-effect IGLS pooling example
(2 effects in 1 cluster → pooled MD 5.4705882, model SE 0.6507914), a
two-identical-cluster log-OR case (τ²=0, I²=0, back-transform = exp(logOR)),
matrix `inv·A = I`, the Cheung I²_L2 + I²_L3 = I²_total coherence check, the k=1
degrade-to-FE passthrough, the documented RVE-NaN-at-m=1-cluster contract,
empty-input guards, and the conditional continuity-correction rule.

## Caveats

- The IGLS variance estimator is a moment-based three-level scheme; it is a
  transparent teaching/exploration implementation, not a substitute for
  `metafor::rma.mv` REML for confirmatory work. τ²_between, τ²_within, I²_total
  and the RVE SE are reported alongside every estimate.
- **RVE is undefined for a single cluster** (`√(m/(m−1))` → ∞ at m=1), so the
  robust SE is reported as `NaN` when there is only one study; the model-based
  SE remains finite. This is inherent to the sandwich estimator and is left
  in place by design (asserted by the test suite).
- Egger's test is reported only for k≥3. Apache-2.0 licensed.

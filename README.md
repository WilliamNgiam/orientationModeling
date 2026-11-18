# Code accompanying Ngiam & Lee (2026)

**Ngiam, W. X. Q., & Lee, M. D.** (2026). *Model-based evidence for task-dependent representations in psychological similarity and working memory.*  
Preprint: [https://osf.io/preprints/psyarxiv/fm9vz_v2](https://osf.io/preprints/psyarxiv/fm9vz_v3)

This repository contains the MATLAB and JAGS code used to fit and compare generative psychological models of orientation representations across perceptual reproduction, visual working memory, and similarity comparison. The analyses test whether a single shared psychological representation can account for all three tasks, or whether memory and similarity require task-specific distortions relative to perception.

Empirical data are from [Tomic & Bays (2024)](https://psycnet.apa.org/record/2023-21056-001), packaged as **`data/tomicBays2024.mat`** (see `data/README.md`). Each task is implemented with a wrap-copy generative model and, for several tasks, an alternative von Mises generative model—each with its own JAGS specification (`.txt`) and MATLAB driver (`.m`).

Models write MCMC output under per-task `storage/` (gitignored), task figures under `figures/`, and paper-level analyses under `analysis/`.

## Requirements

- **MATLAB** (scripts tested with standard toolboxes; parallel chains use the Parallel Computing Toolbox when `doParallel = 1`)
- **JAGS 4.x**
- **[trinity](https://github.com/joachimvandekerckhove/trinity)** (`callbayes`, `codatable`, `grtable`, `gelmanrubin`, etc.)
- Lab plotting utilities in `supportingFiles/` (`setFigure`, `Raxes`, `findKeepChains`, `pantoneColors.mat`)

Add the repo root and `supportingFiles/` to the MATLAB path, or run drivers from their task folders as documented below.

## Data

| Location | Contents |
|----------|----------|
| `data/tomicBays2024.mat` | **Canonical** perceptual (`dp`), memory (`dm`), and similarity (`ds`) structs |
| `data/*.csv` | Raw task exports from the authors |
| `data/parse_tomicBays2024.m` | Builds `tomicBays2024.mat` from CSV (+ `ds` from a local parse_1 product) |

All task scripts load `tomicBays2024` (via aliases such as `tomicBaysPerception2024`). See `data/README.md` for the CSV angle convention and the discarded mis-parsed `tomicBays` packaging. MCMC chain counts, burn-in, thinning, and keep rules for each analysis are in [`MCMC.md`](MCMC.md).

## Repository layout

```
orientationModeling/
├── data/                    # tomicBays2024.mat and raw CSVs
├── supportingFiles/         # Shared plotting and chain-filter helpers
├── jags-vonMises/           # JAGS module: x ~ dvonmises(mu, kappa)
├── perceptualReproduction/  # Single-stimulus reproduction
├── memoryReproduction/      # Working-memory reproduction (swap + no-swap)
├── similarityComparison/    # Pairwise similarity choice (AB vs CD)
├── jonesPewsey/             # Jones–Pewsey drivers + JAGS module
├── commonRepresentation/    # Joint representation across three tasks
├── clusterRepresentation/   # Clustered deviations from perceptual μ
└── analysis/                # Paper figures and trial-level diagnostics
```

Each task folder typically contains:

- `*_jags.txt` — JAGS model definition  
- `*.m` — MATLAB driver (`preLoad` loads cached `storage/*.mat` when available)  
- `storage/` — Saved chains (gitignored; produced by fitting)  
- `figures/` — Model-specific output figures  

### Wrap-copy generative models

Orientations are modeled on a half-circle, `[0, π)`, because line stimuli are equivalent under 180° rotation. Distances, reproduction errors, and similarity comparisons must therefore respect **modulus arithmetic on the half-circle** (wrap at 0 and π), not ordinary arithmetic on the real line.

JAGS does not provide a half-circle Gaussian. The wrap-copy construction implements the required geometry for the **Gaussian (normal-error) generative model**: on each trial a latent index selects among three copies of each mental orientation, `μ`, `μ − π`, and `μ + π`. Mental samples are drawn from a Gaussian centered on the chosen copy; distances use `min(|x₁ − x₂|, π − |x₁ − x₂|)`. That discrete wrap choice is what makes normal errors behave correctly on `[0, π)`.

| Task folder | Driver | JAGS model | What it models |
|-------------|--------|------------|----------------|
| `perceptualReproduction/` | `perceptualReproduction.m` | `perceptualReproduction_jags.txt` | Reproduce one orientation |
| `memoryReproduction/` | `memoryReproduction.m` | `memoryReproduction_jags.txt` | Recall target among foils (with swap process) |
| `memoryReproduction/` | `memoryReproductionNoSwap.m` | `memoryReproductionNoSwap_jags.txt` | Same without swap errors |
| `similarityComparison/` | `similarityComparison.m` | `similarityComparison_jags.txt` | Choose more similar pair (censored decision) |
| `commonRepresentation/` | `commonRepresentation.m` | `commonRepresentation_jags.txt` | One μ map for all three tasks |
| `clusterRepresentation/` | `clusterRepresentation.m` | `clusterRepresentation_jags.txt` | Perceptual μ plus task-specific shifts (γ) |

`similarityComparison/similarityWrapChainInits.m` supplies censored trial initial values for the similarity driver.

**Cluster γ modes.** On `tomicBays2024`, cluster posteriors are multimodal in `(γ_pm, γ_ps)`. Notes and DIC guidance live in `clusterRepresentation.m` (near the Savage–Dickey block). Prefer the near-zero `γ_pm` mode for paper figures; do not pool disagreeing modes.

### Von Mises generative models (alternative specification)

These are separate generative models, not reparameterizations of the wrap-copy models. Mental variables are von Mises on the doubled circle (`[0, 2π)`), mapped to half-circle distances without latent ±π wraps. Each task has its own `*_vonMises_jags.txt` paired with a matching driver:

- `perceptualReproduction/perceptualReproduction_vonMises.m`
- `memoryReproduction/memoryReproduction_vonMises.m`, `memoryReproductionNoSwap_vonMises.m`
- `similarityComparison/similarityComparison_vonMises.m`
- `commonRepresentation/commonRepresentation_vonMises.m`

The repository includes the modified von Mises JAGS module under `jags-vonMises/` (fork of [yeagle/jags-vonmises](https://github.com/yeagle/jags-vonmises); see that folder’s README). After clone:

```bash
cd jags-vonMises && make
```

Then set `JAGS_LIBS` to `jags-vonMises/` (or `make install` system-wide) and pass `'modules', {'vonmises'}` through `callbayes` (see comments in each von Mises driver).

**Sampling limitations.** Von Mises directions in JAGS can suffer from poor MCMC exploration when the posterior wraps the circle—scalar samplers treat `[0, 2π)` as an interval with endpoints, so chains may fail to move between modes near 0 and `2π`. That issue is discussed for the yeagle module in [this Cross Validated thread](https://stats.stackexchange.com/questions/459521/jags-circular-distribution-sampling-issues). In our hands this limits reliable inference most clearly for **similarity comparison**, where each trial conditions on four mental samples under **censoring** (`dinterval`), multiplying the effect of poor circular mixing. The wrap-copy Gaussian models in this repository are the primary fits used in the paper; von Mises drivers are included as an alternative specification, not as a drop-in replacement with equal performance.

Representation plots in both model families use **boundary-aware credible intervals** on the half-circle (error bars split across 0 and π when the interval wraps).

### Jones–Pewsey generative models

Jones–Pewsey drivers, figures, and the custom JAGS module live under `jonesPewsey/` (same level as `perceptualReproduction/`). The wrap-copy circular-normal and von Mises drivers stay in their original task folders. Unlike EM fits that search `psi` on a discrete grid, this module treats `psi` as a continuous parameter (`x ~ dJonesPewsey(mu, kappa, psi)`).

```bash
cd jonesPewsey/jags-jonesPewsey && make
```

Then set `JAGS_LIBS` to `jonesPewsey/jags-jonesPewsey/` and pass `'modules', {'jonespewsey'}` through `callbayes` (see comments in each JP driver). See `jonesPewsey/README.md`.

### Stan

This repository does **not** include Stan code. We attempted Stan implementations but did not obtain workable inference for the full paper models:

- **Similarity comparison:** the wrap-copy model uses **censoring** (`dinterval` on the similarity decision). Stan does not handle this structure in the same way as JAGS. We also tried a trembling-hand variant without censoring, but Stan’s need to represent intermediate mental samples made that approach impractical at scale.
- **Memory reproduction:** the swap model relies on a **discrete latent mixture** over recall identity (`ξ`). A standard Stan implementation with mixture marginalization did not converge in long runs with heavy thinning.
- **Perceptual reproduction:** a Stan port with the built-in von Mises would likely be feasible, but we did not pursue it once similarity and memory proved problematic.

We would welcome Stan implementations from others who can address these constraints.

## Quick start

### 1. Fit a single task model

```matlab
cd /path/to/orientationModeling/perceptualReproduction
perceptualReproduction          % set preLoad = false in script to refit
```

Equivalent entry points exist in each task folder. After a successful run, expect `storage/<modelName>_tomicBays2024_jags.mat` and figures under `figures/`.

### 2. Fit joint-task models (optional standalone)

`commonRepresentation` and `clusterRepresentation` fit **all three tasks in one JAGS run**. They load `data/tomicBays2024.mat` and write their own `storage/*.mat` files—they do **not** read posteriors from the single-task fits.

```matlab
cd ../commonRepresentation
commonRepresentation

cd ../clusterRepresentation
clusterRepresentation
```

Fit the single-task models first only if you need their chains for `analysis/drawFigures.m` (e.g. `fourScatter`) or for side-by-side comparison with the joint models.

### 3. Generate cross-task figures

```matlab
cd ../analysis
drawFigures   % edit analysisList at top of script to enable analyses
```

## Analysis scripts (`analysis/`)

| Script | Purpose |
|--------|---------|
| `drawFigures.m` | Multi-panel paper / SOI figures; enable cases via `analysisList` |

**`drawFigures.m` analysis options** (uncomment in `analysisList`):

| Name | Description |
|------|-------------|
| `fourScatter` | Perceptual, memory, similarity, common representation scatters |
| `clusterSavageDickeyEntropification` | Savage–Dickey for cluster+entrop γ |
| `statisticalSummaries` | Text summaries under `analysis/results/` |
| `threeScatter` / `twoRepresentations` / `descriptiveAdequacy*` / `swapExamples` | Optional / exploratory (commented by default) |
| `clusterSavageDickey` | Plain-cluster SD (needs Mode-A keep on 2024; commented until available) |

Outputs go to `analysis/figures/` (and `analysis/results/` for summaries). Most analyses load precomputed chains from `../<task>/storage/`.

| Script | Purpose |
|--------|---------|
| `findDiscriminatingSimilarityTrials.m` | Example trials where compressed vs physical μ maps disagree |

## Known gaps (tomicBays2024)

- **commonRepresentation** on 2024: in progress / use legacy `commonRepresentation_tomicBays_jags` as a temporary placeholder in `drawFigures` only.
- **similarityComparison** on 2024: no `storage/*_tomicBays2024_*` yet (drivers default to 2024; refit needed).
- **plain clusterRepresentation** on 2024: no agreed keep set across γ modes; DIC favors near-zero `γ_pm` (see `clusterRepresentation.m`). Cluster+entrop weekend keep is available.
- **von Mises / Jones–Pewsey** on 2024: drivers updated; refits may still be needed.

## Gitignored paths

`.gitignore` excludes `**/storage/`, `**/tmp/`, mis-parsed `data/tomicBays.mat` / `parse_1.m`, transitional `data/tomicBays_v2.mat`, ad-hoc weekend/DIC runners, and `analysis/drawFigures_legacyTomicBays.m`. Clone the repo and run drivers to regenerate MCMC files locally, or copy `storage/*.mat` from a completed analysis environment.

## Citation

Please cite the preprint above and Tomic & Bays (2024) when using this code or the data packaging.

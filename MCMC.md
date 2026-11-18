# MCMC schemes

Canonical data for all analyses: **`tomicBays2024`** (`data/tomicBays2024.mat` from `data/parse_tomicBays2024.m`). Storage files are gitignored under each task’s `storage/` folder.

Engine is **JAGS** via [trinity](https://github.com/joachimvandekerckhove/trinity) `callbayes` throughout. Parallel chains use MATLAB’s Parallel Computing Toolbox when `doParallel = 1`.

Two reusable sampling styles appear below:

1. **Fixed schedule** — one `callbayes` run with fixed chains / burn-in / samples / thin; post-hoc chain keep via `findKeepChains` (or a mode-aware selector).
2. **`orientationMcmcProtocol`** (`supportingFiles/orientationMcmcProtocol.m`) — chunked sampling with wrap alignment of location parameters, lost-chain restart, and a stop rule on all monitored diagnostics. Keeps **all** chains that finish (no `findKeepChains` surgery).

Default protocol knobs (unless a driver overrides or applies a preset):

| Knob | Default |
|------|---------|
| `nChains` | 8 |
| `nBurninAdapt` (first chunk) | 7000 |
| `nBurninRestart` (later chunks / restarts) | 2000 |
| `nThin` | 10 |
| `nDrawsWanted` (joint) | 10 000 → ~1250 kept draws/chain at 8 chains |
| `nSamplesChunk` | 250 |
| `nSamplesMax` (per chain) | 2500 |
| `rHatMax` | 1.05 |
| `essMin` / `essMinTail` | 400 |

**Stop rule (protocol):** max R̂ < `rHatMax`, min ESS and min tail ESS ≥ `essMin`, joint draws ≥ `nDrawsWanted`, and ≥ 4 chains. If the per-chain cap is hit first, all chains are saved anyway and diagnostics are printed.

**Presets** (override the defaults above):

| Preset | Where used | Knobs |
|--------|------------|-------|
| `applyLightMcmcCfg` | `commonRepresentationEntropificationMCMC` | 8 chains; burn-in 4000 / restart 1500; thin 5; 6 000 joints; cap 1500/chain; R̂ 1.1; ESS 200 |
| `applyClusterFastMcmcCfg` | `clusterRepresentationMCMC`, `clusterRepresentationEntropificationMCMC` | 4 chains; burn-in 2500 / restart 800; thin 2; 4 000 joints; cap 1000/chain; R̂ 1.1; ESS 150 |

**Weekend / select waves** (ad-hoc `runWeekend*` / `supportingFiles/weekendRunJob.m` etc. are local/gitignored; shared defs and Select drivers are tracked): fixed per-wave schedule with batch growth until enough agreeing chains. Shared defaults in `supportingFiles/weekendJobDefs.m`: burn-in 3000, samples 1500, thin 3, batch 16, max 96 chains, **min keep 4 / max keep 8**, R̂ max 1.1 on job-specific fields. Mode-aware keep (`selectClusterModeChains`) for cluster jobs; R̂-field keep (`selectConvergedChains`) for common. Paper storage tags use main model names with `*_tomicBays2024_*` (no `Weekend` in filenames).

---

## Perceptual reproduction

| | Plain driver | Automated MCMC |
|--|--------------|----------------|
| **Driver** | `perceptualReproduction/perceptualReproduction.m` | `perceptualReproduction/perceptualReproductionMCMC.m` |
| **JAGS** | `perceptualReproduction_jags.txt` | `perceptualReproductionMCMC_jags.txt` |
| **Scheme** | Fixed: 12 chains, burn-in 2000, 2000 samples, thin 5 | `orientationMcmcProtocol` **full defaults** |
| **Keep / stop** | `findKeepChains(sigma)`, R̂ ≤ 1.1, ≥2 chains | Protocol stop rule; size field `sigma` |
| **Storage tag** | `perceptualReproduction_tomicBays2024_jags.mat` | `perceptualReproductionMCMC_tomicBays2024_jags.mat` |

Paper analyses use the **MCMC** product when available.

---

## Memory reproduction (with swaps)

| | Plain driver | Automated MCMC |
|--|--------------|----------------|
| **Driver** | `memoryReproduction/memoryReproduction.m` | `memoryReproduction/memoryReproductionMCMC.m` |
| **JAGS** | `memoryReproduction_jags.txt` | `memoryReproductionMCMC_jags.txt` |
| **Scheme** | Fixed: 12 / 2000 / 2000 / thin 5 | Same numeric protocol as full defaults (**inline** chunked loop; does not call `orientationMcmcProtocol`) |
| **Keep / stop** | `findKeepChains(sigma_1)`, R̂ ≤ 1.1 | Protocol stop rule; size field `sigma_1` |
| **Storage tag** | `memoryReproduction_tomicBays2024_jags.mat` | `memoryReproductionMCMC_tomicBays2024_jags.mat` |

---

## Memory reproduction (no swap)

| | Plain driver | Automated MCMC |
|--|--------------|----------------|
| **Driver** | `memoryReproduction/memoryReproductionNoSwap.m` | `memoryReproduction/memoryReproductionNoSwapMCMC.m` |
| **JAGS** | `memoryReproductionNoSwap_jags.txt` | `memoryReproductionNoSwapMCMC_jags.txt` |
| **Scheme** | Fixed: 12 / 2000 / 2000 / thin 5 | `orientationMcmcProtocol` full defaults |
| **Keep / stop** | `findKeepChains(sigma_1)`, R̂ ≤ 1.1 | Protocol stop rule |
| **Storage tag** | `memoryReproductionNoSwap_tomicBays2024_jags.mat` | `memoryReproductionNoSwapMCMC_tomicBays2024_jags.mat` |

---

## Similarity comparison

| | Plain | Entropification | Automated MCMC |
|--|-------|-----------------|----------------|
| **Driver** | `similarityComparison/similarityComparison.m` | `similarityComparison/similarityComparisonEntropification.m` | `similarityComparison/similarityComparisonMCMC.m` |
| **JAGS** | `similarityComparison_jags.txt` | `similarityComparisonEntropification_jags.txt` | `similarityComparisonMCMC_jags.txt` |
| **Data default** | `tomicBaysSimilarity2024` → `tomicBays2024` | same | same |
| **Scheme** | Fixed: 12 / 2000 / 2000 / thin 5 | `orientationMcmcProtocol` full defaults; extra diag on `w` | `orientationMcmcProtocol` full defaults |
| **Keep / stop** | `findKeepChains(sigma)`, R̂ ≤ 1.1 | Protocol stop rule; size `sigma` | Protocol stop rule |
| **Storage tag** | `similarityComparison_tomicBays2024_jags.mat` | `similarityComparisonEntropification_tomicBays2024_jags.mat` | `similarityComparisonMCMC_tomicBays2024_jags.mat` |

Similarity `ds` is unchanged across parses, so existing fits are tagged `*_tomicBays2024_*` (renamed from `*_tomicBays_*`).

Censored-trial inits: `similarityComparison/similarityWrapChainInits.m`.

---

## Common representation

| | Plain | Entropification |
|--|-------|-----------------|
| **Driver (short)** | `commonRepresentation/commonRepresentation.m` | `commonRepresentation/commonRepresentationEntropification.m` |
| **Driver (automated)** | `commonRepresentation/commonRepresentationMCMC.m` | `commonRepresentation/commonRepresentationEntropificationMCMC.m` |
| **JAGS** | `commonRepresentation_jags.txt` | `commonRepresentationEntropification_jags.txt` |
| **Short scheme** | Fixed: 4 / 1000 / 2000 / thin 1; keep `findKeepChains(sigmaS)`, R̂ ≤ 1.1 | Prefer the `*MCMC` / weekend product |
| **MCMC scheme** | Protocol **full defaults**; size `sigmaS` | Protocol + **`applyLightMcmcCfg`**; size `sigmaS` |
| **Weekend / select** | Waves: 3000 / 1500 / thin 3; batch 16 up to 96; keep 4–8 with R̂ ≤ 1.1 on `{sigmaP, sigmaM_1, sigmaM_2, sigmaS}`. **Paper keep `[1 3 5 8 10 11 12]`** (worst R̂ ≈ 1.005) | Same wave schedule; **Paper keep `[1 2 8 11 15 16]`** (worst R̂ ≈ 1.037) |
| **Storage tags** | `commonRepresentation*_tomicBays2024_jags.mat` | `commonRepresentationEntropification*_tomicBays2024_jags.mat` |

---

## Cluster representation

| | Plain | Entropification |
|--|-------|-----------------|
| **Driver (short)** | `clusterRepresentation/clusterRepresentation.m` | Prefer entrop / Select / weekend |
| **Driver (automated)** | `clusterRepresentation/clusterRepresentationMCMC.m` | `clusterRepresentation/clusterRepresentationEntropificationMCMC.m` |
| **Select driver** | `clusterRepresentationSelectMCMC.m` | `clusterRepresentationEntropificationSelectMCMC.m` |
| **JAGS** | `clusterRepresentation_jags.txt` | `clusterRepresentationEntropification_jags.txt` |
| **Short scheme** | Fixed: 8 / 1000 / 2000 / thin 5 (**no automatic keep** on 2024 — γ multimodal) | — |
| **MCMC scheme** | Protocol + **`applyClusterFastMcmcCfg`** (gamma-focused monitors) | Same fast preset; extras include `w` when monitored |
| **Select schedule** | Fixed: **20** chains, burn-in 2500, 1500 samples, thin 3; `selectClusterModeChains` min 3 / max 8, R̂ ≤ 1.1 on `{gammaPM, gammaPS, gammaPrior}` | Same; R̂ fields also include `w` |
| **Wave runners** | 3000 / 1500 / thin 3; batch 16 up to 96; keep 4–8; mode-aware select on γ fields | Same (+ `w`); **entropification keep is the current usable paper set** → `clusterRepresentationEntropification_tomicBays2024_jags.mat` |
| **Storage tags** | `clusterRepresentation*_tomicBays2024_jags.mat` / `*SelectMCMC_*` | `clusterRepresentationEntropification*_tomicBays2024_jags.mat` |

### γ modes and DIC

On `tomicBays2024`, `(γ_pm, γ_ps)` is multimodal (recurrent modes near `(0, 1.5)` and `(2, 1.5)`). Short DIC comparisons with the Trinity `dic` module favor the near-zero `γ_pm` mode. **Do not pool disagreeing modes** for Savage–Dickey or paper figures; prefer Mode A (near-zero `γ_pm`). Notes live in `clusterRepresentation.m`. Paper SD uses `loadModeAClusterChains` on `clusterRepresentation_tomicBays2024` (Mode A region + R̂ ≤ 1.1; typical keep of 5 near `(0.12, 1.48)`).

---

## Von Mises / Jones–Pewsey (not paper primary)

- **Von Mises:** fixed schedules similar to plain wrap-copy drivers (typically 12 / 2000 / 2000 / thin 5; common uses 8 / 1000 / 2000 / thin 1) with module `{'vonmises'}`. Keep via `findKeepChains` at R̂ ≈ 1.05–1.1. Drivers: `*_vonMises.m` under each task folder.
- **Jones–Pewsey:** `jonesPewsey/*JP.m` use `orientationMcmcProtocol` full defaults with module `{'jonespewsey'}`.

---

## How to choose a driver

| Goal | Use |
|------|-----|
| Quick exploratory fit / legacy figures | Plain `*.m` fixed schedule |
| Paper posterior / Savage–Dickey sample size | `*MCMC.m` → protocol (full, light, or cluster-fast as wired) |
| Hard multimodal / stubborn R̂ (cluster, common) | Weekend / Select wave runners (local) until ≥4 agreeing chains |
| Entropification weight `w` | Matching `*Entropification*` driver (protocol or weekend) |

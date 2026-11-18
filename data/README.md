CSV versions of the Tomic & Bays (2024) orientation data.

## Packaging

**`parse_tomicBays2024.m`** builds **`tomicBays2024.mat`** (`ds`, `dp`, `dm`).

Tomic’s CSVs store perception/memory **target**, **stimulus**, and **response** as `2 * orientation` on `[-π, π)`. Those columns must be divided by 2 before folding onto the half-circle:

```matlab
mod(round(x, decimals) / 2, pi)
```

Memory **nontarget_*** columns are already on `[-π/2, π/2)` — fold with `mod` only (no `/2`). Similarity angles use the same half-circle fold as that task (no `/2` step in this builder).

Rebuild with:

```matlab
cd data
parse_tomicBays2024   % writes tomicBays2024.mat
```

All model drivers load `tomicBays2024`.

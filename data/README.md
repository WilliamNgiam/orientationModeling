CSV versions of the Tomic & Bays (2024) orientation data.

## Canonical packaging

**`parse_tomicBays2024.m`** builds **`tomicBays2024.mat`** (`ds`, `dp`, `dm`):

- Similarity `ds` matches the fold used for that task (reused from a local `parse_1` product).
- Perception `dp` and memory `dm` divide target/stimulus/response by 2 before `mod(., pi)` (Tomic stored those as `2 * orientation`; memory nontargets were already on `[-pi/2, pi/2)`).

All model drivers default to `tomicBays2024`. Rebuild with:

```matlab
cd data
parse_tomicBays2024   % writes tomicBays2024.mat
```

## Legacy / mis-parsed (local only; gitignored)

| File | Role |
|------|------|
| `parse_1.m` | Uniform `mod(x, pi)` fold → `tomicBays.mat` (incorrect for perception/memory) |
| `parse_2.m` | Old name for the 2024 builder; superseded by `parse_tomicBays2024.m` |
| `tomicBays.mat` | Mis-parsed product of `parse_1` |
| `tomicBays_v2.mat` | Historical alias of `tomicBays2024.mat` (jobs now write `*_tomicBays2024_*` tags) |

Only `parse_tomicBays2024.m` and `tomicBays2024.mat` belong on GitHub.

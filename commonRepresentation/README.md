### Common representation model

The model assumes a single latent representation generates the perceptual reproduction, memory reproduction, and similarity comparison data. It uses the same decision processes as for the individual models.

`commonRepresentation_jags.txt` implements the common representation model as a graphical model in JAGS.

`commonRepresentation.m` (and `*MCMC` / weekend drivers) apply the model to data. Default data is **`tomicBays2024`** from [Tomic & Bays (2024)](https://psycnet.apa.org/record/2023-21056-001). Paper keep: `storage/commonRepresentation_tomicBays2024_jags.mat` (keep `[1 3 5 8 10 11 12]`). See [`../MCMC.md`](../MCMC.md).

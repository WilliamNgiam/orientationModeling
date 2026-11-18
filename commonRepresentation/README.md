### Common representation model

The model assumes a single latent representation generates the perceptual reproduction, memory reproduction, and similarity comparison data. It uses the same decision processes as for the individual models.

`commonRepresentation_jags.txt` implements the common representation model as a graphical model in JAGS.

`commonRepresentation.m` and `commonRepresentationMCMC.m` (plus entropification variants) apply the model to data. Default data is **`tomicBays2024`** from [Tomic & Bays (2024)](https://psycnet.apa.org/record/2023-21056-001). The `*MCMC` drivers use the chunked protocol in [`../MCMC.md`](../MCMC.md) (paper keep `[1 3 5 8 10 11 12]`).

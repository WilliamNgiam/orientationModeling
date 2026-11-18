### Similarity comparison model

On each trial, two pairs of oriented line stimuli are presented. The response is a binary decision about which of the two pairs has more similar stimuli.

`similarityComparison_jags.txt` implements the similarity comparison model as a graphical model in JAGS. There is no need for the latent mixture approach to implementing the circular normal, because ...

`similarityComparison.m` (and entrop / MCMC / von Mises variants) apply the model to data. Default data is **`tomicBays2024`** from [Tomic & Bays (2024)](https://psycnet.apa.org/record/2023-21056-001). Similarity `ds` is unchanged across parses, so storage uses `*_tomicBays2024_*` tags. See [`../MCMC.md`](../MCMC.md).

A trial-subsetting variant of the driver lives in `old/similarityComparison_withSubsetting.m` (see `old/README.md`).

We were not able to implement the model in Stan, because of the use of censoring. In an attempt to produce a Stan implementation, we explored an alternative version of the model that did not use censoring, but instead used an error-of-execution (trembling hand) approach to make the deterministic decision rule stochastic. We were again not able to implement this model scalably in Stan, because it insists on tracking the intermediate mental samples.  We would be pleased to hear from anybody who is able to provide a Stan implementation of either version of the model.

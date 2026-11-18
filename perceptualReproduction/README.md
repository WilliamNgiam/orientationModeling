### Perceptual reproduction model

On each trial an oriented line stimulus is presented and the participant uses an interface to reproduce it.

`perceptualReproduction_jags.txt` implements the perceptual reproduction model as a graphical model in JAGS.

`perceptualReproduction.m` and `perceptualReproductionMCMC.m` (plus `perceptualReproduction_vonMises.m`) apply the model to data. Default `dataList` entry is `tomicBaysPerception2024` → `../data/tomicBays2024.mat`. Empirical data: [Tomic & Bays (2024)](https://psycnet.apa.org/record/2023-21056-001). See [`../data/README.md`](../data/README.md) and [`../MCMC.md`](../MCMC.md).

The JAGS script implements the circular normal distribution using a latent variable that allows any stimulus location to be shifted by +π or −π. We attempted to use the von Mises module [here](https://github.com/yeagle/jags-vonmises) for JAGS, but believe it does not work well. It did not produce sensible results even for very simple inference problems. We would be pleased to hear from anybody who is able to provide a working von Mises implementation of the model in JAGS.

Given difficulties in producing a Stan implementation of the memory reproduction and similarity comparison models (as documented in their READMEs), we did not pursue a Stan implementation of the perceptual reproduction model, although we think it would be straightforward using the inbuilt von Mises distribution.

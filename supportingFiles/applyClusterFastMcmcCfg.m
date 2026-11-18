function cfg = applyClusterFastMcmcCfg(cfg)
%APPLYCLUSTERFASTMCMCCFG Knobs for 2 cluster fits in ~6h (~3h each).
%
% Assumes only gammaPM/gammaPS/gammaPrior are monitored. 4 chains, burn-in
% 2.5e3, thin 2, ~4k joint draws, essMin 150, Rhat 1.1, cap 1000/chain.

cfg.nChains = 4;
cfg.nBurninAdapt = 2500;
cfg.nBurninRestart = 800;
cfg.nThin = 2;
cfg.nDrawsWanted = 4e3;
cfg.nSamplesChunk = 250;
cfg.nSamplesMax = 1000;
cfg.rHatMax = 1.1;
cfg.essMin = 150;
end

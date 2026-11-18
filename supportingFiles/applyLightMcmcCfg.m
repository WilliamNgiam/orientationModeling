function cfg = applyLightMcmcCfg(cfg)
%APPLYLIGHTMCMCCFG Knobs aimed at 3 fits in ~24h wall (~8h each).
%
% Heavier than the failed 4-chain/thin-2 pass; lighter than full defaults
% (burn-in 7e3, thin 10, ESS 400, 10k joints, cap 2500).
%
% 8 chains, burn-in 4e3, thin 5, ~6k joint draws, essMin 200, Rhat 1.1,
% per-chain cap 1500.

cfg.nChains = 8;
cfg.nBurninAdapt = 4000;
cfg.nBurninRestart = 1500;
cfg.nThin = 5;
cfg.nDrawsWanted = 6e3;
cfg.nSamplesChunk = 250;
cfg.nSamplesMax = 1500;
cfg.rHatMax = 1.1;
cfg.essMin = 200;
end

%% Cluster representation + similarity entropification MCMC
%
% Same model as clusterRepresentationEntropification.
% Sampler: orientationMcmcProtocol. Reuses clusterRepresentationEntropification_jags.txt.
% Protocol knobs via applyClusterFastMcmcCfg (~3h; gamma-only monitors).
%
% dataList: 'tomicBays2024' (parse_tomicBays2024) or 'tomicBays' (parse_1).
% Storage: storage/clusterRepresentationEntropificationMCMC_<dataName>_jags.mat

clear; close all;
preLoad = false;
printFigures = false;

thisDir = fileparts(mfilename('fullpath'));
cd(thisDir);
addpath(fullfile(thisDir, '..', 'supportingFiles'));

modelDir = './';
modelName = 'clusterRepresentationEntropificationMCMC';
jagsBase = 'clusterRepresentationEntropification';
engine = 'jags';

dataList = {...
   'tomicBays2024', ...
   % 'tomicBays', ...
   };

pi = 3.1415;
load(fullfile(thisDir, '..', 'supportingFiles', 'pantoneColors.mat'), 'pantone')
fontSize = 18;
CI = [2.5 97.5];

for dataIdx = 1:numel(dataList)
   dataName = dataList{dataIdx};
   switch dataName

      case {'tomicBays', 'tomicBays2024'}
         dataDir = fullfile(thisDir, '..', 'data');
         load(fullfile(dataDir, dataName), 'dp', 'dm', 'ds');

         nStimuli = dp.nStimuli;

         yP = dp.response;
         yM = dm.response;
         yS = ds.response;

         nPTrials = length(yP);
         sP = dp.sIdx;

         nMTrials = length(yM);
         [~, ~, setSize] = unique(dm.setSize, 'stable');
         sM = [dm.tIdx dm.nIdx];
         sM(isnan(sM)) = 1;
         [~, maxPresented] = size(sM);

         nSTrials = length(yS);
         a = ds.aIdx;
         b = ds.bIdx;
         c = ds.cIdx;
         d = ds.dIdx;
   end

   % Only gammas needed for inference (same as clusterRepresentationMCMC).
   params = {'gammaPM', 'gammaPS', 'gammaPrior'};
   stimuli = dp.stimuli(:);

   data = struct(...
      'yP'          , yP          , ...
      'yM'          , yM          , ...
      'yS'          , yS          , ...
      'nStimuli'    , nStimuli    , ...
      'nPTrials'    , nPTrials    , ...
      'nMTrials'    , nMTrials    , ...
      'nSTrials'    , nSTrials    , ...
      'sP'          , sP          , ...
      'sM'          , sM          , ...
      'setSize'     , setSize     , ...
      'maxPresented', maxPresented, ...
      'a'           , a           , ...
      'b'           , b           , ...
      'c'           , c           , ...
      'd'           , d           );

   xAinit = clampHalfLineSample(ds.a(:), pi);
   xBinit = clampHalfLineSample(ds.b(:), pi);
   xCinit = clampHalfLineSample(ds.c(:), pi);
   xDinit = clampHalfLineSample(ds.d(:), pi);

   cfg = struct();
   cfg.modelFile = sprintf('%s/%s_%s.txt', modelDir, jagsBase, engine);
   cfg.modelName = modelName;
   cfg.data = data;
   cfg.params = params;
   cfg.nStimuli = nStimuli;
   cfg.stimuli = stimuli;
   cfg.makeInits = @()makeInits(stimuli, pi, nPTrials, nMTrials, nSTrials, ...
      xAinit, xBinit, xCinit, xDinit);
   cfg.lastStateInits = @(chains)lastStateInits(chains, stimuli, pi, ...
      nPTrials, nMTrials, nSTrials, xAinit, xBinit, xCinit, xDinit);
   cfg.sizeField = 'gammaPM';
   cfg.piVal = pi;
   cfg.muPrefix = 'muP';
   cfg.alignCompanions = {};
   cfg.extraDiagPrefixes = {'gammaPM', 'gammaPS', 'gammaPrior'};
   cfg.restartJitter = @(t)jitterRestart(t, pi, nPTrials, nMTrials, nSTrials, ...
      xAinit, xBinit, xCinit, xDinit);
   cfg = applyClusterFastMcmcCfg(cfg);

   fileName = sprintf('%s_%s_%s.mat', modelName, dataName, engine);

   if preLoad && isfile(fullfile('storage', fileName))
      fprintf('Loading pre-stored samples for model %s on data %s\n', modelName, dataName);
      load(fullfile('storage', fileName), 'chains', 'stats', 'diagnostics', 'info');
      cfg.chains = chains;
      if exist('stats', 'var'), cfg.stats = stats; end
      if exist('info', 'var'), cfg.info = info; end
      [chains, stats, diagnostics, info] = orientationMcmcProtocol(cfg);
   else
      tic;
      [chains, stats, diagnostics, info] = orientationMcmcProtocol(cfg);
      fprintf('%s took %f seconds!\n', upper(engine), toc);
      fprintf('Saving samples for model %s on data %s\n', modelName, dataName);
      if ~isfolder('storage')
         mkdir('storage');
      end
      save(fullfile('storage', fileName), 'chains', 'stats', 'diagnostics', 'info', '-v7.3');
   end

   if ~diagnostics.converged
      warning('clusterRepresentationEntropificationMCMC:diagnostics', ...
         'Posterior diagnostics did not pass; summaries still use all %d chains.', ...
         size(chains.gammaPM, 2));
   end

   gammaPM = codatable(chains, 'gammaPM', @mean);
   boundsPM = prctile(chains.gammaPM(:), CI);
   fprintf('Posterior mean of gammaPM is %1.3f, with 95%% CI (%1.3f, %1.3f)\n', gammaPM, boundsPM);
   gammaPS = codatable(chains, 'gammaPS', @mean);
   boundsPS = prctile(chains.gammaPS(:), CI);
   fprintf('Posterior mean of gammaPS is %1.3f, with 95%% CI (%1.3f, %1.3f)\n', gammaPS, boundsPS);
   gammaPrior = codatable(chains, 'gammaPrior', @mean);
   boundsPrior = prctile(chains.gammaPrior(:), CI);
   fprintf('Posterior mean of gammaPrior is %1.3f, with 95%% CI (%1.3f, %1.3f)\n', gammaPrior, boundsPrior);

   if printFigures
      error('clusterRepresentationEntropificationMCMC:figures', ...
         'Figures disabled for remote runs (setFigure crashes headless).');
   end
end

function init = makeInits(stimuli, piVal, nPTrials, nMTrials, nSTrials, xAinit, xBinit, xCinit, xDinit)
nStimuli = numel(stimuli);
muP = nan(3 * nStimuli, 1);
muP(2:nStimuli) = stimuli(2:end) + 0.05 * randn(nStimuli - 1, 1);
muP(2:nStimuli) = min(max(muP(2:nStimuli), 1e-4), piVal - 1e-4);
muP(2) = min(max(muP(2), 1e-3), piVal/2 - 1e-3);
delta = 0.1 * randn(nStimuli, 2);
delta = min(max(delta, -piVal/2 + 1e-4), piVal/2 - 1e-4);
init = struct( ...
   'muP', muP, ...
   'sigmaP', 0.1 + 0.2 * rand, ...
   'sigmaM', 0.2 + 0.7 * rand(1, 2), ...
   'sigmaS', 0.15 + 0.25 * rand, ...
   'w', 0.5 + 3.5 * rand, ...
   'gammaPM', 0.3 + 0.7 * rand, ...
   'gammaPS', 0.3 + 0.7 * rand, ...
   'gammaPrior', 0.3 + 0.7 * rand, ...
   'delta', delta, ...
   'omega3', jitterSimplex([0.85 0.075 0.075]), ...
   'omega6', jitterSimplex([0.55 0.09 0.09 0.09 0.09 0.09]), ...
   'nu', ones(nPTrials, 1), ...
   'xi', ones(nMTrials, 2), ...
   'nuS', ones(nSTrials, 1), ...
   'xA', xAinit, ...
   'xB', xBinit, ...
   'xC', xCinit, ...
   'xD', xDinit);
end

function inits = lastStateInits(chains, stimuli, piVal, nPTrials, nMTrials, nSTrials, xAinit, xBinit, xCinit, xDinit)
% Gammas are monitored; other nodes are re-drawn (not in coda).
nCh = size(chains.gammaPM, 2);
inits = cell(1, nCh);
for c = 1:nCh
   init = makeInits(stimuli, piVal, nPTrials, nMTrials, nSTrials, ...
      xAinit, xBinit, xCinit, xDinit);
   init.gammaPM = max(chains.gammaPM(end, c), 1e-3);
   init.gammaPS = max(chains.gammaPS(end, c), 1e-3);
   init.gammaPrior = max(chains.gammaPrior(end, c), 1e-3);
   inits{c} = init;
end
end

function t = jitterRestart(t, piVal, nPTrials, nMTrials, nSTrials, xAinit, xBinit, xCinit, xDinit)
t.sigmaP = min(max(t.sigmaP + 0.02 * randn, 0.05), piVal - 1e-3);
t.sigmaM = min(max(t.sigmaM + 0.02 * randn(size(t.sigmaM)), 0.05), piVal - 1e-3);
t.sigmaS = min(max(t.sigmaS + 0.02 * randn, 0.05), piVal - 1e-3);
t.w = min(max(t.w + 0.3 * randn, 1e-3), 19.9);
t.gammaPM = max(t.gammaPM + 0.05 * randn, 1e-3);
t.gammaPS = max(t.gammaPS + 0.05 * randn, 1e-3);
t.gammaPrior = max(t.gammaPrior + 0.05 * randn, 1e-3);
t.delta = min(max(t.delta + 0.02 * randn(size(t.delta)), -piVal/2 + 1e-4), piVal/2 - 1e-4);
t.omega3 = jitterSimplex(t.omega3);
t.omega6 = jitterSimplex(t.omega6);
t.nu = ones(nPTrials, 1);
t.xi = ones(nMTrials, 2);
t.nuS = ones(nSTrials, 1);
t.xA = xAinit;
t.xB = xBinit;
t.xC = xCinit;
t.xD = xDinit;
end

function p = jitterSimplex(p)
p = p(:)' .* exp(0.15 * randn(size(p)));
p = max(p, 1e-6);
p = p / sum(p);
end

function x = clampHalfLineSample(x, piVal)
x = mod(x, piVal);
x = min(max(x, 0), piVal - 1e-10);
end

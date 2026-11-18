%% Cluster + similarity entropification: many chains, keep good muS mode
%
% Same selection workflow as clusterRepresentationSelectMCMC.
% Reuses clusterRepresentationEntropification_jags.txt.
% Monitors: gammaPM, gammaPS, gammaPrior, muP, muM, muS, w.
%
% Storage: storage/clusterRepresentationEntropificationSelectMCMC_<data>_jags.mat
% printFigures must stay false for headless / matlab-batch runs.

clear; close all;
preLoad = true;
printFigures = false;  % never open figures remotely (setFigure crashes headless)

thisDir = fileparts(mfilename('fullpath'));
cd(thisDir);
addpath(fullfile(thisDir, '..', 'supportingFiles'));

modelDir = './';
modelName = 'clusterRepresentationEntropificationSelectMCMC';
jagsBase = 'clusterRepresentationEntropification';
engine = 'jags';

dataList = {'tomicBays2024'};

nChains    = 20;
nBurnin    = 2500;
nSamples   = 1500;
nThin      = 3;
doParallel = 1;
minKeep    = 3;
maxKeep    = 8;
rHatMax    = 1.1;

pi = 3.1415;
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

   params = {'gammaPM', 'gammaPS', 'gammaPrior', 'muP', 'muM', 'muS', 'w'};
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

   fileName = sprintf('%s_%s_%s.mat', modelName, dataName, engine);
   modelFile = sprintf('%s/%s_%s.txt', modelDir, jagsBase, engine);

   if preLoad && isfile(fullfile('storage', fileName))
      fprintf('Loading pre-stored samples for model %s on data %s\n', modelName, dataName);
      load(fullfile('storage', fileName), 'chainsAll', 'chains', 'keepChains', ...
         'modeScores', 'modeInfo', 'stats', 'diagnostics', 'info');
   else
      tic;
      fprintf(['Sampling %s: %d chains, burn-in %d, keep %d, thin %d ', ...
         '(monitors: gammas + muP/muM/muS + w)\n'], ...
         modelName, nChains, nBurnin, nSamples, nThin);
      [stats, chainsAll, diagnostics, info] = callbayes(engine, ...
         'model'           , modelFile                                 , ...
         'data'            , data                                      , ...
         'outputname'      , 'samples'                                 , ...
         'init'            , @()makeInits(stimuli, pi, nPTrials, nMTrials, nSTrials, ...
                                xAinit, xBinit, xCinit, xDinit)        , ...
         'datafilename'    , modelName                                 , ...
         'initfilename'    , modelName                                 , ...
         'scriptfilename'  , modelName                                 , ...
         'logfilename'     , sprintf('/tmp/%s', modelName)             , ...
         'nchains'         , nChains                                   , ...
         'nburnin'         , nBurnin                                   , ...
         'nsamples'        , nSamples                                  , ...
         'monitorparams'   , params                                    , ...
         'thin'            , nThin                                     , ...
         'workingdir'      , sprintf('/tmp/%s', modelName)             , ...
         'verbosity'       , 0                                         , ...
         'saveoutput'      , true                                      , ...
         'allowunderscores', 1                                         , ...
         'parallel'        , doParallel                                );
      fprintf('%s took %.1f seconds\n', upper(engine), toc);

      chainsAll = dropWrapMus(chainsAll, nStimuli, {'muP', 'muM', 'muS'});
      [chainsAll, nFlips] = alignClusterMuChains(chainsAll, stimuli, pi, ...
         {'muP', 'muM', 'muS'});
      fprintf('Wrap-aligned %d chain-by-bin shifts\n', nFlips);

      [keepChains, modeScores, modeInfo] = selectClusterModeChains( ...
         chainsAll, stimuli, pi, ...
         'minKeep', minKeep, 'maxKeep', maxKeep, 'rHatMax', rHatMax, ...
         'rHatFields', {'gammaPM', 'gammaPS', 'gammaPrior', 'w'});

      if isempty(keepChains)
         warning('clusterRepresentationEntropificationSelectMCMC:noKeep', ...
            '%s Saving chainsAll only.', modeInfo.message);
         chains = struct();
      else
         chains = subsetChains(chainsAll, keepChains);
         fprintf('Rhat on kept chains (gammas, w):\n');
         for g = {'gammaPM', 'gammaPS', 'gammaPrior', 'w'}
            r = gelmanrubin(chains.(g{1}), 0, 1, 'rhat');
            fprintf('  %s Rhat = %.3f\n', g{1}, r);
         end
      end

      if ~isfolder('storage')
         mkdir('storage');
      end
      fprintf('Saving samples for model %s on data %s\n', modelName, dataName);
      save(fullfile('storage', fileName), ...
         'chainsAll', 'chains', 'keepChains', 'modeScores', 'modeInfo', ...
         'stats', 'diagnostics', 'info', 'stimuli', 'params', '-v7.3');
   end

   if isempty(keepChains) || ~isfield(chains, 'gammaPM')
      warning('clusterRepresentationEntropificationSelectMCMC:skipSummary', ...
         'No agreed keep set; skip gamma/muS summaries.');
      continue
   end

   gammaPM = mean(chains.gammaPM(:));
   boundsPM = prctile(chains.gammaPM(:), CI);
   fprintf('Kept posterior mean gammaPM = %1.3f, 95%% CI (%1.3f, %1.3f)\n', ...
      gammaPM, boundsPM);
   gammaPS = mean(chains.gammaPS(:));
   boundsPS = prctile(chains.gammaPS(:), CI);
   fprintf('Kept posterior mean gammaPS = %1.3f, 95%% CI (%1.3f, %1.3f)\n', ...
      gammaPS, boundsPS);
   gammaPrior = mean(chains.gammaPrior(:));
   boundsPrior = prctile(chains.gammaPrior(:), CI);
   fprintf('Kept posterior mean gammaPrior = %1.3f, 95%% CI (%1.3f, %1.3f)\n', ...
      gammaPrior, boundsPrior);
   wMean = mean(chains.w(:));
   wBounds = prctile(chains.w(:), CI);
   fprintf('Kept posterior mean w = %1.3f, 95%% CI (%1.3f, %1.3f)\n', ...
      wMean, wBounds);

   muSmean = nan(nStimuli, 1);
   for j = 1:nStimuli
      muSmean(j) = circularMeanHalf(chains.(sprintf('muS_%d', j))(:), pi);
   end
   fprintf('Kept muS circular means (for paper comparison):\n');
   fprintf('  %s\n', mat2str(muSmean', 3));

   if printFigures
      error('clusterRepresentationEntropificationSelectMCMC:figures', ...
         'Figures disabled for remote runs; set printFigures only on a display.');
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

function chains = subsetChains(chainsAll, keep)
chains = struct();
fn = fieldnames(chainsAll);
for i = 1:numel(fn)
   chains.(fn{i}) = chainsAll.(fn{i})(:, keep);
end
end

function chains = dropWrapMus(chains, nStimuli, prefixes)
fn = fieldnames(chains);
for i = 1:numel(fn)
   name = fn{i};
   for p = 1:numel(prefixes)
      token = [prefixes{p} '_'];
      if numel(name) <= numel(token) || ~strncmp(name, token, numel(token))
         continue
      end
      idx = str2double(name(numel(token)+1:end));
      if ~isnan(idx) && idx == floor(idx) && idx > nStimuli
         chains = rmfield(chains, name);
      end
   end
end
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

function mu0 = circularMeanHalf(vals, piVal)
phi = 2 * vals(:);
phi0 = atan2(mean(sin(phi)), mean(cos(phi)));
mu0 = mod(phi0 / 2, piVal);
end

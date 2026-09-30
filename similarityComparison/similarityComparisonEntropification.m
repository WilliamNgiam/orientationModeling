%% Similarity comparison with Grünwald entropification
%
% Same wrap-copy mental samples as similarityComparisonMCMC. The censored
% decision is replaced by a logistic in the signed distance difference,
% with one shared weight w. Sampler: orientationMcmcProtocol (8 chains
% kept in full, data-informed inits, wrap alignment, split-Rhat / ESS,
% 10k joint draws).
%
% Storage: storage/similarityComparisonEntropification_tomicBays_jags.mat

clear; close all;
preLoad = true;
printFigures = false;

thisDir = fileparts(mfilename('fullpath'));
cd(thisDir);
addpath(fullfile(thisDir, '..', 'supportingFiles'));

modelDir = './';
modelName = 'similarityComparisonEntropification';
engine = 'jags';

dataList = {...
   'tomicBaysSimilarity', ...
   };

pi = 3.1415;
load(fullfile(thisDir, '..', 'supportingFiles', 'pantoneColors.mat'), 'pantone')
fontSize = 18;
CI = [2.5 97.5];

for dataIdx = 1:numel(dataList)
   dataName = dataList{dataIdx};
   switch dataName

      case 'tomicBaysSimilarity'
         dataDir = fullfile(thisDir, '..', 'data');
         dataName = 'tomicBays';
         load(fullfile(dataDir, dataName), 'ds');

         a = ds.aIdx;
         b = ds.bIdx;
         c = ds.cIdx;
         d = ds.dIdx;
         y = ds.response;
         nTrials = ds.nTrials;
         nStimuli = ds.nStimuli;
   end

   params = {'mu', 'sigma', 'w'};
   stimuli = ds.stimuli(:);

   data = struct(...
      'a'        , a        , ...
      'b'        , b        , ...
      'c'        , c        , ...
      'd'        , d        , ...
      'y'        , y        , ...
      'nStimuli' , nStimuli , ...
      'nTrials'  , nTrials  );

   xAinit = clampHalfLineSample(ds.a(:), pi);
   xBinit = clampHalfLineSample(ds.b(:), pi);
   xCinit = clampHalfLineSample(ds.c(:), pi);
   xDinit = clampHalfLineSample(ds.d(:), pi);
   nuInit = ones(nTrials, 1);

   cfg = struct();
   cfg.modelFile = sprintf('%s/%s_%s.txt', modelDir, modelName, engine);
   cfg.modelName = modelName;
   cfg.data = data;
   cfg.params = params;
   cfg.nStimuli = nStimuli;
   cfg.stimuli = stimuli;
   cfg.makeInits = @()makeInits(stimuli, pi, xAinit, xBinit, xCinit, xDinit, nuInit);
   cfg.lastStateInits = @(chains)lastStateInits(chains, nStimuli, pi, ...
      xAinit, xBinit, xCinit, xDinit, nuInit);
   cfg.sizeField = 'sigma';
   cfg.piVal = pi;
   cfg.extraDiagPrefixes = {'w'};
   cfg.restartJitter = @(t)jitterRestart(t, pi, xAinit, xBinit, xCinit, xDinit, nuInit);

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
      warning('similarityComparisonEntropification:diagnostics', ...
         'Posterior diagnostics did not pass; summaries still use all %d chains.', ...
         size(chains.sigma, 2));
   end

   sigma = codatable(chains, 'sigma', @mean);
   bounds = prctile(chains.sigma(:), CI);
   fprintf('Posterior mean of sigma is %1.3f, with 95%% CI (%1.3f, %1.3f)\n', sigma, bounds);

   wMean = codatable(chains, 'w', @mean);
   wBounds = prctile(chains.w(:), CI);
   fprintf('Posterior mean of w is %1.3f, with 95%% CI (%1.3f, %1.3f)\n', wMean, wBounds);

   F = figure; clf; hold on;
   setFigure(F, [0.2 0.2 0.4 0.4], '');

   mu = nan(ds.nStimuli, 1);
   muBounds = nan(ds.nStimuli, 2);
   for idx = 1:ds.nStimuli
      vals = chains.(sprintf('mu_%d', idx))(:);
      [mu(idx), muBounds(idx, :)] = summarizeHalfCircle(vals, CI);
   end
   muTruth = ds.stimuli;

   cla; hold on;
   set(gca, ...
      'xlim'       , [0 pi]    , ...
      'xtick'      , [0 pi/4 pi/2 3*pi/4 pi]   , ...
      'xticklabelrot', 0, ...
      'xticklabel' , {'$0$', '$\frac{\pi}{4}$', '$\frac{\pi}{2}$', '$\frac{3\pi}{4}$', '$\pi$'}, ...
      'ylim'       , [0 pi]    , ...
      'ytick'      , [0 pi/4 pi/2 3*pi/4 pi]   , ...
      'yticklabel' , {'$0$', '$\frac{\pi}{4}$', '$\frac{\pi}{2}$', '$\frac{3\pi}{4}$', '$\pi$'}, ...
      'ticklabelinterpreter', 'latex', ...
      'box'        , 'off'     , ...
      'tickdir'    , 'out'     , ...
      'layer'      , 'top'     , ...
      'ticklength' , [0.02 0]  , ...
      'clipping'   , 'off'     , ...
      'fontsize'   , fontSize  );
   axis square;
   ylabel('Psychological', 'fontsize', fontSize);
   xlabel('Physical', 'fontsize', fontSize);
   moveAxis(gca, [1 1 0.95 0.95], [0 0.025 0 0]);
   Raxes(gca, 0.02, 0.01);

   for i = pi/4:pi/4:3*pi/4
      plot([i i], [0 pi], '-', 'color', pantone.GlacierGray);
      plot([0 pi], [i i], '-', 'color', pantone.GlacierGray);
   end

   for idx = 1:ds.nStimuli
      if muBounds(idx, 1) > muBounds(idx, 2)
         plot(muTruth(idx) * [1 1], [0 muBounds(idx, 2)], '-', ...
            'color', pantone.ClassicBlue, 'linewidth', 1);
         plot(muTruth(idx) * [1 1], [muBounds(idx, 1) pi], '-', ...
            'color', pantone.ClassicBlue, 'linewidth', 1);
      else
         plot(muTruth(idx) * [1 1], muBounds(idx, :), '-', ...
            'color', pantone.ClassicBlue, 'linewidth', 1);
      end
      plot(muTruth(idx), mu(idx), 'o', ...
         'markerfacecolor', pantone.ClassicBlue, ...
         'markeredgecolor', 'w', 'linewidth', 0.5, 'markersize', 4);
   end
   plot([0 pi], [0 pi], '-', 'color', pantone.AuroraRed, 'linewidth', 0.5);

   if printFigures
      if ~isfolder('figures')
         mkdir('figures');
      end
      figBase = sprintf('figures/%s_%s', dataName, modelName);
      print([figBase '.png'], '-dpng');
      print([figBase '.eps'], '-depsc');
   end

end

function init = makeInits(stimuli, piVal, xAinit, xBinit, xCinit, xDinit, nuInit)
nStimuli = numel(stimuli);
mu = stimuli + 0.05 * randn(nStimuli, 1);
mu = min(max(mu, 1e-4), piVal - 1e-4);
mu(1) = 1e-7;
mu(2) = min(max(mu(2), 1e-3), piVal/2 - 1e-3);
init = struct( ...
   'mu', mu, ...
   'sigma', 0.15 + 0.25 * rand, ...
   'w', 0.5 + 3.5 * rand, ...
   'xA', xAinit, ...
   'xB', xBinit, ...
   'xC', xCinit, ...
   'xD', xDinit, ...
   'nu', nuInit);
end

function inits = lastStateInits(chains, nStimuli, piVal, xAinit, xBinit, xCinit, xDinit, nuInit)
nCh = size(chains.sigma, 2);
inits = cell(1, nCh);
for c = 1:nCh
   mu = zeros(nStimuli, 1);
   for j = 1:nStimuli
      fname = sprintf('mu_%d', j);
      if isfield(chains, fname)
         mu(j) = chains.(fname)(end, c);
      end
   end
   mu(1) = min(max(mu(1), 0), 1e-6);
   mu(2) = min(max(mu(2), 1e-3), piVal/2 - 1e-3);
   for j = 3:nStimuli
      mu(j) = min(max(mu(j), 1e-4), piVal - 1e-4);
   end
   inits{c} = struct( ...
      'mu', mu, ...
      'sigma', min(max(chains.sigma(end, c), 1e-3), piVal - 1e-3), ...
      'w', min(max(chains.w(end, c), 1e-3), 19.9), ...
      'xA', xAinit, ...
      'xB', xBinit, ...
      'xC', xCinit, ...
      'xD', xDinit, ...
      'nu', nuInit);
end
end

function t = jitterRestart(t, piVal, xAinit, xBinit, xCinit, xDinit, nuInit)
t.sigma = min(max(t.sigma + 0.02 * randn, 0.05), piVal - 1e-3);
t.w = min(max(t.w + 0.3 * randn, 1e-3), 19.9);
t.xA = xAinit;
t.xB = xBinit;
t.xC = xCinit;
t.xD = xDinit;
t.nu = nuInit;
end

function x = clampHalfLineSample(x, piVal)
x = mod(x, piVal);
x = min(max(x, 0), piVal - 1e-10);
end

function [muMean, bounds] = summarizeHalfCircle(vals, CI)
vals = vals(:);
phi = mod(2 * vals, 2 * pi);
phi0 = atan2(mean(sin(phi)), mean(cos(phi)));
if phi0 < 0
   phi0 = phi0 + 2 * pi;
end
rel = angle(exp(1i * (phi - phi0)));
relBounds = prctile(rel, CI);
muMean = mod(phi0 / 2, pi);
bounds = mod(muMean + relBounds / 2, pi);
end

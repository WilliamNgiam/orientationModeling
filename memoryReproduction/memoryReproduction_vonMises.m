%% Memory reproduction model (von Mises)
%
% Von Mises rewrite of orientationModeling/memoryReproduction.
% Uses doubled-angle responses (yCirc = 2*y) instead of latent +/- pi wraps.
%
% Requires on the MATLAB path:
%   callbayes, codatable, grtable, gelmanrubin (trinity), get_matrix_from_coda, moveAxis
%   supportingFiles: setFigure.m, Raxes.m, findKeepChains.m, pantoneColors.mat
%
% Requires the von Mises JAGS module in ../jags-vonMises (see README).

clear; close all;
preLoad = true;
printFigures = true;

thisDir = fileparts(mfilename('fullpath'));
cd(thisDir);
addpath(fullfile(thisDir, '..', 'supportingFiles'));

jagsModuleDir = fullfile(thisDir, '..', 'jags-vonMises');
jagsModuleFile = fullfile(jagsModuleDir, 'vonmises.so');
if ~isfile(jagsModuleFile)
   error(['Build the JAGS von Mises module first: cd %s && make\n' ...
          'Expected: %s'], jagsModuleDir, jagsModuleFile);
end
setenv('JAGS_LIBS', jagsModuleDir);

modelDir = './';
modelName = 'memoryReproduction_vonMises';
engine = 'jags';

dataList = {...
   'tomicBaysMemory'; ...
   };

pi = 3.141592653589793;
load(fullfile(thisDir, '..', 'supportingFiles', 'pantoneColors.mat'), 'pantone')
fontSize = 18;
CI = [2.5 97.5];

for dataIdx = 1:numel(dataList)
   dataName = dataList{dataIdx};
   switch dataName

      case 'tomicBaysMemory'
         dataDir = fullfile(thisDir, '..', 'data');
         dataName = 'tomicBays';
         load(fullfile(dataDir, dataName), 'dm');

         [~, ~, setSize] = unique(dm.setSize, 'stable');
         y = dm.response;
         yCirc = mod(2 * y, 2 * pi);
         nTrials = dm.nTrials;
         nStimuli = dm.nStimuli;
         s = [dm.tIdx dm.nIdx];
         s(isnan(s)) = 1;
         [~, maxPresented] = size(s);

         for t = 1:nTrials
            vals = dm.nontarget(t, 1:(dm.setSize(t)-1));
            [~, srt] = sort([0 min(abs(dm.target(t)-vals), pi-abs(dm.target(t)-vals))], 'ascend');
            s(t, 1:dm.setSize(t)) = s(t, srt);
         end
   end

   params = {'mu', 'sigma', 'xi', 'yCircP', 'omega3', 'omega6', 'omega3prior', 'omega6prior'};

   nChains    = 12;
   nBurnin    = 2e3;
   nSamples   = 2e3;
   nThin      = 5;
   doParallel = 1;

   data = struct(...
      'setSize'     , setSize     , ...
      's'           , s           , ...
      'yCirc'       , yCirc       , ...
      'nStimuli'    , nStimuli    , ...
      'nTrials'     , nTrials     , ...
      'maxPresented', maxPresented);

   generator = @()struct('sigma', rand(1, 2) * pi);

   fileName = sprintf('%s_%s_%s.mat', modelName, dataName, engine);

   if preLoad && isfile(fullfile('storage', fileName))
      fprintf('Loading pre-stored samples for model %s on data %s\n', modelName, dataName);
      load(fullfile('storage', fileName), 'chains', 'stats', 'diagnostics', 'info');
   else
      tic;
      [stats, chains, diagnostics, info] = callbayes(engine, ...
         'model'           , sprintf('%s/%s_%s.txt', modelDir, modelName, engine), ...
         'data'            , data                                      , ...
         'outputname'      , 'samples'                                 , ...
         'init'            , generator                                 , ...
         'datafilename'    , modelName                                 , ...
         'initfilename'    , modelName                                 , ...
         'scriptfilename'  , modelName                                 , ...
         'logfilename'     , sprintf('/tmp/%s', modelName)             , ...
         'nchains'         , nChains                                   , ...
         'nburnin'         , nBurnin                                   , ...
         'nsamples'        , nSamples                                  , ...
         'monitorparams'   , params                                    , ...
         'thin'            , nThin                                     , ...
         'modules'         , {'vonmises'}                              , ...
         'workingdir'      , sprintf('/tmp/%s', modelName)             , ...
         'verbosity'       , 0                                         , ...
         'saveoutput'      , true                                      , ...
         'allowunderscores', true                                      , ...
         'parallel'        , doParallel                                );
      fprintf('%s took %f seconds!\n', upper(engine), toc);

      disp('Convergence statistics:')
      grtable(chains, 1.05)

      disp('Descriptive statistics for all chains:')
      codatable(chains);

      fprintf('Saving samples for model %s on data %s\n', modelName, dataName);
      if ~isfolder('storage')
         mkdir('storage');
      end
      save(fullfile('storage', fileName), 'chains', 'stats', 'diagnostics', 'info', '-v7.3');
   end

   [keepChains, rHat] = findKeepChains(chains.sigma_1, 2, 1.05);
   fields = fieldnames(chains);
   for i = 1:numel(fields)
      chains.(fields{i}) = chains.(fields{i})(:, keepChains);
   end

   sigma3 = codatable(chains, 'sigma_1', @mean);
   bounds3 = prctile(chains.sigma_1(:), CI);
   fprintf('Posterior mean of sigma for set size 3 is %1.3f, with 95%% CI (%1.3f, %1.3f)\n', sigma3, bounds3);

   sigma6 = codatable(chains, 'sigma_2', @mean);
   bounds6 = prctile(chains.sigma_2(:), CI);
   fprintf('Posterior mean of sigma for set size 6 is %1.3f, with 95%% CI (%1.3f, %1.3f)\n', sigma6, bounds6);

   omega3 = get_matrix_from_coda(chains, 'omega3', @mean);
   omega6 = get_matrix_from_coda(chains, 'omega6', @mean);
   fprintf('Posterior means for omega for set size 3 are (%1.3f, %1.3f, %1.3f)\n', omega3);
   fprintf('Posterior means for omega for set size 6 are (%1.3f, %1.3f, %1.3f, %1.3f, %1.3f, %1.3f)\n', omega6);

   for i = 1:3
      bounds = prctile(chains.(sprintf('omega3_%d', i))(:), CI);
      fprintf('Posterior mean of omega_%d for set size 3 is %1.3f, with 95%% CI (%1.3f, %1.3f)\n', i, omega3(i), bounds);
   end
   for i = 1:6
      bounds = prctile(chains.(sprintf('omega6_%d', i))(:), CI);
      fprintf('Posterior mean of omega_%d for set size 6 is %1.3f, with 95%% CI (%1.3f, %1.3f)\n', i, omega6(i), bounds);
   end

   F = figure; clf; hold on;
   setFigure(F, [0.2 0.2 0.4 0.4], '');

   mu = nan(dm.nStimuli, 1);
   muBounds = nan(dm.nStimuli, 2);
   for idx = 1:dm.nStimuli
      vals = chains.(sprintf('mu_%d', idx))(:);
      [mu(idx), muBounds(idx, :)] = summarizeHalfCircle(vals, CI);
   end
   muTruth = dm.stimuli;

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

   for idx = 1:dm.nStimuli
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

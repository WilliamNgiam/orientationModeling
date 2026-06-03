%% Common representation model (von Mises)
%
% Von Mises rewrite of orientationModeling/commonRepresentation.
% Shared mu across perceptual, memory, and similarity tasks; doubled-angle
% responses for reproduction tasks and dvonmises mental samples for similarity.
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
addpath(fullfile(thisDir, '..'));

jagsModuleDir = fullfile(thisDir, '..', 'jags-vonMises');
jagsModuleFile = fullfile(jagsModuleDir, 'vonmises.so');
if ~isfile(jagsModuleFile)
   error(['Build the JAGS von Mises module first: cd %s && make\n' ...
          'Expected: %s'], jagsModuleDir, jagsModuleFile);
end
setenv('JAGS_LIBS', jagsModuleDir);

modelDir = './';
modelName = 'commonRepresentation_vonMises';
engine = 'jags';

dataList = {...
   'tomicBays'; ...
   };

pi = 3.141592653589793;
load(fullfile(thisDir, '..', 'supportingFiles', 'pantoneColors.mat'), 'pantone')
fontSize = 18;
CI = [2.5 97.5];

for dataIdx = 1:numel(dataList)
   dataName = dataList{dataIdx};
   switch dataName

      case 'tomicBays'
         dataDir = fullfile(thisDir, '..', 'data');
         dataName = 'tomicBays';
         load(fullfile(dataDir, dataName), 'dp', 'dm', 'ds');

         nStimuli = dp.nStimuli;

         yP = dp.response;
         yM = dm.response;
         yS = ds.response;

         yPCirc = mod(2 * yP, 2 * pi);
         yMCirc = mod(2 * yM, 2 * pi);

         nPTrials = length(yP);
         sP = dp.sIdx;

         nMTrials = length(yM);
         [~, ~, setSize] = unique(dm.setSize, 'stable');
         sM = [dm.tIdx dm.nIdx];
         sM(isnan(sM)) = 1;
         [~, maxPresented] = size(sM);

         for t = 1:dm.nTrials
            vals = dm.nontarget(t, 1:(dm.setSize(t)-1));
            [~, srt] = sort([0 min(abs(dm.target(t)-vals), pi-abs(dm.target(t)-vals))], 'ascend');
            sM(t, 1:dm.setSize(t)) = sM(t, srt);
         end

         nSTrials = length(yS);
         a = ds.aIdx;
         b = ds.bIdx;
         c = ds.cIdx;
         d = ds.dIdx;
   end

   params = {'mu', 'sigmaP', 'sigmaM', 'sigmaS', 'ySP', 'xi', 'omega3', 'omega6'};

   nChains    = 8;
   nBurnin    = 1e3;
   nSamples   = 2e3;
   nThin      = 1;
   doParallel = 1;

   data = struct(...
      'yPCirc'      , yPCirc    , ...
      'yMCirc'      , yMCirc    , ...
      'yS'          , yS        , ...
      'nStimuli'    , nStimuli  , ...
      'nPTrials'    , nPTrials  , ...
      'nMTrials'    , nMTrials  , ...
      'nSTrials'    , nSTrials  , ...
      'sP'          , sP        , ...
      'sM'          , sM        , ...
      'setSize'     , setSize   , ...
      'maxPresented', maxPresented, ...
      'a'           , a         , ...
      'b'           , b         , ...
      'c'           , c         , ...
      'd'           , d         );

   for t = 1:nSTrials
      if yS(t) == 0
         xAinit(t) = 0.6; xBinit(t) = 0.7;
         xCinit(t) = 0.6; xDinit(t) = 0.8;
      else
         xAinit(t) = 0.6; xBinit(t) = 0.8;
         xCinit(t) = 0.6; xDinit(t) = 0.7;
      end
   end

   generator = @()struct(...
      'sigmaP', rand * pi, ...
      'sigmaS', rand * pi, ...
      'sigmaM', rand(1, 2) * pi, ...
      'phiA', mod(2 * xAinit, 2 * pi), ...
      'phiB', mod(2 * xBinit, 2 * pi), ...
      'phiC', mod(2 * xCinit, 2 * pi), ...
      'phiD', mod(2 * xDinit, 2 * pi));

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

   [keepChains, rHat] = findKeepChains(chains.sigmaS, 2, 1.05);
   fields = fieldnames(chains);
   for i = 1:numel(fields)
      chains.(fields{i}) = chains.(fields{i})(:, keepChains);
   end

   sigmaP = codatable(chains, 'sigmaP', @mean);
   bounds = prctile(chains.sigmaP(:), CI);
   fprintf('Posterior mean of sigma perception is %1.3f, with 95%% CI (%1.3f, %1.3f)\n', sigmaP, bounds);

   sigma3 = codatable(chains, 'sigmaM_1', @mean);
   bounds3 = prctile(chains.sigmaM_1(:), CI);
   fprintf('Posterior mean of sigma memory for set size 3 is %1.3f, with 95%% CI (%1.3f, %1.3f)\n', sigma3, bounds3);
   sigma6 = codatable(chains, 'sigmaM_2', @mean);
   bounds6 = prctile(chains.sigmaM_2(:), CI);
   fprintf('Posterior mean of sigma memory for set size 6 is %1.3f, with 95%% CI (%1.3f, %1.3f)\n', sigma6, bounds6);

   F = figure; clf; hold on;
   setFigure(F, [0.2 0.2 0.4 0.4], '');

   [mu, muBounds] = summarizeMuSnaked(chains, 'mu', dp.nStimuli, CI, ...
      'highStimuliFilterFrom', 67);
   muTruth = dp.stimuli;

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

   for idx = 1:dp.nStimuli
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

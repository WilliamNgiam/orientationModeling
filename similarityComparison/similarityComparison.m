%% Similarity comparison model (Gaussian, mu wrap copies)

clear; close all;
preLoad = true;
printFigures = true;

% graphical model script
modelDir = './';
modelName = 'similarityComparison';
engine = 'jags';

% data sets
dataList = {...
   'tomicBaysSimilarity'; ...
   };

%% constants
pi = 3.1415;
load pantoneColors pantone
fontSize = 18;
CI = [2.5 97.5];

% loop over data
for dataIdx = 1:numel(dataList)
   dataName = dataList{dataIdx};
   switch dataName

      case 'tomicBaysSimilarity'
         dataDir = '../data/';
         dataName = 'tomicBays';
         load([dataDir dataName], 'ds');

         a = ds.aIdx;
         b = ds.bIdx;
         c = ds.cIdx;
         d = ds.dIdx;
         y = ds.response;
         nTrials = ds.nTrials;
         nStimuli = ds.nStimuli;
   end

   %% sampling from graphical model
   % parameters to monitor
   params = {'mu', 'sigma'};

   % MCMC properties
   nChains    = 12;     % number of MCMC chains
   nBurnin    = 2e3;    % number of discarded burn-in samples
   nSamples   = 2e3;    % number of collected samples
   nThin      = 5;      % number of samples between those collected
   doParallel = 1;      % whether MATLAB parallel toolbox parallizes chains

   % assign MATLAB variables to the observed nodes
   data = struct(...
      'a'        , a        , ...
      'b'        , b        , ...
      'c'        , c        , ...
      'd'        , d        , ...
      'y'        , y        , ...
      'nStimuli' , nStimuli , ...
      'nTrials'  , nTrials  );

   % censoring initial values so data have likelihood on first sample
   xAinit = zeros(nTrials, 1);
   xBinit = zeros(nTrials, 1);
   xCinit = zeros(nTrials, 1);
   xDinit = zeros(nTrials, 1);
   for t = 1:nTrials
      [xAinit(t), xBinit(t), xCinit(t), xDinit(t)] = censorSimilarityWrapTrialInits( ...
         ds.a(t), ds.b(t), ds.c(t), ds.d(t), y(t), pi);
   end

   stim = ds.stimuli(:);
   nuInit = ones(nTrials, 1);
   generator = @()similarityWrapChainInits(stim, ...
      xAinit, xBinit, xCinit, xDinit, nuInit, true);

   fileName = sprintf('%s_%s_%s.mat', modelName, dataName, engine);

   if preLoad && isfile(sprintf('storage/%s', fileName))
      fprintf('Loading pre-stored samples for model %s on data %s\n', modelName, dataName);
      load(sprintf('storage/%s', fileName), 'chains', 'stats', 'diagnostics', 'info');
   else
      tic; % start clock
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
         'workingdir'      , sprintf('/tmp/%s', modelName)             , ...
         'verbosity'       , 0                                         , ...
         'saveoutput'      , true                                      , ...
         'parallel'        , doParallel                                );
      fprintf('%s took %f seconds!\n', upper(engine), toc); % show timing

      % convergence of each parameter
      disp('Convergence statistics:')
      grtable(chains, 1.05)

      % basic descriptive statistics
      disp('Descriptive statistics for all chains:')
      codatable(chains);

      fprintf('Saving samples for model %s on data %s\n', modelName, dataName);
      if ~isfolder('storage')
         !mkdir storage
      end
      save(sprintf('storage/%s', fileName), 'chains', 'stats', 'diagnostics', 'info', '-v7.3');

   end

   % just convergent enough chains
   [keepChains, rHat] = findKeepChains(chains.sigma, 2, 1.1);
   fields = fieldnames(chains);
   for i = 1:numel(fields)
      chains.(fields{i}) = chains.(fields{i})(:, keepChains);
   end

   % posterior summary sigma
   sigma = codatable(chains, 'sigma', @mean);
   bounds = prctile(chains.sigma(:), CI);
   fprintf('Posterior mean of sigma is %1.3f, with 95%% CI (%1.3f, %1.3f)\n', sigma, bounds);

   % inferred representation
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
      plot([i i], [0 pi], '-', ...
         'color', pantone.GlacierGray);
      plot([0 pi], [i i], '-', ...
         'color', pantone.GlacierGray);
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
   plot([0 pi], [0 pi], '-', ...
      'color', pantone.AuroraRed, 'linewidth', 0.5);

   % print
   if printFigures
      if ~isfolder('figures')
         !mkdir figures
      end
      print(sprintf('figures/%s_%s.png', dataName, modelName), '-dpng');
      print(sprintf('figures/%s_%s.eps', dataName, modelName), '-depsc');
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

function [xA, xB, xC, xD] = censorSimilarityWrapTrialInits(a, b, c, d, yObs, piVal)
%CENSORSIMILARITYWRAPTRIALINITS Inits for wrap-copy similarity models.
%
% Distance matches similarityComparison_jags.txt:
%   d = min(|x1-x2|, pi-|x1-x2|)
% y=0 => dAB < dCD; y=1 => dAB > dCD.

epsSep = 1e-4;
maxIter = 20;

xA = clampHalfLineSample(a, piVal);
xB = clampHalfLineSample(b, piVal);
xC = clampHalfLineSample(c, piVal);
xD = clampHalfLineSample(d, piVal);

for iter = 1:maxIter
   dAB = halfCircleDist(xA, xB, piVal);
   dCD = halfCircleDist(xC, xD, piVal);
   margin = dAB - dCD;
   ok = (yObs == 0 && margin < -epsSep) || (yObs == 1 && margin > epsSep);
   if ok
      break;
   end
   if yObs == 0
      need = max(dAB - dCD + epsSep, epsSep);
      [xA, xB] = moveHalfLinePairCloser(xA, xB, need, piVal);
   else
      need = max(dCD - dAB + epsSep, epsSep);
      [xC, xD] = moveHalfLinePairCloser(xC, xD, need, piVal);
   end
end

margin = halfCircleDist(xA, xB, piVal) - halfCircleDist(xC, xD, piVal);
ok = (yObs == 0 && margin < 0) || (yObs == 1 && margin > 0);
if ~ok
   error('censorSimilarityWrapTrialInits: trial inits still inconsistent with y.');
end

xA = clampHalfLineSample(xA, piVal);
xB = clampHalfLineSample(xB, piVal);
xC = clampHalfLineSample(c, piVal);
xD = clampHalfLineSample(d, piVal);
end

function x = clampHalfLineSample(x, piVal)
x = mod(x, piVal);
x = min(max(x, 0), piVal - 1e-10);
end

function d = halfCircleDist(u, v, piVal)
d = min(abs(u - v), piVal - abs(u - v));
end

function [u, v] = moveHalfLinePairCloser(u, v, reduction, piVal)
if reduction <= 0
   return;
end
diff = v - u;
if abs(diff) <= piVal / 2
   dir = sign(diff);
   if dir == 0
      dir = 1;
   end
else
   dir = -sign(diff);
   if dir == 0
      dir = 1;
   end
end
half = reduction / 2;
u = clampHalfLineSample(u + dir * half, piVal);
v = clampHalfLineSample(v - dir * half, piVal);
end

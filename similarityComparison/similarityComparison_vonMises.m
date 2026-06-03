% Similarity compatison model

clear; close all;
preLoad = true;
printFigures = true;

% graphical model script
modelDir = './';
modelName = 'similarityComparison_vonMises';
engine = 'jags';
setenv('JAGS_LIBS', fullfile(fileparts(mfilename('fullpath')), '..', 'jags-vonMises'));

% data sets
dataList = {...
   'tomicBaysSimilarity'; ...
   };

%% constants (must match similarityComparison_vonMises_symmetry_jags.txt)
pi = 3.141592653589793;
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
   nChains    = 8;     % number of MCMC chains
   nBurnin    = 0;   % number of discarded burn-in samples
   nSamples   = 50;   % number of collected samples
   nThin      = 1;    % number of samples between those collected
   doParallel = 1;     % whether MATLAB parallel toolbox parallizes chains

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
      [xAinit(t), xBinit(t), xCinit(t), xDinit(t)] = censorSimilarityTrialInits( ...
         ds.a(t), ds.b(t), ds.c(t), ds.d(t), y(t), pi);
   end

   % generator for initialization
   % (note intialization of mu, which encourages the
   % better log-likelihood representation)
   muInit = ds.stimuli(:);
   generator = @()struct(...
      'xAtmp', xAinit, ...
      'xBtmp', xBinit, ...
      'xCtmp', xCinit, ...
      'xDtmp', xDinit, ...
      'mu', muInit);

   fileName = sprintf('%s_%s_%s.mat', modelName, dataName, engine);
   if preLoad && isfile(sprintf('storage/%s', fileName))
      fprintf('Loading pre-stored samples for model %s on data %s\n', modelName, dataName);
      load(sprintf('storage/%s', fileName), 'chains', 'stats', 'diagnostics', 'info');
   else
      tic; % start clock
      [stats, chains, diagnostics, info] = callbayes(engine, ...
         'model'           , sprintf('%s/%s_%s.txt', modelDir, modelName, engine)   , ...   , ...
         'data'            , data                                      , ...
         'outputname'      , 'samples'                                 , ...
         'init'            , generator                                 , ...
         'datafilename'    , modelName                                 , ...
         'initfilename'    , modelName                                 , ...
         'scriptfilename'  , modelName                                 , ...
         'logfilename'     , sprintf('/tmp/%s', modelName)              , ...
         'nchains'         , nChains                                   , ...
         'nburnin'         , nBurnin                                   , ...
         'nsamples'        , nSamples                                  , ...
         'monitorparams'   , params                                    , ...
         'thin'            , nThin                                     , ...
         'workingdir'      , sprintf('/tmp/%s', modelName)              , ...
         'verbosity'       , 0                                         , ...
         'saveoutput'      , true                                      , ...
         'modules'         , {'vonmises'}                                , ...
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

   % log likelihood
   if contains(params, 'yp')
      LL = sum(y.*log(yp)) + sum((1-y).*log(1-yp));
      fprintf('Log-likelihood = %1.4f\n', LL);
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

   mu = codatable(chains, 'mu', @mean);
   muBounds = nan(ds.nStimuli, 2);
   for idx = 1:ds.nStimuli
      muBounds(idx, :) = prctile(chains.(sprintf('mu_%d', idx))(:), CI);
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
      'layer'      , 'top'     , ...
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
      plot(muTruth(idx)*ones(1, 2), muBounds(idx, :),  '-', ...
         'color', pantone.ClassicBlue, ...
         'linewidth', 1);
      plot(muTruth(idx), mu(idx),  'o', ...
         'markerfacecolor', pantone.ClassicBlue, ...
         'markeredgecolor', 'w', ...
         'linewidth', 0.5, ...
         'markersize', 4);

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

function [xA, xB, xC, xD] = censorSimilarityTrialInits(a, b, c, d, yObs, piVal)
%CENSORSIMILARITYTRIALINITS Minimal xAtmp inits consistent with JAGS dinterval(y).
%
% xAtmp are doubled-circle samples; distance matches
% similarityComparison_vonMises_symmetry_jags.txt:
%   d = min(|phi1-phi2|, 2*pi-|phi1-phi2|) / 2
% y=0 => dAB < dCD; y=1 => dAB > dCD.

pi2 = 2 * piVal;
epsSep = 1e-4;
maxIter = 20;

xA = angleToDoubledSample(a, pi2);
xB = angleToDoubledSample(b, pi2);
xC = angleToDoubledSample(c, pi2);
xD = angleToDoubledSample(d, pi2);

for iter = 1:maxIter
   dAB = doubledHalfCircleDist(xA, xB, pi2);
   dCD = doubledHalfCircleDist(xC, xD, pi2);
   margin = dAB - dCD;
   ok = (yObs == 0 && margin < -epsSep) || (yObs == 1 && margin > epsSep);
   if ok
      break;
   end
   if yObs == 0
      need = max(dAB - dCD + epsSep, epsSep);
      [xA, xB] = moveDoubledPairCloser(xA, xB, need, pi2);
   else
      need = max(dCD - dAB + epsSep, epsSep);
      [xC, xD] = moveDoubledPairCloser(xC, xD, need, pi2);
   end
end

margin = doubledHalfCircleDist(xA, xB, pi2) - doubledHalfCircleDist(xC, xD, pi2);
ok = (yObs == 0 && margin < 0) || (yObs == 1 && margin > 0);
if ~ok
   error('censorSimilarityTrialInits: trial inits still inconsistent with y.');
end

xA = clampDoubled(xA, pi2);
xB = clampDoubled(xB, pi2);
xC = clampDoubled(xC, pi2);
xD = clampDoubled(xD, pi2);
end

function phi = angleToDoubledSample(angle, pi2)
% Map half-line physical angle to doubled-circle sample in [0, 2*pi).
phi = mod(2 * angle, pi2);
phi = min(phi, pi2 - 1e-10);
end

function phi = clampDoubled(phi, pi2)
phi = mod(phi, pi2);
phi = min(phi, pi2 - 1e-10);
end

function d = doubledHalfCircleDist(u, v, pi2)
% Match JAGS: min(|phi1-phi2|, 2*pi-|phi1-phi2|) / 2.
d = min(abs(u - v), pi2 - abs(u - v)) / 2;
end

function [u, v] = moveDoubledPairCloser(u, v, halfReduction, pi2)
% Reduce half-circle distance by halfReduction on the doubled circle.
if halfReduction <= 0
   return;
end
circReduction = 2 * halfReduction;
diff = v - u;
if abs(diff) <= pi2 / 2
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
half = circReduction / 2;
u = clampDoubled(u + dir * half, pi2);
v = clampDoubled(v - dir * half, pi2);
end
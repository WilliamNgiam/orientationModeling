function jobs = weekendJobDefs(rootDir, muSRef)
%WEEKENDJOBDEFS Job structs for runWeekendV2Mcmc (priority order).
% Order: clusterEntrop → common → cluster. commonEntrop deferred.

jobs = {};

% --- cluster entropification (often already PASSed; reuse keep) ---
j = baseJob('clusterRepresentation', 'clusterRepresentationEntropification', ...
   'clusterRepresentationEntropificationWeekend', rootDir);
j.kind = 'cluster';
j.params = {'gammaPM', 'gammaPS', 'gammaPrior', 'muP', 'muM', 'muS', 'w'};
j.alignPrefixes = {'muP', 'muM', 'muS'};
j.rHatFields = {'gammaPM', 'gammaPS', 'gammaPrior', 'w'};
j.referenceMu = muSRef;
j.makeInits = @initsClusterEntrop;
jobs{end+1} = j; %#ok<*AGROW>

% --- common ---
j = baseJob('commonRepresentation', 'commonRepresentation', ...
   'commonRepresentationWeekend', rootDir);
j.kind = 'common';
j.params = {'mu', 'sigmaP', 'sigmaM', 'sigmaS'};
j.alignPrefixes = {'mu'};
j.rHatFields = {'sigmaP', 'sigmaM_1', 'sigmaM_2', 'sigmaS'};
j.referenceMu = [];
j.makeInits = @initsCommon;
jobs{end+1} = j;

% --- cluster (after common; commonEntrop waits) ---
j = baseJob('clusterRepresentation', 'clusterRepresentation', ...
   'clusterRepresentationWeekend', rootDir);
j.kind = 'cluster';
j.params = {'gammaPM', 'gammaPS', 'gammaPrior', 'muP', 'muM', 'muS'};
j.alignPrefixes = {'muP', 'muM', 'muS'};
j.rHatFields = {'gammaPM', 'gammaPS', 'gammaPrior'};
j.referenceMu = muSRef;
j.makeInits = @initsCluster;
jobs{end+1} = j;
end

function j = baseJob(folder, jagsBase, tag, rootDir)
j = struct();
j.name = tag;
j.folder = folder;
j.tag = tag;
j.modelFile = fullfile(rootDir, folder, [jagsBase '_jags.txt']);
% Paper product names omit "Weekend" (tag still used for job logging).
mainTag = regexprep(tag, 'Weekend$', '');
j.outMat = fullfile(rootDir, folder, 'storage', ...
   sprintf('%s_tomicBays2024_jags.mat', mainTag));
j.nBurnin = 3000;
% Match SelectMCMC wave length (1500) so seeded chains append cleanly
j.nSamples = 1500;
j.nThin = 3;
j.batchSize = 16;
j.maxChains = 96;
j.minKeep = 4;
j.maxKeep = 8;
j.rHatMax = 1.1;
j.wave = 0;
j.state = 'pending';
j.nChainsTotal = 0;
j.keepChains = [];
j.note = '';
end

function init = initsCluster(dp)
piVal = dp.piVal;
stimuli = dp.stimuli;
nStimuli = numel(stimuli);
muP = nan(3 * nStimuli, 1);
muP(2:nStimuli) = stimuli(2:end) + 0.05 * randn(nStimuli - 1, 1);
muP(2:nStimuli) = min(max(muP(2:nStimuli), 1e-4), piVal - 1e-4);
muP(2) = min(max(muP(2), 1e-3), piVal/2 - 1e-3);
delta = min(max(0.1 * randn(nStimuli, 2), -piVal/2 + 1e-4), piVal/2 - 1e-4);
init = struct( ...
   'muP', muP, ...
   'sigmaP', 0.1 + 0.2 * rand, ...
   'sigmaM', 0.2 + 0.7 * rand(1, 2), ...
   'sigmaS', 0.15 + 0.25 * rand, ...
   'gammaPM', 0.3 + 0.7 * rand, ...
   'gammaPS', 0.3 + 0.7 * rand, ...
   'gammaPrior', 0.3 + 0.7 * rand, ...
   'delta', delta, ...
   'omega3', jitterSimplex([0.85 0.075 0.075]), ...
   'omega6', jitterSimplex([0.55 0.09 0.09 0.09 0.09 0.09]), ...
   'nu', ones(dp.nPTrials, 1), ...
   'xi', ones(dp.nMTrials, 2), ...
   'xA', dp.xAinit, 'xB', dp.xBinit, 'xC', dp.xCinit, 'xD', dp.xDinit);
end

function init = initsClusterEntrop(dp)
init = initsCluster(dp);
init.w = 0.5 + 3.5 * rand;
init.nuS = ones(dp.nSTrials, 1);
init.xA = dp.xAinitE;
init.xB = dp.xBinitE;
init.xC = dp.xCinitE;
init.xD = dp.xDinitE;
end

function init = initsCommon(dp)
piVal = dp.piVal;
stimuli = dp.stimuli;
nStimuli = numel(stimuli);
mu = nan(3 * nStimuli, 1);
mu(2:nStimuli) = stimuli(2:end) + 0.05 * randn(nStimuli - 1, 1);
mu(2:nStimuli) = min(max(mu(2:nStimuli), 1e-4), piVal - 1e-4);
mu(2) = min(max(mu(2), 1e-3), piVal/2 - 1e-3);
init = struct( ...
   'mu', mu, ...
   'sigmaP', 0.1 + 0.2 * rand, ...
   'sigmaM', 0.2 + 0.7 * rand(1, 2), ...
   'sigmaS', 0.15 + 0.25 * rand, ...
   'omega3', jitterSimplex([0.85 0.075 0.075]), ...
   'omega6', jitterSimplex([0.55 0.09 0.09 0.09 0.09 0.09]), ...
   'nu', ones(dp.nPTrials, 1), ...
   'xi', ones(dp.nMTrials, 2), ...
   'xA', dp.xAinit, 'xB', dp.xBinit, 'xC', dp.xCinit, 'xD', dp.xDinit);
end

function init = initsCommonEntrop(dp)
init = initsCommon(dp);
init.w = 0.5 + 3.5 * rand;
init.nuS = ones(dp.nSTrials, 1);
init.xA = dp.xAinitE;
init.xB = dp.xBinitE;
init.xC = dp.xCinitE;
init.xD = dp.xDinitE;
end

function p = jitterSimplex(p)
p = p(:)' .* exp(0.15 * randn(size(p)));
p = max(p, 1e-6);
p = p / sum(p);
end

function [muRef, info] = makeClusterMuSReference(rootDir, piVal)
%MAKECLUSTERMUSREFERENCE Paper-like bent muS from v1 tomicBays storage.
%
% Uses keepChains [11 15 19] when those columns exist (paper script). Saves
% supportingFiles/clusterMuSReference_tomicBays.mat for reuse.

if nargin < 2
   piVal = 3.1415;
end
refPath = fullfile(rootDir, 'supportingFiles', 'clusterMuSReference_tomicBays.mat');
v1Path = fullfile(rootDir, 'clusterRepresentation', 'storage', ...
   'clusterRepresentation_tomicBays_jags.mat');
dataPath = fullfile(rootDir, 'data', 'tomicBays.mat');

info = struct('source', v1Path, 'keepChains', [], 'nStimuli', 0);

if isfile(refPath)
   S = load(refPath, 'muRef', 'info');
   if isfield(S, 'muRef') && ~isempty(S.muRef)
      muRef = S.muRef(:);
      if isfield(S, 'info'), info = S.info; end
      return
   end
end

if ~isfile(v1Path)
   error('makeClusterMuSReference:missing', 'Missing v1 storage %s', v1Path);
end
D = load(dataPath, 'dp');
nStimuli = D.dp.nStimuli;
info.nStimuli = nStimuli;

S = load(v1Path, 'chains');
chains = S.chains;
nCh = size(chains.gammaPM, 2);
paperKeep = [11 15 19];
if max(paperKeep) <= nCh
   keep = paperKeep;
else
   keep = 1:nCh;
end
info.keepChains = keep;

muRef = nan(nStimuli, 1);
for j = 1:nStimuli
   fname = sprintf('muS_%d', j);
   if ~isfield(chains, fname)
      error('makeClusterMuSReference:muS', 'Missing %s in v1 chains', fname);
   end
   muRef(j) = circularMeanHalf(chains.(fname)(:, keep), piVal);
end

save(refPath, 'muRef', 'info');
end

function mu0 = circularMeanHalf(vals, piVal)
phi = 2 * vals(:);
phi0 = atan2(mean(sin(phi)), mean(cos(phi)));
mu0 = mod(phi0 / 2, piVal);
end

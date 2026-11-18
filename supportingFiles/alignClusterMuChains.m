function [chains, nFlips] = alignClusterMuChains(chains, stimuli, piVal, prefixes)
%ALIGNCLUSTERMUCHAINS Per-chain wrap shift of location params toward stimuli.
%
% prefixes: cell of stems, e.g. {'muP','muM','muS'}. First stem is the
% reference for choosing the shift; companions get the same shift.

if nargin < 4 || isempty(prefixes)
   prefixes = {'muP', 'muM', 'muS'};
end
prefixes = prefixes(:)';
muPrefix = prefixes{1};
companions = prefixes(2:end);
nStimuli = numel(stimuli);
sizeField = pickSizeField(chains, muPrefix);
nFlips = 0;
nCh = size(chains.(sizeField), 2);
for j = 2:nStimuli
   fname = sprintf('%s_%d', muPrefix, j);
   if ~isfield(chains, fname)
      continue
   end
   X = chains.(fname);
   sLoc = stimuli(j);
   for c = 1:nCh
      m = circularMeanHalf(X(:, c), piVal);
      cands = [m, m - piVal, m + piVal];
      [~, k] = min(abs(cands - sLoc));
      if k == 1
         continue
      end
      nFlips = nFlips + 1;
      shift = cands(k) - m;
      X(:, c) = X(:, c) + shift;
      for p = 1:numel(companions)
         cname = sprintf('%s_%d', companions{p}, j);
         if isfield(chains, cname)
            chains.(cname)(:, c) = chains.(cname)(:, c) + shift;
         end
      end
   end
   chains.(fname) = X;
end
end

function mu0 = circularMeanHalf(vals, piVal)
phi = 2 * vals;
phi0 = atan2(mean(sin(phi)), mean(cos(phi)));
mu0 = mod(phi0 / 2, piVal);
end

function sizeField = pickSizeField(chains, muPrefix)
cands = {'gammaPM', 'sigmaP', 'sigmaS', 'sigmaM_1', 'w'};
for i = 1:numel(cands)
   if isfield(chains, cands{i})
      sizeField = cands{i};
      return
   end
end
fname = sprintf('%s_2', muPrefix);
if isfield(chains, fname)
   sizeField = fname;
   return
end
fn = fieldnames(chains);
if isempty(fn)
   error('alignClusterMuChains:sizeField', 'chains struct is empty');
end
sizeField = fn{1};
end

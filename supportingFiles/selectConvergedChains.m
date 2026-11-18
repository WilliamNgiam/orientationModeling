function [keepChains, scores, info] = selectConvergedChains(chains, stimuli, piVal, varargin)
%SELECTCONVERGEDCHAINS For common-style models: mu-to-physical band + Rhat.
%
% Optional: 'muPrefix' ('mu'), 'minKeep' (4), 'maxKeep' (8), 'rHatMax' (1.1),
% 'rHatFields' ({'sigmaP','sigmaM_1','sigmaM_2','sigmaS'}), 'slack', 'absSlack'.

p = inputParser;
addParameter(p, 'muPrefix', 'mu', @(s)ischar(s) || isstring(s));
addParameter(p, 'minKeep', 4, @(x)isnumeric(x) && x >= 2);
addParameter(p, 'maxKeep', 8, @(x)isnumeric(x) && x >= 2);
addParameter(p, 'slack', 0.35, @isnumeric);
addParameter(p, 'absSlack', 0.15, @isnumeric);
addParameter(p, 'rHatMax', 1.1, @isnumeric);
addParameter(p, 'rHatFields', {'sigmaP', 'sigmaM_1', 'sigmaM_2', 'sigmaS'}, @iscell);
parse(p, varargin{:});
opt = p.Results;
opt.muPrefix = char(opt.muPrefix);
opt.maxKeep = max(opt.maxKeep, opt.minKeep);
stimuli = stimuli(:);
nStimuli = numel(stimuli);
nCh = size(chains.(opt.rHatFields{1}), 2);
if isfield(chains, 'gammaPM')
   nCh = size(chains.gammaPM, 2);
elseif isfield(chains, 'sigmaS')
   nCh = size(chains.sigmaS, 2);
elseif isfield(chains, 'sigmaP')
   nCh = size(chains.sigmaP, 2);
end

scores = zeros(1, nCh);
for c = 1:nCh
   acc = 0;
   for j = 1:nStimuli
      fname = sprintf('%s_%d', opt.muPrefix, j);
      if ~isfield(chains, fname)
         acc = inf;
         break
      end
      m = circularMeanHalf(chains.(fname)(:, c), piVal);
      acc = acc + halfCircleDist(m, stimuli(j), piVal)^2;
   end
   scores(c) = acc;
end

[sortedScores, order] = sort(scores);
thresh = sortedScores(1) + max(opt.absSlack, opt.slack * max(sortedScores(1), eps));
band = order(sortedScores <= thresh);
band = band(1:min(numel(band), opt.maxKeep));

info = struct('scores', scores, 'band', band, 'threshold', thresh, ...
   'order', order, 'sortedScores', sortedScores, 'keptRhat', nan, ...
   'message', '', 'rHatFields', {opt.rHatFields});

fprintf('selectConvergedChains: best mu SSE=%.4f; %d in band %s\n', ...
   sortedScores(1), numel(band), mat2str(band));

if numel(band) < opt.minKeep
   keepChains = [];
   info.message = sprintf('Only %d chain(s) in mu band (need %d).', ...
      numel(band), opt.minKeep);
   fprintf('%s\n', info.message);
   return
end

keepChains = [];
bestSize = 0;
bestWorstR = inf;
rHatFields = opt.rHatFields(:);
for n = min(opt.maxKeep, numel(band)):-1:opt.minKeep
   combos = nchoosek(band, n);
   for i = 1:size(combos, 1)
      idx = combos(i, :);
      worst = 0;
      ok = true;
      for f = 1:numel(rHatFields)
         name = rHatFields{f};
         if ~isfield(chains, name)
            continue
         end
         r = gelmanrubin(chains.(name)(:, idx), 0, 1, 'rhat');
         if ~isfinite(r) || r > opt.rHatMax
            ok = false;
            break
         end
         worst = max(worst, r);
      end
      if ~ok, continue; end
      if n > bestSize || (n == bestSize && worst < bestWorstR)
         keepChains = sort(idx);
         bestSize = n;
         bestWorstR = worst;
      end
   end
   if bestSize > 0
      break
   end
end
info.keptRhat = bestWorstR;
if isempty(keepChains)
   info.message = sprintf('No subset of %d+ band chains with Rhat <= %.2f.', ...
      opt.minKeep, opt.rHatMax);
   fprintf('%s\n', info.message);
   return
end
fprintf('selectConvergedChains: keep %s (worst Rhat = %.3f)\n', ...
   mat2str(keepChains), bestWorstR);
end

function mu0 = circularMeanHalf(vals, piVal)
phi = 2 * vals;
phi0 = atan2(mean(sin(phi)), mean(cos(phi)));
mu0 = mod(phi0 / 2, piVal);
end

function d = halfCircleDist(a, b, piVal)
d = abs(a - b);
d = min(d, piVal - d);
end

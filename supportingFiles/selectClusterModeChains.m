function [keepChains, scores, info] = selectClusterModeChains(chains, stimuli, piVal, varargin)
%SELECTCLUSTERMODECHAINS Keep chains in one paper-like similarity mode.
%
% Scores each chain's circular-mean muS against a reference curve (paper
% bent representation) when 'referenceMu' is supplied; otherwise against
% physical stimuli (legacy).
%
% Then: tight score band (no forced mixed modes) + largest subset with
% Rhat <= rHatMax on gamma fields (and optional extras), needing minKeep.
%
% Optional name-value:
%   'muPrefix'    ('muS')
%   'referenceMu' ([])     nStimuli x 1 paper-like target
%   'minKeep'     (4)
%   'maxKeep'     (8)
%   'slack'       (0.35)
%   'absSlack'    (0.15)
%   'rHatMax'     (1.1)
%   'rHatFields'  ({'gammaPM','gammaPS','gammaPrior'})
%   'minBendCorr' (0.3)   corr((muS-s),(ref-s)); ignored if no referenceMu

p = inputParser;
addParameter(p, 'muPrefix', 'muS', @(s)ischar(s) || isstring(s));
addParameter(p, 'referenceMu', [], @isnumeric);
addParameter(p, 'minKeep', 4, @(x)isnumeric(x) && x >= 2);
addParameter(p, 'maxKeep', 8, @(x)isnumeric(x) && x >= 2);
addParameter(p, 'slack', 0.35, @isnumeric);
addParameter(p, 'absSlack', 0.15, @isnumeric);
addParameter(p, 'rHatMax', 1.1, @isnumeric);
addParameter(p, 'rHatFields', {'gammaPM', 'gammaPS', 'gammaPrior'}, @iscell);
addParameter(p, 'minBendCorr', 0.3, @isnumeric);
parse(p, varargin{:});
opt = p.Results;
opt.muPrefix = char(opt.muPrefix);
opt.maxKeep = max(opt.maxKeep, opt.minKeep);
stimuli = stimuli(:);
useRef = ~isempty(opt.referenceMu);
if useRef
   opt.referenceMu = opt.referenceMu(:);
end

nStimuli = numel(stimuli);
nCh = size(chains.gammaPM, 2);
scores = zeros(1, nCh);
bendCorr = nan(1, nCh);
for c = 1:nCh
   muHat = nan(nStimuli, 1);
   ok = true;
   for j = 1:nStimuli
      fname = sprintf('%s_%d', opt.muPrefix, j);
      if ~isfield(chains, fname)
         ok = false;
         break
      end
      muHat(j) = circularMeanHalf(chains.(fname)(:, c), piVal);
   end
   if ~ok
      scores(c) = inf;
      continue
   end
   if useRef
      n = min(nStimuli, numel(opt.referenceMu));
      acc = 0;
      for j = 1:n
         acc = acc + halfCircleDist(muHat(j), opt.referenceMu(j), piVal)^2;
      end
      scores(c) = acc;
      dHat = signedHalfDiff(muHat(1:n), stimuli(1:n), piVal);
      dRef = signedHalfDiff(opt.referenceMu(1:n), stimuli(1:n), piVal);
      if std(dHat) > 1e-8 && std(dRef) > 1e-8
         bendCorr(c) = corr(dHat, dRef);
      end
   else
      acc = 0;
      for j = 1:nStimuli
         acc = acc + halfCircleDist(muHat(j), stimuli(j), piVal)^2;
      end
      scores(c) = acc;
   end
end

% Optional bend filter before banding
eligible = 1:nCh;
if useRef && isfinite(opt.minBendCorr)
   okBend = bendCorr >= opt.minBendCorr;
   okBend(~isfinite(bendCorr)) = false;
   if any(okBend)
      eligible = find(okBend);
   end
end

[sortedScores, orderLocal] = sort(scores(eligible));
order = eligible(orderLocal);
thresh = sortedScores(1) + max(opt.absSlack, opt.slack * max(sortedScores(1), eps));
band = order(sortedScores <= thresh);
band = band(1:min(numel(band), opt.maxKeep));

info = struct();
info.scores = scores;
info.bendCorr = bendCorr;
info.sortedScores = sortedScores;
info.order = order;
info.threshold = thresh;
info.band = band;
info.muPrefix = opt.muPrefix;
info.useReference = useRef;
info.minKeep = opt.minKeep;
info.maxKeep = opt.maxKeep;
info.rHatMax = opt.rHatMax;
info.rHatFields = opt.rHatFields;
info.keptRhat = nan;
info.bandRhat = nan;
info.message = '';

fprintf(['selectClusterModeChains: best SSE = %.4f (%s); thresh = %.4f; ', ...
   '%d in band: %s\n'], sortedScores(1), tern(useRef, 'vs paper muS', 'vs physical'), ...
   thresh, numel(band), mat2str(band));

if numel(band) < opt.minKeep
   keepChains = [];
   info.message = sprintf(['Only %d chain(s) in score band (need %d).'], ...
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
      if ~ok
         continue
      end
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
   worstBand = 0;
   for f = 1:numel(rHatFields)
      name = rHatFields{f};
      if ~isfield(chains, name), continue; end
      r = gelmanrubin(chains.(name)(:, band), 0, 1, 'rhat');
      worstBand = max(worstBand, r);
   end
   info.bandRhat = worstBand;
   info.message = sprintf(['No subset of %d+ band chains has Rhat <= %.2f ', ...
      '(band Rhat ~ %.3f).'], opt.minKeep, opt.rHatMax, worstBand);
   fprintf('%s\n', info.message);
   return
end

fprintf('selectClusterModeChains: keep %s (worst Rhat = %.3f; SSE %s; bendCorr %s)\n', ...
   mat2str(keepChains), bestWorstR, mat2str(scores(keepChains), 3), ...
   mat2str(bendCorr(keepChains), 2));
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

function d = signedHalfDiff(a, b, piVal)
% Smallest signed difference a-b on the half-circle, in (-pi/2, pi/2].
d = mod(a - b + piVal/2, piVal) - piVal/2;
end

function s = tern(c, a, b)
if c, s = a; else, s = b; end
end

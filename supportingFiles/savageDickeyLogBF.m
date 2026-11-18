function [logBF, details] = savageDickeyLogBF(posterior, priorDraws, critical)
%SAVAGEDICKEYLOGBF log BF_01 = log f(θ0|data) - log f(θ0) via MV KDE.
%
% BF_01 > 1 favors the point null; BF_01 < 1 favors the alternative.
% Densities are floored so results stay finite; details.censored is true
% when a floor was binding (interpret as a bound, e.g. log BF_01 < -100).

posterior = double(posterior);
priorDraws = double(priorDraws);
if isvector(posterior)
   posterior = posterior(:);
end
if isvector(priorDraws)
   priorDraws = priorDraws(:);
end
critical = double(critical(:)');

[fPosterior, bwPost] = kdeAt(posterior, critical);
[fPrior, bwPrior] = kdeAt(priorDraws, critical);

% Floor relative to prior scale so we never hit log(0)
floorPost = max(fPrior * exp(-100), realmin('double'));
floorPrior = realmin('double');
censoredPost = fPosterior < floorPost;
censoredPrior = fPrior < floorPrior;
fPostUse = max(fPosterior, floorPost);
fPriorUse = max(fPrior, floorPrior);

logBF = log(fPostUse) - log(fPriorUse);
% Keep in a printable range for exp()
logBF = max(min(logBF, 100), -100);

details = struct();
details.fPosterior = fPosterior;
details.fPrior = fPrior;
details.fPosteriorUsed = fPostUse;
details.fPriorUsed = fPriorUse;
details.bwPosterior = bwPost;
details.bwPrior = bwPrior;
details.censored = censoredPost || censoredPrior;
details.censoredPosterior = censoredPost;
details.censoredPrior = censoredPrior;
details.critical = critical;
end

function [f, bw] = kdeAt(draws, critical)
n = size(draws, 1);
d = size(draws, 2);
sigma = std(draws, 0, 1);
bw = sigma * (4 / (n * (d + 4))) ^ (1 / (d + 4));
bw(bw <= 0 | ~isfinite(bw)) = 1e-3;
% Near-boundary nulls (gamma >= 0): nudge bandwidth up a bit for stability
bw = max(bw, 1e-3);
f = mvksdensity(draws, critical, 'Bandwidth', bw);
if ~isfinite(f) || f < 0
   f = 0;
end
end

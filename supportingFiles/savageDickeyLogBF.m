function [logBF, details] = savageDickeyLogBF(posterior, priorDraws, critical)
%SAVAGEDICKEYLOGBF log BF_01 = log f(θ0|data) - log f(θ0) via MV KDE.
%
% Convention (Savage–Dickey):
%   BF_01 = f_post(θ0) / f_prior(θ0)
%   log BF_01 > 0  => favors the POINT NULL (H0)
%   log BF_01 < 0  => favors the ALTERNATIVE (H1)
%   log BF_10 = -log BF_01
%
% Densities are evaluated in log space (log-sum-exp of Gaussian kernels) so
% results stay finite even when the posterior mass is far from θ0 (avoids
% -Inf from underflow in mvksdensity/ksdensity).

posterior = double(posterior);
priorDraws = double(priorDraws);
if isvector(posterior)
   posterior = posterior(:);
end
if isvector(priorDraws)
   priorDraws = priorDraws(:);
end
critical = double(critical(:)');

[logFPost, bwPost, fPostLin] = logKdeAt(posterior, critical);
[logFPrior, bwPrior, fPriorLin] = logKdeAt(priorDraws, critical);

logBF = logFPost - logFPrior;

details = struct();
details.fPosterior = exp(logFPost);          % may underflow to 0 in double
details.fPrior = exp(logFPrior);
details.fPosteriorLinearKDE = fPostLin;      % classical ksdensity (may be 0)
details.fPriorLinearKDE = fPriorLin;
details.logFPosterior = logFPost;
details.logFPrior = logFPrior;
details.bwPosterior = bwPost;
details.bwPrior = bwPrior;
details.logBF01 = logBF;
details.logBF10 = -logBF;
details.favorsNull = logBF > 0;
details.favorsAlternative = logBF < 0;
details.critical = critical;
details.censored = false;
details.censoredPosterior = false;
details.censoredPrior = false;
end

function [logF, bw, fLin] = logKdeAt(draws, critical)
%LOGKDEAT log density at critical via Gaussian product kernel (Silverman bw).
n = size(draws, 1);
d = size(draws, 2);
sigma = std(draws, 0, 1);
bw = sigma * (4 / (n * (d + 4))) ^ (1 / (d + 4));
bw(bw <= 0 | ~isfinite(bw)) = 1e-3;
bw = max(bw, 1e-3);

% Classical linear KDE (for diagnostics; may underflow to 0)
try
   if d == 1
      fLin = ksdensity(draws, critical, 'Bandwidth', bw);
   else
      fLin = mvksdensity(draws, critical, 'Bandwidth', bw);
   end
   if ~isfinite(fLin) || fLin < 0
      fLin = 0;
   end
catch
   fLin = 0;
end

% Log-sum-exp of kernels: log f(x) = logsumexp_i log K_h(x-x_i) - log(n)
logNorm = -sum(log(bw)) - 0.5 * d * log(2 * pi);
z = (draws - critical) ./ bw;           % n x d
logKi = logNorm - 0.5 * sum(z.^2, 2);  % n x 1
m = max(logKi);
logF = m + log(sum(exp(logKi - m))) - log(n);
if ~isfinite(logF)
   % Pathological: fall back to nearest-sample kernel (still finite)
   [~, j] = min(sum(z.^2, 2));
   logF = logKi(j) - log(n);
end
end

function writeSummariesStandalone
root = '/home/mdlee/Dropbox/GitHub/orientationModeling';
addpath(fullfile(root, 'supportingFiles'));
outFile = fullfile(root, 'analysis', 'results', 'statisticalSummaries.txt');
fid = fopen(outFile, 'w');
assert(fid > 0);
CIbounds = [2.5 97.5];
scale = 4;

% --- Perceptual ---
chains = loadKeep(fullfile(root,'perceptualReproduction/storage/perceptualReproductionMCMC_tomicBays2024_jags.mat'), 'sigma');
fprintf(fid, '\n----------\nPerceptual Reproduction\n\n');
writeSigma(fid, chains.sigma(:), 'sigma');

% --- Memory ---
S = load(fullfile(root,'memoryReproduction/storage/memoryReproductionMCMC_tomicBays2024_jags.mat'), 'chains');
chains = S.chains;
[keep, ~] = findKeepChains(chains.sigma_1, 2, 1.1);
chains = subset(chains, keep);
fprintf(fid, '\n----------\nMemory Reproduction\n\n');
writeSigma(fid, chains.sigma_1(:), 'sigma for set size 3');
writeSigma(fid, chains.sigma_2(:), 'sigma for set size 6');
omega3 = get_matrix_from_coda(chains, 'omega3', @mean);
omega6 = get_matrix_from_coda(chains, 'omega6', @mean);
fprintf(fid, 'Posterior means for omega for set size 3 are (%1.3f, %1.3f, %1.3f)\n', omega3);
fprintf(fid, 'Posterior means for omega for set size 6 are (%1.3f, %1.3f, %1.3f, %1.3f, %1.3f, %1.3f)\n', omega6);
for i = 1:3
   b = prctile(chains.(sprintf('omega3_%d',i))(:), CIbounds);
   fprintf(fid, 'Posterior mean of omega_%d for set size 3 is %1.3f, with 95%% CI (%1.3f, %1.3f)\n', i, omega3(i), b);
end
for i = 1:6
   b = prctile(chains.(sprintf('omega6_%d',i))(:), CIbounds);
   fprintf(fid, 'Posterior mean of omega_%d for set size 6 is %1.3f, with 95%% CI (%1.3f, %1.3f)\n', i, omega6(i), b);
end
[logBF3,~] = savageDickeyLogBF( ...
   [chains.omega3_1(:) chains.omega3_2(:) chains.omega3_3(:)], ...
   [chains.omega3prior_1(:) chains.omega3prior_2(:) chains.omega3prior_3(:)], [1 0 0]);
[logBF6,~] = savageDickeyLogBF( ...
   [chains.omega6_1(:) chains.omega6_2(:) chains.omega6_3(:) chains.omega6_4(:) chains.omega6_5(:) chains.omega6_6(:)], ...
   [chains.omega6prior_1(:) chains.omega6prior_2(:) chains.omega6prior_3(:) chains.omega6prior_4(:) chains.omega6prior_5(:) chains.omega6prior_6(:)], ...
   [1 0 0 0 0 0]);
fprintf(fid, 'Savage-Dickey log BF_01 (positive favors null / no-swap; negative favors alternative / swap).\n');
fprintf(fid, '%s\n', sent('memory set-size 3 no-swap null', logBF3));
fprintf(fid, '%s\n', sent('memory set-size 6 no-swap null', logBF6));

% --- Similarity ---
chains = loadKeep(fullfile(root,'similarityComparison/storage/similarityComparison_tomicBays2024_jags.mat'), 'sigma');
fprintf(fid, '\n----------\nSimilarity Comparison\n\n');
writeSigma(fid, chains.sigma(:), 'sigma');

% --- Similarity entropification ---
chains = loadKeep(fullfile(root,'similarityComparison/storage/similarityComparisonEntropification_tomicBays2024_jags.mat'), 'sigma');
fprintf(fid, '\n----------\nSimilarity Comparison Entropification\n\n');
writeSigma(fid, chains.sigma(:), 'sigma');
if isfield(chains,'w')
   writeSigma(fid, chains.w(:), 'w');
end

% --- Cluster Mode A ---
fprintf(fid, '\n----------\nCluster Bayes Factor (Mode A / DIC-preferred)\n\n');
chains = loadModeA({ ...
   fullfile(root,'clusterRepresentation/storage/clusterRepresentation_tomicBays2024_jags.mat'), ...
   fullfile(root,'clusterRepresentation/storage/clusterRepresentationSelectMCMC_tomicBays2024_jags.mat')});
writeGamma(fid, chains, CIbounds);
writeClusterBF(fid, chains, scale);

% --- Cluster entrop ---
S = load(fullfile(root,'clusterRepresentation/storage/clusterRepresentationEntropification_tomicBays2024_jags.mat'), 'chains');
chains = S.chains;
fprintf(fid, '\n----------\nCluster Entropification (paper keep)\n\n');
writeGamma(fid, chains, CIbounds);
if isfield(chains,'w'), writeSigma(fid, chains.w(:), 'w'); end
writeClusterBF(fid, chains, scale);

fclose(fid);
fprintf('Wrote %s\n', outFile);
end

function chains = loadKeep(matPath, field)
S = load(matPath, 'chains');
chains = S.chains;
[keep, ~] = findKeepChains(chains.(field), 2, 1.1);
chains = subset(chains, keep);
fprintf('  %s keep=%s\n', matPath, mat2str(keep));
end

function chains = subset(chains, keep)
fn = fieldnames(chains);
for i = 1:numel(fn)
   chains.(fn{i}) = chains.(fn{i})(:, keep);
end
end

function writeSigma(fid, x, name)
b = prctile(x, [2.5 97.5]);
fprintf(fid, 'Posterior mean of %s is %1.3f, with 95%% CI (%1.3f, %1.3f)\n', name, mean(x), b);
end

function writeGamma(fid, chains, CI)
writeSigma(fid, chains.gammaPM(:), 'gamma_pm');
writeSigma(fid, chains.gammaPS(:), 'gamma_ps');
end

function writeClusterBF(fid, chains, scale)
fprintf(fid, ['Savage-Dickey log BF_01 (positive favors null / common map; ', ...
   'negative favors alternative / separate maps). ', ...
   'Densities via Gaussian KDE in log space (bandwidth x%d).\n'], scale);
Xpost = [chains.gammaPS(:) chains.gammaPM(:)];
g0 = chains.gammaPrior(:);
Xpri = [g0, g0(randperm(numel(g0)))];
fprintf(fid, '%s\n', sent('joint null (gamma_pm, gamma_ps)=(0,0)', kdeBF(Xpost, Xpri, [0 0], scale)));
fprintf(fid, '%s\n', sent('perceptual-memory null (gamma_pm=0)', kdeBF(chains.gammaPM(:), chains.gammaPrior(:), 0, scale)));
fprintf(fid, '%s\n', sent('perceptual-similarity null (gamma_ps=0)', kdeBF(chains.gammaPS(:), chains.gammaPrior(:), 0, scale)));
end

function logBF = kdeBF(post, prior, crit, scale)
if isvector(post), post = post(:); end
if isvector(prior), prior = prior(:); end
crit = double(crit(:)');
logBF = logKde(post, crit, scale) - logKde(prior, crit, scale);
end

function logF = logKde(draws, critical, scale)
n = size(draws,1); d = size(draws,2);
bw = std(draws,0,1) * (4/(n*(d+4)))^(1/(d+4));
bw(bw<=0|~isfinite(bw)) = 1e-3;
bw = max(bw,1e-3)*scale;
logNorm = -sum(log(bw)) - 0.5*d*log(2*pi);
z = (draws - critical)./bw;
logKi = logNorm - 0.5*sum(z.^2,2);
m = max(logKi);
logF = m + log(sum(exp(logKi-m))) - log(n);
end

function str = sent(label, logBF01)
if logBF01 > 0, fav = 'NULL (common map / no-swap)';
elseif logBF01 < 0, fav = 'ALTERNATIVE (separate maps / swap)';
else, fav = 'neither'; end
str = sprintf('  %s: log BF_01 = %.2f (favors %s); equivalently log BF_10 = %.2f', ...
   label, logBF01, fav, -logBF01);
end

function chains = loadModeA(mats)
modeA_pmMax = 0.8; modeA_psLo = 0.8; modeA_psHi = 2.2;
modeA_center = [0.1, 1.5];
minKeep = 4; maxKeep = 8; rHatMax = 1.1;
rHatFields = {'gammaPM','gammaPS','gammaPrior'};
for i = 1:numel(mats)
   if ~isfile(mats{i}), continue; end
   S = load(mats{i});
   if isfield(S,'chainsAll') && isfield(S.chainsAll,'gammaPM')
      src = S.chainsAll;
   elseif isfield(S,'chains') && isfield(S.chains,'gammaPM')
      src = S.chains;
   else
      continue
   end
   mPM = median(src.gammaPM,1); mPS = median(src.gammaPS,1);
   cand = find(mPM < modeA_pmMax & mPS >= modeA_psLo & mPS <= modeA_psHi);
   fprintf('Mode A cand %s: n=%d\n', mats{i}, numel(cand));
   if numel(cand) < minKeep, continue; end
   dCenter = (mPM(cand) - modeA_center(1)).^2 + (mPS(cand) - modeA_center(2)).^2;
   [~, ord] = sort(dCenter);
   band = cand(ord);
   if numel(band) > 12, band = band(1:12); end
   keep = []; worstR = inf;
   for n = min(maxKeep, numel(band)):-1:minKeep
      combos = nchoosek(band, n);
      bestWorst = inf; bestIdx = [];
      for c = 1:size(combos, 1)
         idx = combos(c, :);
         worst = 0; ok = true;
         for f = 1:numel(rHatFields)
            r = gelmanrubin(src.(rHatFields{f})(:, idx), 0, 1, 'rhat');
            if ~isfinite(r) || r > rHatMax, ok = false; break; end
            worst = max(worst, r);
         end
         if ok && worst < bestWorst
            bestWorst = worst; bestIdx = sort(idx);
         end
      end
      if ~isempty(bestIdx)
         keep = bestIdx; worstR = bestWorst; break;
      end
   end
   if numel(keep) >= minKeep
      chains = subset(src, keep);
      fprintf('Mode A keep=%s worstR=%.3f\n', mat2str(keep), worstR);
      return
   end
end
error('writeSummariesStandalone:noModeA', 'No Mode A keep found');
end

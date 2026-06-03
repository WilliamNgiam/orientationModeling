%% Example trials: compressed vs physical agreement with diverging model probs
%
% Diagonal: mu = physical orientation. Compressed: posterior mu1 (mixture).
% Plots 10 trials where compressed agrees (physical wrong) and 10 where
% physical agrees (compressed wrong), each with |P_compressed(AB) - P_diagonal(AB)| large.

clear; close all;

thisDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(thisDir);
piVal = pi;
pi2 = 2 * piVal;

storageFile = fullfile(rootDir, 'similarityComparison', 'storage', ...
   'similarityComparison_mixture_tomicBays_jags.mat');
if ~isfile(storageFile)
   error(['Missing %s\nRun similarityComparison_mixture in similarityComparison/ ' ...
      'or copy similarityComparison_mixture_tomicBays_jags.mat into that storage/ folder.'], ...
      storageFile);
end
load(storageFile, 'chains');
load(fullfile(rootDir, 'data', 'tomicBays.mat'), 'ds');
load(fullfile(rootDir, 'supportingFiles', ...
   'pantoneColors.mat'), 'pantone');

nTrials = ds.nTrials;
nStimuli = ds.nStimuli;
y = ds.response(:);
aIdx = ds.aIdx(:);
bIdx = ds.bIdx(:);
cIdx = ds.cIdx(:);
dIdx = ds.dIdx(:);

% Ensure mode 1 is the compressed (lower high-stimulus endpoint) curve
m1End = mean(chains.(sprintf('mu1_%d', nStimuli))(:));
m2End = mean(chains.(sprintf('mu2_%d', nStimuli))(:));
if m2End < m1End
   chains = swapMixtureChainsLocal(chains, nStimuli);
end

[mu1Draws, sigmaDraws] = stackCompressedPosterior(chains, nStimuli);
nPost = size(mu1Draws, 1);
nRep = 250;
sigmaMean = mean(sigmaDraws);

mu1Mean = mean(mu1Draws, 1);
fprintf('Computing compressed choice probabilities...\n');
pABcompressed = trialChoiceProbFromMuMap(mu1Mean, aIdx, bIdx, cIdx, dIdx, ...
   nTrials, sigmaMean, nRep, pi2, piVal);

% Diagonal: mu equals physical angle at each bar (identity map)
fprintf('Computing diagonal (physical mu) choice probabilities...\n');
pABdiagonal = trialChoiceProbPhysical(ds.a(:), ds.b(:), ds.c(:), ds.d(:), ...
   nTrials, sigmaMean, nRep, pi2, piVal);

% Refine compressed probs for top candidates
choseAB = y == 0;
physCloserAB = physicalABCloser(ds, piVal);
physicalExplains = physCloserAB == choseAB;

compressedExplains = (pABcompressed >= 0.5) == choseAB;
diagonalExplains = (pABdiagonal >= 0.5) == choseAB;
probDiff = abs(pABcompressed - pABdiagonal);

minProbDiff = 0.25;
pCorrectCompressed = choseAB .* pABcompressed + (~choseAB) .* (1 - pABcompressed);
pCorrectDiagonal = choseAB .* pABdiagonal + (~choseAB) .* (1 - pABdiagonal);

maskCompressed = ~physicalExplains & compressedExplains & probDiff >= minProbDiff;
maskPhysical = physicalExplains & ~compressedExplains & probDiff >= minProbDiff;

scoreCompressed = nan(nTrials, 1);
scoreCompressed(maskCompressed) = probDiff(maskCompressed) + pCorrectCompressed(maskCompressed);

scorePhysical = nan(nTrials, 1);
scorePhysical(maskPhysical) = probDiff(maskPhysical) + pCorrectDiagonal(maskPhysical);

T = table((1:nTrials)', probDiff, scoreCompressed, scorePhysical, ...
   pABcompressed, pABdiagonal, y, physCloserAB, ...
   physicalExplains, compressedExplains, diagonalExplains, ...
   'VariableNames', {'trial', 'probDiff', 'score_compressed', 'score_physical', ...
   'pAB_compressed', 'pAB_diagonal', 'response', 'physCloserAB', ...
   'physicalExplains', 'compressedExplains', 'diagonalExplains'});

outDir = fullfile(thisDir, 'figures');
if ~isfolder(outDir)
   mkdir(outDir);
end
writetable(T, fullfile(outDir, 'similarityRepresentationTrials_all.csv'));

fprintf('Compressed agrees, physical wrong (|DeltaP|>=%.2f): %d trials\n', ...
   minProbDiff, sum(maskCompressed));
fprintf('Physical agrees, compressed wrong (|DeltaP|>=%.2f): %d trials\n', ...
   minProbDiff, sum(maskPhysical));

nPick = 10;
examplesPerFigure = 5;
pAB = [pABcompressed, pABdiagonal];

sets = { ...
   struct('name', 'compressedAgrees', 'mask', maskCompressed, 'scoreVar', 'score_compressed', ...
   'summary', 'Physical wrong; compressed agrees', 'prefix', 'similarityCompressedAgrees_probDiff'), ...
   struct('name', 'physicalAgrees', 'mask', maskPhysical, 'scoreVar', 'score_physical', ...
   'summary', 'Physical agrees; compressed wrong', 'prefix', 'similarityPhysicalAgrees_probDiff')};

for s = 1:numel(sets)
   cfg = sets{s};
   sub = T(cfg.mask, :);
   sub.score = sub.(cfg.scoreVar);
   picked = pickDiverseTrials(sub, nPick);
   if height(picked) < nPick
      warning('%s: only %d qualifying trials (requested %d).', ...
         cfg.name, height(picked), nPick);
   end
   writetable(picked, fullfile(outDir, [cfg.prefix '_picked.csv']));
   fprintf('\n%s:\n', cfg.name);
   disp(picked(:, {'trial', 'score', 'probDiff', 'pAB_compressed', 'pAB_diagonal', 'response'}));

   nFigures = ceil(height(picked) / examplesPerFigure);
   for figIdx = 1:nFigures
      i0 = (figIdx - 1) * examplesPerFigure + 1;
      i1 = min(figIdx * examplesPerFigure, height(picked));
      batch = picked(i0:i1, :);
      figPath = fullfile(outDir, sprintf('%s_%02d.png', cfg.prefix, figIdx));
      plotDiscriminatingSimilarityExamples(batch, ds, pAB, pantone, piVal, figPath, cfg.summary);
      fprintf('Wrote %s\n', figPath);
   end
end

%% --- local functions ---

function physCloserAB = physicalABCloser(ds, piVal)
n = ds.nTrials;
physCloserAB = false(n, 1);
for t = 1:n
   dAB = angDist(ds.a(t), ds.b(t), piVal);
   dCD = angDist(ds.c(t), ds.d(t), piVal);
   physCloserAB(t) = dAB < dCD;
end
end

function d = angDist(a, b, piVal)
a = mod(a, piVal);
b = mod(b, piVal);
d = min(abs(a - b), piVal - abs(a - b));
end

function pAB = trialChoiceProbFromMuMap(mu, aIdx, bIdx, cIdx, dIdx, nTrials, sigma, nRep, pi2, piVal)
pAB = nan(nTrials, 1);
kap = 1 / (sigma^2 + 1e-8);
for t = 1:nTrials
   muA = mu(aIdx(t));
   muB = mu(bIdx(t));
   muC = mu(cIdx(t));
   muD = mu(dIdx(t));
   phiA = sampleDoubledVonMises(muA, kap, nRep, piVal, pi2);
   phiB = sampleDoubledVonMises(muB, kap, nRep, piVal, pi2);
   phiC = sampleDoubledVonMises(muC, kap, nRep, piVal, pi2);
   phiD = sampleDoubledVonMises(muD, kap, nRep, piVal, pi2);
   pAB(t) = mean(similarityHalfCircleDist(phiA, phiB, pi2) < ...
      similarityHalfCircleDist(phiC, phiD, pi2));
end
end

function pAB = trialChoiceProbPhysical(angA, angB, angC, angD, nTrials, sigma, nRep, pi2, piVal)
pAB = nan(nTrials, 1);
kap = 1 / (sigma^2 + 1e-8);
for t = 1:nTrials
   phiA = sampleDoubledVonMises(angA(t), kap, nRep, piVal, pi2);
   phiB = sampleDoubledVonMises(angB(t), kap, nRep, piVal, pi2);
   phiC = sampleDoubledVonMises(angC(t), kap, nRep, piVal, pi2);
   phiD = sampleDoubledVonMises(angD(t), kap, nRep, piVal, pi2);
   pAB(t) = mean(similarityHalfCircleDist(phiA, phiB, pi2) < ...
      similarityHalfCircleDist(phiC, phiD, pi2));
end
end

function chainsOut = swapMixtureChainsLocal(chains, nStimuli)
chainsOut = chains;
for j = 1:nStimuli
   n1 = sprintf('mu1_%d', j);
   n2 = sprintf('mu2_%d', j);
   tmp = chainsOut.(n1);
   chainsOut.(n1) = chainsOut.(n2);
   chainsOut.(n2) = tmp;
end
end

function [mu1, sigma] = stackCompressedPosterior(chains, nStimuli)
sigma = chains.sigma(:);
nDraws = numel(sigma);
mu1 = nan(nDraws, nStimuli);
for j = 1:nStimuli
   mu1(:, j) = chains.(sprintf('mu1_%d', j))(:);
end
end

function phi = sampleDoubledVonMises(mu, kappa, n, piVal, pi2)
meanPhi = mod(2 * mu, pi2);
if kappa < 1e-6
   phi = mod(meanPhi + pi2 * rand(n, 1), pi2);
   return;
end
tau = 1 + sqrt(1 + 4 * kappa^2);
rho = (tau - sqrt(2 * tau)) / (2 * kappa);
r = (1 + rho^2) / (2 * rho);
phi = nan(n, 1);
filled = 0;
while filled < n
   m = n - filled;
   U1 = rand(m, 1);
   U2 = rand(m, 1);
   z = cos(piVal * U1);
   f = (1 + r * z) ./ (r + z);
   c = kappa * (r - f);
   ok = (c .* (2 - c) - U2 > 0) | (log(max(U2, eps)) < c);
   cnt = sum(ok);
   if cnt == 0
      continue;
   end
   zOk = z(ok);
   theta0 = acos(min(max((1 + r * zOk) ./ (r + zOk), -1), 1));
   theta0 = sign(rand(cnt, 1) - 0.5) .* theta0;
   phi(filled + (1:cnt)) = mod(meanPhi + theta0, pi2);
   filled = filled + cnt;
end
end

function d = similarityHalfCircleDist(phi1, phi2, pi2)
d = min(abs(phi1 - phi2), pi2 - abs(phi1 - phi2)) / 2;
end

function picked = pickDiverseTrials(T, nPick)
T = sortrows(T, 'score', 'descend');
picked = T(1, :);
usedT = picked.trial(1);
rem = T(2:end, :);
while height(picked) < nPick && ~isempty(rem)
   best = 1;
   bestScore = -inf;
   for i = 1:height(rem)
      bonus = 0.05 * ~ismember(rem.trial(i), usedT);
      s = rem.score(i) * (1 + bonus);
      if s > bestScore
         bestScore = s;
         best = i;
      end
   end
   picked = [picked; rem(best, :)]; %#ok<AGROW>
   usedT(end+1) = rem.trial(best);
   rem(best, :) = [];
end
end

function plotDiscriminatingSimilarityExamples(picked, ds, pAB, pantone, piVal, savePath, summaryLine)
if nargin < 8
   summaryLine = '';
end
n = height(picked);
nCols = 2;
nRows = n;
fig = figure('Color', 'w', 'Position', [40 40 720 140 + 118 * nRows], 'Units', 'pixels');
tlo = tiledlayout(nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');
tlo.TileIndexing = 'rowmajor';

xLim = [-1.12 1.12];
yLim = [-0.02 1.06];
rayLen = 0.92;
labelFs = 9;

for k = 1:n
   tr = picked.trial(k);
   angA = mod(ds.a(tr), piVal);
   angB = mod(ds.b(tr), piVal);
   angC = mod(ds.c(tr), piVal);
   angD = mod(ds.d(tr), piVal);
   choseAB = ds.response(tr) == 0;
   physCloserAB = angDist(ds.a(tr), ds.b(tr), piVal) < angDist(ds.c(tr), ds.d(tr), piVal);
   pCompressed = pAB(tr, 1);
   pDiagonal = pAB(tr, 2);

   ax1 = nexttile(tlo);
   hold(ax1, 'on');
   setupSemicircleAxes(ax1, xLim, yLim, pantone);
   drawSemicircleRay(ax1, angA, pantone.ClassicBlue, '-', 2.4, rayLen);
   drawSemicircleRay(ax1, angB, pantone.ClassicBlue, '-', 2.4, rayLen);
   drawSemicircleRay(ax1, angC, pantone.AuroraRed, '-', 2.4, rayLen);
   drawSemicircleRay(ax1, angD, pantone.AuroraRed, '-', 2.4, rayLen);
   labelRay(ax1, angA, 'A', rayLen, labelFs);
   labelRay(ax1, angB, 'B', rayLen, labelFs);
   labelRay(ax1, angC, 'C', rayLen, labelFs);
   labelRay(ax1, angD, 'D', rayLen, labelFs);
   title(ax1, sprintf('trial %d', tr), 'FontSize', 10);

   ax2 = nexttile(tlo);
   axis(ax2, 'off');
   if choseAB
      choiceStr = 'Participant: AB closer';
      col = pantone.ClassicBlue;
   else
      choiceStr = 'Participant: CD closer';
      col = pantone.AuroraRed;
   end
   if physCloserAB
      physStr = 'Physical: AB closer';
   else
      physStr = 'Physical: CD closer';
   end
   text(ax2, 0, 0.88, choiceStr, 'Units', 'normalized', ...
      'FontSize', 11, 'Color', col, 'FontWeight', 'bold');
   text(ax2, 0, 0.70, physStr, 'Units', 'normalized', 'FontSize', 10, ...
      'Color', [0.4 0.4 0.4]);
   text(ax2, 0, 0.50, sprintf('Compressed: P(choose AB) = %.2f', pCompressed), ...
      'Units', 'normalized', 'FontSize', 10);
   text(ax2, 0, 0.32, sprintf('Diagonal (physical \\mu): P(choose AB) = %.2f', pDiagonal), ...
      'Units', 'normalized', 'FontSize', 10);
   if ~isempty(summaryLine)
      text(ax2, 0, 0.12, summaryLine, 'Units', 'normalized', ...
         'FontSize', 9, 'Color', [0.35 0.35 0.35]);
   end
end

exportgraphics(fig, savePath, 'Resolution', 150);
end

function labelRay(ax, angle, lbl, rayLen, fontSize)
r = rayLen * 1.08;
x = r * cos(angle);
y = r * sin(angle);
text(ax, x, y, lbl, 'HorizontalAlignment', 'center', ...
   'VerticalAlignment', 'middle', 'FontSize', fontSize, 'FontWeight', 'bold');
end

function setupSemicircleAxes(ax, xLim, yLim, pantone)
set(ax, 'xlim', xLim, 'ylim', yLim, 'xtick', [], 'ytick', [], ...
   'box', 'off', 'clipping', 'off', 'layer', 'top');
axis(ax, 'equal');
axis(ax, 'off');
th = linspace(0, pi, 120);
plot(ax, cos(th), sin(th), '-', 'Color', pantone.GlacierGray, 'LineWidth', 0.9);
plot(ax, xLim, [0 0], '-', 'Color', pantone.GlacierGray, 'LineWidth', 0.9);
plot(ax, 0, 0, 'o', 'markerfacecolor', pantone.GlacierGray, ...
   'markeredgecolor', 'w', 'markersize', 2.5);
end

function drawSemicircleRay(ax, angle, col, ls, lw, len)
angle = mod(angle, pi);
plot(ax, [0 len * cos(angle)], [0 len * sin(angle)], ...
   'Color', col, 'LineStyle', ls, 'LineWidth', lw);
end

% figures and non-modeling analyses for Tomic & Bays (canonical parse).
%
% DATA: tomicBays2024 (parse_tomicBays2024).
%

clear;
close all;

printFigures = true;

analysisList = {...
   % --- paper ---
   % 'fourScatter'; ...
   % 'clusterSavageDickey'; ...      % full [0,3]; DIC Mode A (~0, 1.5) + R̂ keep
   % 'clusterSavageDickeyEntropification'; ...  % zoom [0,1]
   'statisticalSummaries'; ...
   % 'twoRepresentations'; ...
   % 'swapExamples'; ...
   % --- supplementary (SI) ---
   % 'memoryNoSwap'; ...                 % μ scatter, no-swap memory MCMC
   % 'similarityEntropification'; ...    % μ scatter, similarity + entropification
   % 'descriptiveAdequacyPerceptual'; ...
   % 'descriptiveAdequacyMemory'; ...
   % 'descriptiveAdequacyMemoryNoSwap'; ...  % needs plain NoSwap yP mat
   % 'descriptiveAdequacySimilarity'; ...
   };









% load data
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'supportingFiles'));
dataDir = '../data/';
dataName = 'tomicBays2024';
load([dataDir dataName], 'ds', 'dp', 'dm');

% constants
load pantoneColors pantone;
pi = 3.1415;

% loops over analyses
for analysisIdx = 1:numel(analysisList)
   analysisName = analysisList{analysisIdx};

   switch analysisName
      case 'fourScatter'

         % perceptual/memory/similarity/common: tomicBays2024
         % (common keep [1 3 5 8 10 11 12]).
         fileList = {...
            '../perceptualReproduction/storage/perceptualReproductionMCMC_tomicBays2024_jags'; ...
            '../memoryReproduction/storage/memoryReproductionMCMC_tomicBays2024_jags'; ...
            '../similarityComparison/storage/similarityComparison_tomicBays2024_jags'; ...
            '../commonRepresentation/storage/commonRepresentation_tomicBays2024_jags'};

         fontSize = 20;
         CIbounds = [2.5 97.5];
         labels = {'perceptual', 'memory', 'similarity', 'common'};

         nModels = numel(fileList);
         [nRows, nCols] = subplotArrange(nModels);

         F = figure; clf; hold on;
         setFigure(F, [0.2 0.2 0.6 0.7], '');

         for modelIdx = 1:nModels
            fileName = fileList{modelIdx};
            fprintf('Loading pre-stored samples from file %s\n', fileName);
            % Load chains only (avoids stale function-handle warnings from inits).
            load(sprintf('%s', fileName), 'chains');

            % just keep converged chains (weekend common is already subset)
            switch modelIdx
               case 1 % perceptual
                  [keepChains, rHat] = findKeepChains(chains.sigma, 2, 1.1);
               case 2 % memory
                  [keepChains, rHat] = findKeepChains(chains.sigma_1, 2, 1.1);
               case 3 % similarity
                  [keepChains, rHat] = findKeepChains(chains.sigma, 2, 1.1);
               case 4 % common
                  [keepChains, rHat] = findKeepChains(chains.sigmaS, 2, 1.1);
            end
            fields = fieldnames(chains);
            for i = 1:numel(fields)
               chains.(fields{i}) = chains.(fields{i})(:, keepChains);
            end

            if modelIdx ~=4
               mu = nan(dm.nStimuli, 1);
               muBounds = nan(dm.nStimuli, 2);
               for idx = 1:dm.nStimuli
                  vals = chains.(sprintf('mu_%d', idx))(:);
                  [mu(idx), muBounds(idx, :)] = summarizeHalfCircle(vals, CIbounds);
               end
            else
               mu = nan(dp.nStimuli, 1);
               muBounds = nan(dp.nStimuli, 2);
               for idx = 1:dp.nStimuli
                  vals = chains.(sprintf('mu_%d', idx))(:);
                  if idx > 66 % so hard to get convergent chains here
                     vals = vals(vals > 2);
                  end
                  [mu(idx), muBounds(idx, :)] = summarizeHalfCircle(vals, CIbounds);
               end
            end
            muTruth = dp.stimuli;

            subplot(nRows, nCols, modelIdx); cla; hold on;
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

            text(0, pi, labels{modelIdx}, ...
               'fontsize', fontSize-4, ...
               'fontweight', 'normal', ...
               'vert', 'bot', 'hor', 'lef');
            moveAxis(gca, [1 1 0.95 0.95], [0 0.025 0 0]);
            Raxes(gca, 0.02, 0.01);

            for i = pi/4:pi/4:3*pi/4
               plot([i i], [0 pi], '-', ...
                  'color', pantone.GlacierGray);
               plot([0 pi], [i i], '-', ...
                  'color', pantone.GlacierGray);
            end

            for idx = 1:dp.nStimuli
               if muBounds(idx, 1) > muBounds(idx, 2)
                  plot(muTruth(idx) * [1 1], [0 muBounds(idx, 2)], '-', ...
                     'color', pantone.ClassicBlue, 'linewidth', 1);
                  plot(muTruth(idx) * [1 1], [muBounds(idx, 1) pi], '-', ...
                     'color', pantone.ClassicBlue, 'linewidth', 1);
               else
                  plot(muTruth(idx) * [1 1], muBounds(idx, :), '-', ...
                     'color', pantone.ClassicBlue, 'linewidth', 1);
               end
               plot(muTruth(idx), mu(idx), 'o', ...
                  'markerfacecolor', pantone.ClassicBlue, ...
                  'markeredgecolor', 'w', 'linewidth', 0.5, 'markersize', 4);
            end
            plot([0 pi], [0 pi], '-', ...
               'color', pantone.AuroraRed, 'linewidth', 0.5);
         end

         [~, T(1)] = suplabel('Physical','x');
         [~, T(2)] = suplabel('Psychological','y');
         set(T, 'fontsize', fontSize+2, 'fontweight', 'normal');
         set(T(2), 'vert', 'top');


      case 'twoRepresentations'

         % relative directory structure in github repository
         fileList = {...
            '../perceptualReproduction/storage/perceptualReproductionMCMC_tomicBays2024_jags'; ...
            '../similarityComparison/storage/similarityComparison_tomicBays2024_jags'; ...
            };

         fontSize = 20;

         nModels = numel(fileList);
         [nRows, nCols] = subplotArrange(nModels);

         F = figure; clf; hold on;
         setFigure(F, [0.2 0.2 0.6 0.4], '');

         for modelIdx = 1:nModels
            fileName = fileList{modelIdx};

            fprintf('Loading pre-stored samples from file %s\n', fileName);
            % chains only — avoids stale similarityWrapChainInits handle warnings
            load(sprintf('%s', fileName), 'chains');

            % just keep converged chains
            switch modelIdx
               case 1 % perceptual
                  [keepChains, rHat] = findKeepChains(chains.sigma, 2, 1.1);
               case 2 % similarity
                  [keepChains, rHat] = findKeepChains(chains.sigma, 2, 1.1);
            end
            fields = fieldnames(chains);
            for i = 1:numel(fields)
               chains.(fields{i}) = chains.(fields{i})(:, keepChains);
            end

            mu = codatable(chains, 'mu', @mean);
            muTruth = ds.stimuli;

            xLim = [-1 1];
            subplot(nRows, nCols, modelIdx); cla; hold on;
            set(gca, ...
               'xlim'       , xLim    , ...
               'xtick'      , []   , ...
               'xticklabelrot', 0, ...
               'ylim'       , xLim    , ...
               'ytick'      , []  , ...
               'box'        , 'off'     , ...
               'tickdir'    , 'out'     , ...
               'layer'      , 'top'     , ...
               'ticklength' , [0.01 0]  , ...
               'layer'      , 'top'     , ...
               'clipping'   , 'off'     , ...
               'fontsize'   , fontSize  );
            moveAxis(gca, [1 1 0.9 1], [0 0 0 0]);
            axis square;
            axis off
            text(-1.15, 1.185, lower(char(64+modelIdx)), ...
               'fontsize', fontSize, 'fontweight', 'bold');

            r = [1.05 1.15];
            plot(xLim, [0 0], '-', ...
               'color', pantone.GlacierGray);
            plot(0, 0, 'o', ...
               'markerfacecolor', pantone.GlacierGray, ...
               'markeredgecolor', 'w', ...
               'markersize', 4);

            for idx = 1:ds.nStimuli
               plot([cos(muTruth(idx)) r(1)*cos(mu(idx))], [sin(muTruth(idx)) r(1)*sin(mu(idx))], '-', ...
                  'color', pantone.Titanium);
               plot(cos(muTruth(idx)), sin(muTruth(idx)), 'o', ...
                  'markerfacecolor',  'w', ...
                  'markeredgecolor', 'k', ...
                  'markersize', 5);
               plot(r*cos(mu(idx)), r*sin(mu(idx)), 'k-', ...
                  'linewidth', 2);
            end
         end


      case 'clusterSavageDickey'
         % DIC Mode A (~0, 1.5) with R̂-ok subset (prefer large pools first).
         chains = loadModeAClusterChains({ ...
            '../clusterRepresentation/storage/clusterRepresentation_tomicBays2024_jags', ...
            '../clusterRepresentation/storage/clusterRepresentationSelectMCMC_tomicBays2024_jags'});
         drawClusterSavageDickeyFigure(chains, pantone, [0 3], 1);

      case 'clusterSavageDickeyEntropification'
         % Entrop keep is the paper set (often not Mode A); zoom axes [0,1].
         fileName = '../clusterRepresentation/storage/clusterRepresentationEntropification_tomicBays2024_jags';
         fprintf('Loading keep from %s\n', fileName);
         S = load(fileName, 'chains', 'keepChains');
         chains = S.chains;
         if isfield(S, 'keepChains')
            fprintf('  keep=%s nCh=%d\n', mat2str(S.keepChains), size(chains.gammaPM, 2));
         end
         drawClusterSavageDickeyFigure(chains, pantone, [0 1], 0.5);

      case 'statisticalSummaries'

         fid = fopen('results/statisticalSummaries.txt', 'w');

         % Text-only analysis (does not change printFigures for later cases).

         % relative directory structure in github repository
         fileList = {...
            '../perceptualReproduction/storage/perceptualReproductionMCMC_tomicBays2024_jags'; ...
            '../memoryReproduction/storage/memoryReproductionMCMC_tomicBays2024_jags'; ...
            '../similarityComparison/storage/similarityComparison_tomicBays2024_jags'; ...
            '../clusterRepresentation/storage/clusterRepresentationSelectMCMC_tomicBays2024_jags'};
         % Cluster BF uses DIC Mode A (~0, 1.5); see loadModeAClusterChains.

         CIbounds = [2.5 97.5];

         % PERCEPTUAL
         fileName = fileList{1};
         fprintf('Loading pre-stored samples from file %s\n', fileName);
         load(sprintf('%s', fileName), 'chains');
         [keepChains, rHat] = findKeepChains(chains.sigma, 2, 1.1);
         fields = fieldnames(chains);
         for i = 1:numel(fields)
            chains.(fields{i}) = chains.(fields{i})(:, keepChains);
         end
         fprintf(fid, '\n----------\nPerceptual Reproduction\n\n');

         sigma = codatable(chains, 'sigma', @mean);
         bounds = prctile(chains.sigma(:), CIbounds);
         fprintf(fid, 'Posterior mean of sigma is %1.3f, with 95%% CI (%1.3f, %1.3f)\n', sigma, bounds);

         % MEMORY
         fileName = fileList{2};
         fprintf('\n----------\nLoading pre-stored samples from file %s\n', fileName);
         load(sprintf('%s', fileName), 'chains');
         [keepChains, rHat] = findKeepChains(chains.sigma_1, 2, 1.1);
         fields = fieldnames(chains);
         for i = 1:numel(fields)
            chains.(fields{i}) = chains.(fields{i})(:, keepChains);
         end
         fprintf(fid, '\n----------\nMemory Reproduction\n\n');

         % posterior summary for sigmas
         sigma3 = codatable(chains, 'sigma_1', @mean);
         bounds3 = prctile(chains.sigma_1(:), CIbounds);
         fprintf(fid, 'Posterior mean of sigma for set size 3 is %1.3f, with 95%% CI (%1.3f, %1.3f)\n', sigma3, bounds3);
         sigma6 = codatable(chains, 'sigma_2', @mean);
         bounds6 = prctile(chains.sigma_2(:), CIbounds);
         fprintf(fid, 'Posterior mean of sigma for set size 6 is %1.3f, with 95%% CI (%1.3f, %1.3f)\n', sigma6, bounds6);

         omega3 = get_matrix_from_coda(chains, 'omega3', @mean);
         omega6 = get_matrix_from_coda(chains, 'omega6', @mean);
         fprintf(fid, 'Posterior means for omega for set size 3 are (%1.3f, %1.3f, %1.3f)\n', omega3);
         fprintf(fid, 'Posterior means for omega for set size 6 are (%1.3f, %1.3f, %1.3f, %1.3f, %1.3f, %1.3f)\n', omega6);
         for i = 1:3
            bounds = prctile(chains.(sprintf('omega3_%d', i))(:), CIbounds);
            fprintf(fid, 'Posterior mean of omega_%d for set size 3 is %1.3f, with 95%% CI (%1.3f, %1.3f)\n', i, omega3(i), bounds);
         end
         for i = 1:6
            bounds = prctile(chains.(sprintf('omega6_%d', i))(:), CIbounds);
            fprintf(fid, 'Posterior mean of omega_%d for set size 6 is %1.3f, with 95%% CI (%1.3f, %1.3f)\n', i, omega6(i), bounds);
         end

         % bayes factor to no swap model at critical points
         omega3full = [chains.omega3_1(:) chains.omega3_2(:) chains.omega3_3(:)];
         omega6full = [chains.omega6_1(:) chains.omega6_2(:) chains.omega6_3(:) chains.omega6_4(:) chains.omega6_5(:) chains.omega6_6(:)];
         omega3fullPrior = [chains.omega3prior_1(:) chains.omega3prior_2(:) chains.omega3prior_3(:)];
         omega6fullPrior = [chains.omega6prior_1(:) chains.omega6prior_2(:) chains.omega6prior_3(:) chains.omega6prior_4(:) chains.omega6prior_5(:) chains.omega6prior_6(:)];
         omega3critical = [1 0 0];
         omega6critical = [1 0 0 0 0 0];

         % Log-space Gaussian KDE (same as supportingFiles/savageDickeyLogBF):
         % avoids -Inf when mvksdensity underflows at the no-swap corner.
         [logBF3, ~] = savageDickeyLogBF(omega3full, omega3fullPrior, omega3critical);
         [logBF6, ~] = savageDickeyLogBF(omega6full, omega6fullPrior, omega6critical);

         fprintf(fid, ['Savage-Dickey log BF_01 via multivariate KDE in log space ', ...
            '(positive favors null / no-swap; negative favors alternative / swap).\n']);
         fprintf(fid, '%s\n', formatLogBFSentence('memory set-size 3 no-swap null', logBF3));
         fprintf(fid, '%s\n', formatLogBFSentence('memory set-size 6 no-swap null', logBF6));

         % SIMILARITY
         fileName = fileList{3};
         fprintf('\n----------\nLoading pre-stored samples from file %s\n', fileName);
         load(sprintf('%s', fileName), 'chains');
         [keepChains, rHat] = findKeepChains(chains.sigma, 2, 1.1);
         fields = fieldnames(chains);
         for i = 1:numel(fields)
            chains.(fields{i}) = chains.(fields{i})(:, keepChains);
         end
         fprintf(fid, '\n----------\nSimilarity Comparison\n\n');

         sigma = codatable(chains, 'sigma', @mean);
         bounds = prctile(chains.sigma(:), CIbounds);
         fprintf(fid, 'Posterior mean of sigma is %1.3f, with 95%% CI (%1.3f, %1.3f)\n', sigma, bounds);

         % SIMILARITY + ENTROPIFICATION
         fileName = '../similarityComparison/storage/similarityComparisonEntropification_tomicBays2024_jags';
         fprintf('\n----------\nLoading pre-stored samples from file %s\n', fileName);
         load(sprintf('%s', fileName), 'chains');
         [keepChains, rHat] = findKeepChains(chains.sigma, 2, 1.1); %#ok<ASGLU>
         fields = fieldnames(chains);
         for i = 1:numel(fields)
            chains.(fields{i}) = chains.(fields{i})(:, keepChains);
         end
         fprintf(fid, '\n----------\nSimilarity Comparison Entropification\n\n');
         sigma = codatable(chains, 'sigma', @mean);
         bounds = prctile(chains.sigma(:), CIbounds);
         fprintf(fid, 'Posterior mean of sigma is %1.3f, with 95%% CI (%1.3f, %1.3f)\n', sigma, bounds);
         if isfield(chains, 'w')
            wMean = mean(chains.w(:));
            wBounds = prctile(chains.w(:), CIbounds);
            fprintf(fid, 'Posterior mean of w is %1.3f, with 95%% CI (%1.3f, %1.3f)\n', ...
               wMean, wBounds);
         end

         % CLUSTER (DIC Mode A; prefer large pool first)
         scale = 4;

         fprintf('\n----------\n');
         chains = loadModeAClusterChains({ ...
            '../clusterRepresentation/storage/clusterRepresentation_tomicBays2024_jags', ...
            fileList{4}});
         fprintf(fid, '\n----------\nCluster Bayes Factor (Mode A / DIC-preferred)\n\n');
         writeClusterGammaSummaries(fid, chains, CIbounds);
         writeClusterSavageDickeyBFs(fid, chains, scale);

         % CLUSTER + ENTROPIFICATION (paper keep)
         fprintf('\n----------\n');
         fileName = '../clusterRepresentation/storage/clusterRepresentationEntropification_tomicBays2024_jags';
         fprintf('Loading keep from %s\n', fileName);
         S = load(fileName, 'chains', 'keepChains');
         chains = S.chains;
         if isfield(S, 'keepChains')
            fprintf('  keep=%s nCh=%d\n', mat2str(S.keepChains), size(chains.gammaPM, 2));
         end
         fprintf(fid, '\n----------\nCluster Entropification (paper keep)\n\n');
         writeClusterGammaSummaries(fid, chains, CIbounds);
         if isfield(chains, 'w')
            wMean = mean(chains.w(:));
            wBounds = prctile(chains.w(:), CIbounds);
            fprintf(fid, 'Posterior mean of w is %1.3f, with 95%% CI (%1.3f, %1.3f)\n', ...
               wMean, wBounds);
         end
         writeClusterSavageDickeyBFs(fid, chains, scale);

         fclose(fid);


      case 'descriptiveAdequacyPerceptual'

         if exist('storage/descriptiveAdequacyPerceptualPreSave.mat', 'file')

            load('storage/descriptiveAdequacyPerceptualPreSave', 'credInterval', 'mn', 'binsC');

         else

            CI = [25 75];

            lo = 0; hi = pi; eps = 0.01;
            binsC = lo:eps:hi;
            binsE = lo-eps/2:eps:hi+eps/2;

            fileName = '../perceptualReproduction/storage/perceptualReproductionMCMC_tomicBays2024_jags';
            fprintf('Loading pre-stored samples from file %s\n', fileName);
            load(sprintf('%s', fileName), 'chains');

            dp.response = mod(dp.response, pi);

            [~, ~, rIdx] = histcounts(dp.response, 'binedges', binsE);

            mn = nan(numel(binsC), 1);
            credInterval = nan(numel(binsC), 2);
            for i = 1:length(binsC)
               match = find(rIdx == i);
               nPer = numel(chains.yP_1);   % e.g. 2000*12, or use 2000*12 explicitly
               vals = zeros(numel(match) * nPer, 1);
               idx = 1;
               for j = 1:numel(match)
                  k = match(j);
                  fname = sprintf('yP_%d', k);
                  block = chains.(fname)(:);
                  vals(idx:idx+numel(block)-1) = block;
                  idx = idx + numel(block);
               end
               k = round((binsC(i) - vals) ./ pi);
               vals = vals + k .* pi;
               mn(i) = mean(vals);
               credInterval(i, :) = prctile(vals, CI);
            end

            save('storage/descriptiveAdequacyPerceptualPreSave', 'mn', 'credInterval', 'binsC');

         end

         fontSize = 18;
         scale = 0.25;
         clr = pantone.DuskBlue;

         F = figure; clf; hold on;
         setFigure(F, [0.2 0.2 0.4 0.6], '');

         set(gca, ...
            'xlim'       , [0 pi]    , ...
            'xtick'      , [0 pi/4 pi/2 3*pi/4 pi]   , ...
            'xticklabelrot', 0, ...
            'xticklabel' , {'$0$', '$\frac{\pi}{4}$', '$\frac{\pi}{2}$', '$\frac{3\pi}{4}$', '$\pi$'}, ...
            'ylim'       , [-0.2 pi+0.2]    , ...
            'ytick'      , [ 0 pi/4 pi/2 3*pi/4 pi ]   , ...
            'yticklabel' , {'$0$', '$\frac{\pi}{4}$', '$\frac{\pi}{2}$', '$\frac{3\pi}{4}$', '$\pi$'}, ...
            'ticklabelinterpreter', 'latex', ...
            'box'        , 'off'     , ...
            'tickdir'    , 'out'     , ...
            'layer'      , 'top'     , ...
            'ticklength' , [0.02 0]  , ...
            'layer'      , 'top'     , ...
            'clipping'   , 'off'     , ...
            'fontsize'   , fontSize  );
         % axis square;
         xlabel('Data', 'fontsize', fontSize);
         ylabel('Posterior Predictive', 'fontsize', fontSize);
         moveAxis(gca, [1 1 0.9 0.9], [0.025 0.025 0 0]);
         Raxes(gca, 0.01, 0.03);

         for i = 1:length(binsC)
            plot([binsC(i) binsC(i)], credInterval(i, :), '-', ...
               'color', pantone.Custard);
         end
         for i = 1:length(binsC)

            plot(binsC(i), mn(i), 'o', ...
               'markerfacecolor', pantone.ClassicBlue, ...
               'markersize', 2, ...
               'markeredgecolor', pantone.ClassicBlue);
         end


      case 'descriptiveAdequacyMemory'

         if exist('storage/descriptiveAdequacyMemoryPreSave.mat', 'file')

            load('storage/descriptiveAdequacyMemoryPreSave', 'credInterval', 'mn', 'binsC');

         else

            CI = [25 75];

            lo = 0; hi = pi; eps = 0.01;
            binsC = lo:eps:hi;
            binsE = lo-eps/2:eps:hi+eps/2;

            % Needs yP; memoryReproductionMCMC omits posterior-predictive monitors.
            fileName = '../memoryReproduction/storage/memoryReproduction_tomicBays2024_jags';
            fprintf('Loading pre-stored samples from file %s\n', fileName);
            load(sprintf('%s', fileName), 'chains');

            dm.response = mod(dm.response, pi);

            [~, ~, rIdx] = histcounts(dm.response, 'binedges', binsE);

            mn = nan(numel(binsC), 1);
            credInterval = nan(numel(binsC), 2);
            for i = 1:length(binsC)
               match = find(rIdx == i);
               nPer = numel(chains.yP_1);   % e.g. 2000*12, or use 2000*12 explicitly
               vals = zeros(numel(match) * nPer, 1);
               idx = 1;
               for j = 1:numel(match)
                  k = match(j);
                  fname = sprintf('yP_%d', k);
                  block = chains.(fname)(:);
                  vals(idx:idx+numel(block)-1) = block;
                  idx = idx + numel(block);
               end
               k = round((binsC(i) - vals) ./ pi);
               vals = vals + k .* pi;
               mn(i) = mean(vals);
               credInterval(i, :) = prctile(vals, CI);
            end

            save('storage/descriptiveAdequacyMemoryPreSave', 'mn', 'credInterval', 'binsC');
         end

         fontSize = 18;
         scale = 0.25;
         clr = pantone.DuskBlue;

         F = figure; clf; hold on;
         setFigure(F, [0.2 0.2 0.4 0.6], '');

         set(gca, ...
            'xlim'       , [0 pi]    , ...
            'xtick'      , [0 pi/4 pi/2 3*pi/4 pi]   , ...
            'xticklabelrot', 0, ...
            'xticklabel' , {'$0$', '$\frac{\pi}{4}$', '$\frac{\pi}{2}$', '$\frac{3\pi}{4}$', '$\pi$'}, ...
            'ylim'       , [-0.8 pi+0.8]    , ...
            'ytick'      , [ 0 pi/4 pi/2 3*pi/4 pi ]   , ...
            'yticklabel' , {'$0$', '$\frac{\pi}{4}$', '$\frac{\pi}{2}$', '$\frac{3\pi}{4}$', '$\pi$'}, ...
            'ticklabelinterpreter', 'latex', ...
            'box'        , 'off'     , ...
            'tickdir'    , 'out'     , ...
            'layer'      , 'top'     , ...
            'ticklength' , [0.02 0]  , ...
            'layer'      , 'top'     , ...
            'clipping'   , 'off'     , ...
            'fontsize'   , fontSize  );
         xlabel('Data', 'fontsize', fontSize);
         ylabel('Posterior Predictive', 'fontsize', fontSize);
         moveAxis(gca, [1 1 0.9 0.9], [0.025 0.025 0 0]);
         Raxes(gca, 0.01, 0.03);


         for i = 1:length(binsC)
            plot([binsC(i) binsC(i)], credInterval(i, :), '-', ...
               'color', pantone.Custard);
         end
         for i = 1:length(binsC)

            plot(binsC(i), mn(i), 'o', ...
               'markerfacecolor', pantone.ClassicBlue, ...
               'markersize', 2, ...
               'markeredgecolor', pantone.ClassicBlue);
         end

      case 'memoryNoSwap'
         % SI: physical vs psychological μ under the no-swap memory model.
         fileName = '../memoryReproduction/storage/memoryReproductionNoSwapMCMC_tomicBays2024_jags';
         fprintf('Loading pre-stored samples from file %s\n', fileName);
         load(sprintf('%s', fileName), 'chains');
         [keepChains, ~] = findKeepChains(chains.sigma_1, 2, 1.1);
         fields = fieldnames(chains);
         for i = 1:numel(fields)
            chains.(fields{i}) = chains.(fields{i})(:, keepChains);
         end

         fontSize = 20;
         CIbounds = [2.5 97.5];
         mu = nan(dm.nStimuli, 1);
         muBounds = nan(dm.nStimuli, 2);
         for idx = 1:dm.nStimuli
            vals = chains.(sprintf('mu_%d', idx))(:);
            [mu(idx), muBounds(idx, :)] = summarizeHalfCircle(vals, CIbounds);
         end
         muTruth = dp.stimuli;

         F = figure; clf; hold on;
         setFigure(F, [0.2 0.2 0.35 0.5], '');
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
            'clipping'   , 'off'     , ...
            'fontsize'   , fontSize  );
         axis square;
         xlabel('Physical', 'fontsize', fontSize);
         ylabel('Psychological', 'fontsize', fontSize);
         text(0, pi, 'memory (no swap)', ...
            'fontsize', fontSize-4, ...
            'fontweight', 'normal', ...
            'vert', 'bot', 'hor', 'lef');
         moveAxis(gca, [1 1 0.95 0.95], [0 0.025 0 0]);
         Raxes(gca, 0.02, 0.01);

         for i = pi/4:pi/4:3*pi/4
            plot([i i], [0 pi], '-', 'color', pantone.GlacierGray);
            plot([0 pi], [i i], '-', 'color', pantone.GlacierGray);
         end
         for idx = 1:dp.nStimuli
            if muBounds(idx, 1) > muBounds(idx, 2)
               plot(muTruth(idx) * [1 1], [0 muBounds(idx, 2)], '-', ...
                  'color', pantone.ClassicBlue, 'linewidth', 1);
               plot(muTruth(idx) * [1 1], [muBounds(idx, 1) pi], '-', ...
                  'color', pantone.ClassicBlue, 'linewidth', 1);
            else
               plot(muTruth(idx) * [1 1], muBounds(idx, :), '-', ...
                  'color', pantone.ClassicBlue, 'linewidth', 1);
            end
            plot(muTruth(idx), mu(idx), 'o', ...
               'markerfacecolor', pantone.ClassicBlue, ...
               'markeredgecolor', 'w', 'linewidth', 0.5, 'markersize', 4);
         end
         plot([0 pi], [0 pi], '-', ...
            'color', pantone.AuroraRed, 'linewidth', 0.5);

      case 'similarityEntropification'
         % SI: physical vs psychological μ under similarity + entropification.
         fileName = '../similarityComparison/storage/similarityComparisonEntropification_tomicBays2024_jags';
         fprintf('Loading pre-stored samples from file %s\n', fileName);
         load(sprintf('%s', fileName), 'chains');
         [keepChains, ~] = findKeepChains(chains.sigma, 2, 1.1);
         fields = fieldnames(chains);
         for i = 1:numel(fields)
            chains.(fields{i}) = chains.(fields{i})(:, keepChains);
         end

         fontSize = 20;
         CIbounds = [2.5 97.5];
         mu = nan(ds.nStimuli, 1);
         muBounds = nan(ds.nStimuli, 2);
         for idx = 1:ds.nStimuli
            vals = chains.(sprintf('mu_%d', idx))(:);
            [mu(idx), muBounds(idx, :)] = summarizeHalfCircle(vals, CIbounds);
         end
         muTruth = dp.stimuli;

         F = figure; clf; hold on;
         setFigure(F, [0.2 0.2 0.35 0.5], '');
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
            'clipping'   , 'off'     , ...
            'fontsize'   , fontSize  );
         axis square;
         xlabel('Physical', 'fontsize', fontSize);
         ylabel('Psychological', 'fontsize', fontSize);
         text(0, pi, 'similarity (entropification)', ...
            'fontsize', fontSize-4, ...
            'fontweight', 'normal', ...
            'vert', 'bot', 'hor', 'lef');
         moveAxis(gca, [1 1 0.95 0.95], [0 0.025 0 0]);
         Raxes(gca, 0.02, 0.01);

         for i = pi/4:pi/4:3*pi/4
            plot([i i], [0 pi], '-', 'color', pantone.GlacierGray);
            plot([0 pi], [i i], '-', 'color', pantone.GlacierGray);
         end
         for idx = 1:dp.nStimuli
            if muBounds(idx, 1) > muBounds(idx, 2)
               plot(muTruth(idx) * [1 1], [0 muBounds(idx, 2)], '-', ...
                  'color', pantone.ClassicBlue, 'linewidth', 1);
               plot(muTruth(idx) * [1 1], [muBounds(idx, 1) pi], '-', ...
                  'color', pantone.ClassicBlue, 'linewidth', 1);
            else
               plot(muTruth(idx) * [1 1], muBounds(idx, :), '-', ...
                  'color', pantone.ClassicBlue, 'linewidth', 1);
            end
            plot(muTruth(idx), mu(idx), 'o', ...
               'markerfacecolor', pantone.ClassicBlue, ...
               'markeredgecolor', 'w', 'linewidth', 0.5, 'markersize', 4);
         end
         plot([0 pi], [0 pi], '-', ...
            'color', pantone.AuroraRed, 'linewidth', 0.5);

      case 'descriptiveAdequacyMemoryNoSwap'

         if exist('storage/descriptiveAdequacyMemoryNoSwapPreSave.mat', 'file')

            load('storage/descriptiveAdequacyMemoryNoSwapPreSave', 'credInterval', 'mn', 'binsC');

         else

            CI = [25 75];

            lo = 0; hi = pi; eps = 0.01;
            binsC = lo:eps:hi;
            binsE = lo-eps/2:eps:hi+eps/2;

            % Needs yP; NoSwap MCMC omits it — run memoryReproductionNoSwap.m on 2024.
            fileName = '../memoryReproduction/storage/memoryReproductionNoSwap_tomicBays2024_jags';
            if ~isfile([fileName '.mat'])
               error('drawFigures:noSwapYPmissing', ...
                  ['Missing %s (needs yP). Use SI case ''memoryNoSwap'' for the ', ...
                  'μ scatter from NoSwap MCMC, or fit memoryReproductionNoSwap.m.'], ...
                  [fileName '.mat']);
            end
            fprintf('Loading pre-stored samples from file %s\n', fileName);
            load(sprintf('%s', fileName), 'chains');

            dm.response = mod(dm.response, pi);

            [~, ~, rIdx] = histcounts(dm.response, 'binedges', binsE);

            mn = nan(numel(binsC), 1);
            credInterval = nan(numel(binsC), 2);
            for i = 1:length(binsC)
               match = find(rIdx == i);
               nPer = numel(chains.yP_1);   % e.g. 2000*12, or use 2000*12 explicitly
               vals = zeros(numel(match) * nPer, 1);
               idx = 1;
               for j = 1:numel(match)
                  k = match(j);
                  fname = sprintf('yP_%d', k);
                  block = chains.(fname)(:);
                  vals(idx:idx+numel(block)-1) = block;
                  idx = idx + numel(block);
               end
               k = round((binsC(i) - vals) ./ pi);
               vals = vals + k .* pi;
               mn(i) = mean(vals);
               credInterval(i, :) = prctile(vals, CI);
            end

            save('storage/descriptiveAdequacyMemoryNoSwapPreSave', 'mn', 'credInterval', 'binsC');
         end

         fontSize = 18;
         scale = 0.25;
         clr = pantone.DuskBlue;

         F = figure; clf; hold on;
         setFigure(F, [0.2 0.2 0.4 0.6], '');

         set(gca, ...
            'xlim'       , [0 pi]    , ...
            'xtick'      , [0 pi/4 pi/2 3*pi/4 pi]   , ...
            'xticklabelrot', 0, ...
            'xticklabel' , {'$0$', '$\frac{\pi}{4}$', '$\frac{\pi}{2}$', '$\frac{3\pi}{4}$', '$\pi$'}, ...
            'ylim'       , [-0.8 pi+0.8]    , ...
            'ytick'      , [ 0 pi/4 pi/2 3*pi/4 pi ]   , ...
            'yticklabel' , {'$0$', '$\frac{\pi}{4}$', '$\frac{\pi}{2}$', '$\frac{3\pi}{4}$', '$\pi$'}, ...
            'ticklabelinterpreter', 'latex', ...
            'box'        , 'off'     , ...
            'tickdir'    , 'out'     , ...
            'layer'      , 'top'     , ...
            'ticklength' , [0.02 0]  , ...
            'layer'      , 'top'     , ...
            'clipping'   , 'off'     , ...
            'fontsize'   , fontSize  );
         xlabel('Data', 'fontsize', fontSize);
         ylabel('Posterior Predictive', 'fontsize', fontSize);
         moveAxis(gca, [1 1 0.9 0.9], [0.025 0.025 0 0]);
         Raxes(gca, 0.01, 0.03);


         for i = 1:length(binsC)
            plot([binsC(i) binsC(i)], credInterval(i, :), '-', ...
               'color', pantone.Custard);
         end
         for i = 1:length(binsC)

            plot(binsC(i), mn(i), 'o', ...
               'markerfacecolor', pantone.ClassicBlue, ...
               'markersize', 2, ...
               'markeredgecolor', pantone.ClassicBlue);
         end

      case 'descriptiveAdequacySimilarity'

         nReps = 1e3;

         fileName = '../similarityComparison/storage/similarityComparison_tomicBays2024_jags';
         fprintf('Loading pre-stored samples from file %s\n', fileName);
         load(sprintf('%s', fileName), 'chains');

         if exist('storage/descriptiveAdequacySimilarityPreSave.mat', 'file')

            load('storage/descriptiveAdequacySimilarityPreSave', 'pC');

         else

            nSamples = 1000*8;

            % need to sample posterior predictive because of censoring
            % base it on just posterior means of mu and sigma
            mu = zeros(ds.nStimuli, nSamples);
            sigma = chains.sigma(:);
            for i = 1:ds.nStimuli
               mu(i, :) = chains.(sprintf('mu_%d', i))(:);
            end

            muMean = codatable(chains, 'mu', @mean);
            sigmaMean = codatable(chains, 'sigma', @mean);

            pC = zeros(ds.nTrials, 1);
            for i = 1:ds.nTrials
               yP = zeros(nReps, 1);
               for rep = 1:nReps
                  xA = randn*sigma(ceil(rand*nSamples)) + mu(ds.aIdx(i), ceil(rand*nSamples));
                  xB = randn*sigma(ceil(rand*nSamples)) + mu(ds.bIdx(i), ceil(rand*nSamples));
                  xC = randn*sigma(ceil(rand*nSamples)) + mu(ds.cIdx(i), ceil(rand*nSamples));
                  xD = randn*sigma(ceil(rand*nSamples)) + mu(ds.dIdx(i), ceil(rand*nSamples));

                  dAB = min(abs(xA-xB), pi-abs(xA-xB));
                  dCD = min(abs(xC-xD), pi-abs(xC-xD));

                  yP(rep) = (dAB > dCD);
               end
               if ds.response(i) == 1
                  pC(i) = mean(yP);
               else
                  pC(i) = 1 - mean(yP);
               end
               %mean(pC(1:i))
            end

            save('storage/descriptiveAdequacySimilarityPreSave', 'pC');
         end

         fontSize = 18;
         lo = 0; hi = 1; eps = 0.02; xTick = 0.1;
         clr = pantone.ClassicBlue;

         binsC = lo:eps:hi;
         binsE = lo-eps/2:eps:hi+eps/2;

         F = figure; clf; hold on;
         setFigure(F, [0.2 0.2 0.4 0.6], '');

         set(gca, ...
            'xlim'       , [lo hi]    , ...
            'xtick'      , lo:xTick:hi   , ...
            'xticklabelrot', 0, ...
            'ycolor', 'none', ...
            'box'        , 'off'     , ...
            'tickdir'    , 'out'     , ...
            'layer'      , 'top'     , ...
            'ticklength' , [0.02 0]  , ...
            'layer'      , 'top'     , ...
            'clipping'   , 'off'     , ...
            'fontsize'   , fontSize  );
         xlabel('Posterior Predictive Agreement', 'fontsize', fontSize);
         moveAxis(gca, [1 1 1 1], [0 0.025 0 0]);
         Raxes(gca, 0.01, 0.01);

         count = histcounts(pC, 'binedges', binsE, 'norm', 'prob');
         H = bar(binsC, count);
         set(H, 'facecolor', clr, 'edgecolor', 'w', 'barwidth', 0.9);

         fprintf('\n mean agreement is %1.2f, with interquartile range %1.2f-%1.2f\n', ...
            mean(pC), prctile(pC, [2.5 97.5]));
         fprintf('\n median agreement is %d\n', ...
            round(100*median(pC)));

         fprintf('\n Forced choice agreement is %1.2f\n', mean(pC >= 0.5));

      case 'swapExamples'

         % Need trial-level xi. memoryReproductionMCMC deliberately omits xi/yP
         % monitors; use the plain driver product that stores them.
         fileName = '../memoryReproduction/storage/memoryReproduction_tomicBays2024_jags';
         fprintf('Loading pre-stored samples from file %s\n', fileName);
         load(sprintf('%s', fileName), 'chains');
         if ~isfield(chains, 'xi_1_1')
            error('drawFigures:swapExamplesNoXi', ...
               ['%s has no xi_* fields (MCMC mats omit them). ', ...
               'Re-run memoryReproduction.m on tomicBays2024.'], fileName);
         end

         [keepChains, ~] = findKeepChains(chains.sigma_1, 2, 1.1);
         fields = fieldnames(chains);
         for i = 1:numel(fields)
            chains.(fields{i}) = chains.(fields{i})(:, keepChains);
         end

         [T, picked] = buildSwapExampleTable(dm, chains, pi);
         if ~isfolder('figures')
            mkdir('figures');
         end
         writetable(T, 'figures/swapExampleCandidates.csv');
         writetable(picked, 'figures/swapExamplePicked.csv');

         fprintf('Visually clear candidates: %d (%d setSize=3, %d setSize=6)\n', ...
            height(T), sum(T.setSize == 3), sum(T.setSize == 6));
         fprintf('Picked %d trials (%d setSize=3, %d setSize=6)\n', ...
            height(picked), sum(picked.setSize == 3), sum(picked.setSize == 6));
         if ~isempty(picked)
            disp(picked(:, {'trial', 'setSize', 'pBestFoil', 'bestFoilXi', ...
               'dInferredFoil', 'foilMargin', 'participant'}));
            F = plotSwapExamples(picked, dm, pantone, pi);
            figure(F);
         else
            warning('swapExamples:NoTrials', ...
               'No swap-example trials passed filters; CSV tables written, figure skipped.');
         end

   end

   % print (skip text-only analyses with no figure)
   if printFigures && ~strcmp(analysisName, 'statisticalSummaries')
      if ~isfolder('figures')
         !mkdir figures
      end
      print(sprintf('figures/%s.png', analysisName), '-dpng');
      print(sprintf('figures/%s.eps', analysisName), '-depsc');
   end

end

function writeClusterGammaSummaries(fid, chains, CIbounds)
gPM = mean(chains.gammaPM(:));
gPS = mean(chains.gammaPS(:));
bPM = prctile(chains.gammaPM(:), CIbounds);
bPS = prctile(chains.gammaPS(:), CIbounds);
fprintf(fid, 'Posterior mean of gamma_pm is %1.3f, with 95%% CI (%1.3f, %1.3f)\n', ...
   gPM, bPM);
fprintf(fid, 'Posterior mean of gamma_ps is %1.3f, with 95%% CI (%1.3f, %1.3f)\n', ...
   gPS, bPS);
end

function writeClusterSavageDickeyBFs(fid, chains, scale)
%WRITECLUSTERSAVAGEDICKEYBFS Joint and marginal nulls at gamma = 0.
%
% Reports log BF_01 = log f(θ0|data) - log f(θ0) (Savage–Dickey).
% Positive => favors null (common map); negative => favors alternative.
fprintf(fid, ['Savage-Dickey log BF_01 (positive favors null / common map; ', ...
   'negative favors alternative / separate maps). ', ...
   'Densities via Gaussian KDE in log space (bandwidth x%d).\n'], scale);

Xposterior = [chains.gammaPS(:) chains.gammaPM(:)];
g0 = chains.gammaPrior(:);
Xprior = [g0, g0(randperm(numel(g0)))];
logBF = kdeLogBF(Xposterior, Xprior, [0 0], scale);
fprintf(fid, '%s\n', formatLogBFSentence('joint null (gamma_pm, gamma_ps)=(0,0)', logBF));

logBF = kdeLogBF(chains.gammaPM(:), chains.gammaPrior(:), 0, scale);
fprintf(fid, '%s\n', formatLogBFSentence('perceptual-memory null (gamma_pm=0)', logBF));

logBF = kdeLogBF(chains.gammaPS(:), chains.gammaPrior(:), 0, scale);
fprintf(fid, '%s\n', formatLogBFSentence('perceptual-similarity null (gamma_ps=0)', logBF));
end

function logBF = kdeLogBF(Xposterior, Xprior, Xcrit, scale)
%KDELOGBF log BF_01 via log-space Gaussian KDE (finite; no -Inf underflow).
if isvector(Xposterior), Xposterior = Xposterior(:); end
if isvector(Xprior), Xprior = Xprior(:); end
Xcrit = double(Xcrit(:)');
logFPost = logKdeDensity(Xposterior, Xcrit, scale);
logFPrior = logKdeDensity(Xprior, Xcrit, scale);
logBF = logFPost - logFPrior;
end

function logF = logKdeDensity(draws, critical, scale)
n = size(draws, 1);
d = size(draws, 2);
sigma = std(draws, 0, 1);
bw = sigma * (4 / (n * (d + 4))) ^ (1 / (d + 4));
bw(bw <= 0 | ~isfinite(bw)) = 1e-3;
bw = max(bw, 1e-3) * scale;
logNorm = -sum(log(bw)) - 0.5 * d * log(2 * pi);
z = (draws - critical) ./ bw;
logKi = logNorm - 0.5 * sum(z.^2, 2);
m = max(logKi);
logF = m + log(sum(exp(logKi - m))) - log(n);
if ~isfinite(logF)
   [~, j] = min(sum(z.^2, 2));
   logF = logKi(j) - log(n);
end
end

function logBF = densityLogBF(fPosterior, fPrior)
%DENSITYLOGBF log(fPost/fPrior); ±Inf when a density is zero/nonfinite.
if ~(isfinite(fPosterior) && fPosterior > 0)
   if ~(isfinite(fPrior) && fPrior > 0)
      logBF = NaN;
   else
      logBF = -Inf;
   end
elseif ~(isfinite(fPrior) && fPrior > 0)
   logBF = Inf;
else
   logBF = log(fPosterior) - log(fPrior);
end
end

function str = formatLogBF(logBF, bound)
%FORMATLOGBF Compact printable log BF_01 (used by memory swap BF lines).
if nargin < 2
   bound = 100;
end
if ~isfinite(logBF)
   if isnan(logBF)
      str = 'undefined';
   elseif logBF < 0
      str = sprintf('< -%d', bound);
   else
      str = sprintf('> +%d', bound);
   end
else
   str = sprintf('= %.1f', logBF);
end
end

function str = formatLogBFSentence(label, logBF01)
%FORMATLOGBFSENTENCE Clear null vs alternative wording for summaries.
isMemory = contains(lower(char(label)), 'memory') || contains(lower(char(label)), 'swap');
if isMemory
   nullName = 'NULL / no-swap';
   altName = 'ALTERNATIVE / swap';
else
   nullName = 'NULL / common map';
   altName = 'ALTERNATIVE / separate maps';
end
if ~isfinite(logBF01)
   if isnan(logBF01)
      str = sprintf('  %s: log BF_01 undefined', label);
   elseif logBF01 < 0
      str = sprintf('  %s: log BF_01 -> -Inf (favors %s); log BF_10 -> +Inf', ...
         label, altName);
   else
      str = sprintf('  %s: log BF_01 -> +Inf (favors %s); log BF_10 -> -Inf', ...
         label, nullName);
   end
   return
end
logBF10 = -logBF01;
if logBF01 > 0
   favored = nullName;
elseif logBF01 < 0
   favored = altName;
else
   favored = 'neither (BF_01 = 1)';
end
str = sprintf(['  %s: log BF_01 = %.2f (favors %s); ', ...
   'equivalently log BF_10 = %.2f'], label, logBF01, favored, logBF10);
end

function chains = loadModeAClusterChains(fileNames)
%LOADMODEACLUSTERCHAINS Converged DIC Mode A chains ≈ (0, 1.5).
%
% Mode A region (chain medians): gPM < 0.8 and gPS in [0.8, 2.2].
% Then pick the largest subset with R̂ ≤ 1.1 on gammaPM/gammaPS/gammaPrior
% (same idea as selectClusterModeChains). Does not use muS-band keep alone.
% fileNames: char/string or cell of candidates. Falls back to DIC Mode-A
% seeded run only if that subset also passes R̂.

% Mode A box matches countModeAKeep / merged-pool scan (≈(0,1.5)).
% Main cluster pool has ~62 chains here vs ~57 at Mode B (≈(2,1.5)); R̂≤1.05
% finds none, but R̂≤1.1 yields a keep of 5 near (0.12, 1.48).
modeA_pmMax = 0.8;
modeA_psLo = 0.8;
modeA_psHi = 2.2;
minKeep = 4;
maxKeep = 8;
rHatMax = 1.1;
rHatFields = {'gammaPM', 'gammaPS', 'gammaPrior'};
modeA_center = [0.1, 1.5];

if ischar(fileNames) || isstring(fileNames)
   fileNames = {char(fileNames)};
end

for f = 1:numel(fileNames)
   fileName = fileNames{f};
   matPath = fileName;
   if ~endsWith(matPath, '.mat')
      matPath = [matPath '.mat']; %#ok<AGROW>
   end
   if ~isfile(matPath)
      fprintf('Skipping missing %s\n', matPath);
      continue
   end
   fprintf(['Loading %s\n  Mode A region gPM<%.1f, gPS in [%.1f, %.1f]; ', ...
      'then Rhat<=%.2f, keep %d-%d\n'], ...
      matPath, modeA_pmMax, modeA_psLo, modeA_psHi, rHatMax, minKeep, maxKeep);
   S = load(matPath);

   src = [];
   if isfield(S, 'chainsAll') && isstruct(S.chainsAll) ...
         && isfield(S.chainsAll, 'gammaPM')
      src = S.chainsAll;
   elseif isfield(S, 'chains') && isstruct(S.chains) ...
         && isfield(S.chains, 'gammaPM')
      src = S.chains;
   end
   if isempty(src)
      continue
   end

   mPM = median(src.gammaPM, 1);
   mPS = median(src.gammaPS, 1);
   modeA = find(mPM < modeA_pmMax & mPS >= modeA_psLo & mPS <= modeA_psHi);
   fprintf('  Mode A candidates among %d chains: %s\n', ...
      size(src.gammaPM, 2), mat2str(modeA));
   if numel(modeA) < minKeep
      fprintf('  Too few Mode A chains; trying next source.\n');
      continue
   end
   % Prefer chains closest to DIC Mode A center before R̂ search.
   dCenter = (mPM(modeA) - modeA_center(1)).^2 + (mPS(modeA) - modeA_center(2)).^2;
   [~, ord] = sort(dCenter);
   modeA = modeA(ord);
   fprintf('  Mode A ordered by distance to (%.1f, %.1f): %s\n', ...
      modeA_center(1), modeA_center(2), mat2str(modeA));

   [keep, worstR] = pickRhatOkSubset(src, modeA, minKeep, maxKeep, ...
      rHatMax, rHatFields);
   if isempty(keep)
      fprintf('  No Rhat-ok Mode A subset; trying next source.\n');
      continue
   end
   chains = subsetChainColumns(src, keep);
   fprintf('  Keep %s (worst Rhat=%.3f); medians gPM=%.3f gPS=%.3f\n', ...
      mat2str(keep), worstR, median(chains.gammaPM(:)), median(chains.gammaPS(:)));
   return
end

fprintf('  No Rhat-ok Mode A in candidates; trying DIC compare fallback.\n');
chains = loadDicModeAFallback(minKeep, rHatMax, rHatFields);
end

function chains = loadDicModeAFallback(minKeep, rHatMax, rHatFields)
dicFile = '../clusterRepresentation/storage/clusterGammaModeCompare_tomicBays2024_jags';
fprintf('  Fallback: %s (A_nearZeroPM)\n', dicFile);
C = load(dicFile, 'results');
idx = find(strcmp({C.results.name}, 'A_nearZeroPM'), 1);
if isempty(idx) || ~isfield(C.results(idx), 'chains')
   error('drawFigures:noModeA', ...
      'No Mode A chains in weekend mats or DIC compare fallback.');
end
src = C.results(idx).chains;
nCh = size(src.gammaPM, 2);
[keep, worstR] = pickRhatOkSubset(src, 1:nCh, min(minKeep, nCh), nCh, ...
   rHatMax, rHatFields);
if isempty(keep)
   error('drawFigures:noModeAConverge', ...
      ['DIC Mode A fallback has no subset with Rhat <= %.2f on %s. ', ...
      'Need more Mode A sampling before SD figures.'], ...
      rHatMax, strjoin(rHatFields, ', '));
end
chains = subsetChainColumns(src, keep);
fprintf('  DIC Mode A keep %s (worst Rhat=%.3f); medians gPM=%.3f gPS=%.3f (DIC=%.1f)\n', ...
   mat2str(keep), worstR, median(chains.gammaPM(:)), median(chains.gammaPS(:)), ...
   C.results(idx).DIC);
end

function [keep, worstR] = pickRhatOkSubset(chains, band, minKeep, maxKeep, rHatMax, rHatFields)
% Largest subset of band with gelmanrubin Rhat <= rHatMax on rHatFields.
% `band` should already be preference-ordered (best first). Cap combinatorial
% search at the first maxCand members so nchoosek stays tractable.
keep = [];
worstR = inf;
band = band(:)';
if numel(band) < minKeep
   return
end
maxCand = 12;
if numel(band) > maxCand
   fprintf('  Restricting Rhat search to first %d of %d Mode A candidates\n', ...
      maxCand, numel(band));
   band = band(1:maxCand);
end
for n = min(maxKeep, numel(band)):-1:minKeep
   combos = nchoosek(band, n);
   bestWorst = inf;
   bestIdx = [];
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
         if ~isfinite(r) || r > rHatMax
            ok = false;
            break
         end
         worst = max(worst, r);
      end
      if ok && worst < bestWorst
         bestWorst = worst;
         bestIdx = sort(idx);
      end
   end
   if ~isempty(bestIdx)
      keep = bestIdx;
      worstR = bestWorst;
      return
   end
end
end

function chains = subsetChainColumns(src, idx)
fn = fieldnames(src);
chains = struct();
for i = 1:numel(fn)
   chains.(fn{i}) = src.(fn{i})(:, idx);
end
end

function drawClusterSavageDickeyFigure(chains, pantone, lim, gammaTick)
%DRAWCLUSTERSAVAGEDICKEYFIGURE Prior/posterior γ scatter on lim = [lo hi].
gammaLo = lim(1);
gammaHi = lim(2);
nBin = 40;
gammaE = linspace(gammaLo, gammaHi, nBin + 1);
gammaC = (gammaE(1:end-1) + gammaE(2:end)) / 2;

countPosterior = histcounts2(chains.gammaPM(:), chains.gammaPS(:), ...
   gammaE, gammaE, 'normalization', 'probability');
count = histcounts(chains.gammaPrior(:), gammaE);
if ~isempty(count) && sum(count) > 0
   count(1) = 2 * count(1);
   count = count / sum(count);
else
   count = ones(size(gammaC)) / numel(gammaC);
end
countPrior = count' * count;

fontSize = 18;
scale = 0.5 * (gammaHi - gammaLo) / 3;
colorPrior = pantone.Custard;
colorPosterior = pantone.ClassicBlue;

F = figure; clf; hold on;
setFigure(F, [0.2 0.2 0.4 0.4], '');

set(gca, ...
   'xlim'       , [gammaLo gammaHi]    , ...
   'xtick'      , gammaLo:gammaTick:gammaHi    , ...
   'ylim'       , [gammaLo gammaHi]    , ...
   'ytick'      , gammaLo:gammaTick:gammaHi    , ...
   'box'        , 'off'     , ...
   'tickdir'    , 'out'     , ...
   'layer'      , 'top'     , ...
   'ticklength' , [0.02 0]  , ...
   'clipping'   , 'off'     , ...
   'fontsize'   , fontSize  );
axis square;
xlabel('$\gamma_{\mathrm{pm}}$', 'fontsize', fontSize+4, 'interp', 'latex');
ylabel('$\gamma_{\mathrm{ps}}$', 'fontsize', fontSize+4, 'interp', 'latex');
moveAxis(gca, [1 1 0.9 0.9], [0 0.025 0 0]);
Raxes(gca, 0.01, 0.01);

% Posterior first, prior on top so prior remains visible where they overlap.
for i = 1:length(gammaC)
   for j = 1:length(gammaC)
      val = sqrt(countPosterior(i, j)) * scale;
      if val > 0
         rectangle('position', [gammaC(i)-val/2 gammaC(j)-val/2 val val], ...
            'curvature', [1 1], ...
            'facecolor', colorPosterior, ...
            'edgecolor', colorPosterior);
      end
   end
end
for i = 1:length(gammaC)
   for j = 1:length(gammaC)
      val = sqrt(countPrior(i, j)) * scale;
      if val > 0
         rectangle('position', [gammaC(i)-val/2 gammaC(j)-val/2 val val], ...
            'curvature', [1 1], ...
            'facecolor', colorPrior, ...
            'edgecolor', colorPrior);
      end
   end
end

H(1) = plot(-1, -1, 'o', ...
   'markersize', 6, ...
   'markerfacecolor', colorPrior, ...
   'markeredgecolor', colorPrior);
H(2) = plot(-1, -1, 'o', ...
   'markersize', 6, ...
   'markerfacecolor', colorPosterior, ...
   'markeredgecolor', colorPosterior);

L = legend(H, {'prior', 'posterior'}, ...
   'fontsize', fontSize, ...
   'box', 'off', ...
   'location', 'northeast');
set(L, 'position', get(L, 'position') + [0.075 0.05 0 0]);
end

function [muMean, bounds] = summarizeHalfCircle(vals, CI)
vals = vals(:);
phi = mod(2 * vals, 2 * pi);
phi0 = atan2(mean(sin(phi)), mean(cos(phi)));
if phi0 < 0
   phi0 = phi0 + 2 * pi;
end
rel = angle(exp(1i * (phi - phi0)));
relBounds = prctile(rel, CI);
muMean = mod(phi0 / 2, pi);
bounds = mod(muMean + relBounds / 2, pi);
end

function [T, picked] = buildSwapExampleTable(dm, chains, piVal)
angDist = @(a, b) min(abs(a - b), piVal - abs(a - b));
nTrials = dm.nTrials;
[~, ~, setSizeIdx] = unique(dm.setSize, 'stable');

rows = struct( ...
   'trial', {}, 'participant', {}, 'setSize', {}, ...
   'pSwap', {}, 'pTarget', {}, 'pBestFoil', {}, 'bestXi', {}, 'bestFoilXi', {}, ...
   'visualFoil', {}, 'visualMatch', {}, ...
   'resp', {}, 'target', {}, 'bestFoil', {}, 'bestFoilIdx', {}, ...
   'dTarget', {}, 'dBestFoil', {}, 'dInferredFoil', {}, ...
   'foilMargin', {}, 'score', {});

for t = 1:nTrials
   ss = setSizeIdx(t);
   xiSamples = chains.(sprintf('xi_%d_%d', t, ss))(:);
   pTarget = mean(xiSamples == 1);
   pSwap = 1 - pTarget;
   nCat = dm.setSize(t);
   counts = histcounts(xiSamples, 0.5:(nCat + 0.5));
   pCat = counts / sum(counts);
   [~, bestXi] = max(pCat);
   [pBestFoil, foilSlot] = max(pCat(2:end));
   bestFoilXi = foilSlot + 1;
   inferredFoilNum = bestFoilXi - 1;

   resp = mod(dm.response(t), piVal);
   targ = mod(dm.target(t), piVal);
   nFoils = dm.setSize(t) - 1;
   foilsRaw = mod(dm.nontarget(t, 1:nFoils), piVal);
   foilOrder = sortFoilsByTargetDistance(targ, foilsRaw, piVal);
   numberedFoils = foilsRaw(foilOrder);

   dT = angDist(resp, targ);
   dToNumbered = arrayfun(@(f) angDist(resp, numberedFoils(f)), 1:nFoils);
   [dBestFoil, visualFoil] = min(dToNumbered);
   foilMargin = dT - dBestFoil;
   dInferredFoil = dToNumbered(inferredFoilNum);
   visualMatch = visualFoil == inferredFoilNum;

   if dm.setSize(t) == 3
      minMargin = 0.8;
      maxDInferred = 0.35;
      minPFoil = 0.33;
   elseif bestFoilXi ~= 6
      minMargin = 0.85;
      maxDInferred = 0.05;
      minPFoil = 0.20;
   else
      minMargin = 1.0;
      maxDInferred = 0.30;
      minPFoil = 0.35;
   end

   if foilMargin < minMargin || dInferredFoil > maxDInferred || pBestFoil < minPFoil || ~visualMatch
      continue;
   end

   score = pBestFoil * foilMargin / (dInferredFoil + 0.05);
   rows(end+1) = struct( ...
      'trial', t, ...
      'participant', dm.participant(t), ...
      'setSize', dm.setSize(t), ...
      'pSwap', pSwap, ...
      'pTarget', pTarget, ...
      'pBestFoil', pBestFoil, ...
      'bestXi', bestXi, ...
      'bestFoilXi', bestFoilXi, ...
      'visualFoil', visualFoil, ...
      'visualMatch', visualMatch, ...
      'resp', resp, ...
      'target', targ, ...
      'bestFoil', numberedFoils(visualFoil), ...
      'bestFoilIdx', foilOrder(visualFoil), ...
      'dTarget', dT, ...
      'dBestFoil', dBestFoil, ...
      'dInferredFoil', dInferredFoil, ...
      'foilMargin', foilMargin, ...
      'score', score); %#ok<AGROW>
end

T = struct2table(rows);
T = sortrows(T, 'score', 'descend');

nPick = 6;
T3 = T(T.setSize == 3, :);
T6 = T(T.setSize == 6, :);
pick3 = pickDiverseSwapExamples(T3, nPick);
pick6 = pickSixItemSwapExamples(T6, nPick);
picked = [pick3; pick6];
picked = sortrows(picked, {'setSize', 'score'}, {'ascend', 'descend'});
end

function sub = pickDiverseSwapExamples(T, nPick)
if height(T) < nPick
   sub = T;
   return;
end
sub = T(1, :);
usedP = sub.participant(1);
rem = T(2:end, :);
while height(sub) < nPick && ~isempty(rem)
   best = 1;
   bestScore = -inf;
   for i = 1:height(rem)
      p = rem.participant(i);
      bonus = 0.15 * ~ismember(p, usedP);
      s = rem.score(i) * (1 + bonus);
      if s > bestScore
         bestScore = s;
         best = i;
      end
   end
   sub = [sub; rem(best, :)]; %#ok<AGROW>
   usedP(end+1) = rem.participant(best);
   rem(best, :) = [];
end
end

function pick6 = pickSixItemSwapExamples(T6, nPick)
Tother = T6(T6.bestFoilXi ~= 6, :);
Tother = sortrows(Tother, {'dInferredFoil', 'foilMargin'}, {'ascend', 'descend'});
Tfive = sortrows(T6(T6.bestFoilXi == 6, :), 'score', 'descend');

if isempty(Tother)
   pick6 = pickDiverseSwapExamples(Tfive, nPick);
   return;
end

pickOther = pickDiverseSwapExamples(Tother, 1);
usedP = pickOther.participant;
TfiveRem = Tfive(~ismember(Tfive.participant, usedP), :);
pickFive = pickDiverseSwapExamples(TfiveRem, nPick - 1);
pick6 = [pickOther; pickFive];
end

function fig = plotSwapExamples(picked, dm, pantone, piVal)
n = height(picked);
if n < 1
   error('plotSwapExamples:EmptyTable', 'picked must contain at least one trial.');
end
nCols = 4;
nRows = ceil(n / nCols);
fig = figure('Color', 'w', 'Position', [50 50 820 560], 'Units', 'pixels');
tlo = tiledlayout(nRows, nCols, 'TileSpacing', 'none', 'Padding', 'none');

xLim = [-1.12 1.12];
yLim = [-0.02 1.06];
rayLen = 0.94;
labelFs = 11;
legendFs = 12;
legendHandles = gobjects(3, 1);

for k = 1:n
   tr = picked.trial(k);
   ax = nexttile(tlo);
   hold(ax, 'on');

   targ = mod(dm.target(tr), piVal);
   nFoils = dm.setSize(tr) - 1;
   foilsRaw = mod(dm.nontarget(tr, 1:nFoils), piVal);
   resp = mod(dm.response(tr), piVal);
   foilOrder = sortFoilsByTargetDistance(targ, foilsRaw, piVal);
   inferredFoil = picked.bestFoilXi(k) - 1;

   setupSemicircleAxes(ax, xLim, yLim, pantone);
   drawSemicircleRay(ax, targ, pantone.ClassicBlue, '-', 2.8, rayLen, 0.45);
   drawSemicircleRay(ax, resp, pantone.AuroraRed, '-', 2.8, rayLen, 0.45);
   hFoil = gobjects(1);
   for fIdx = 1:nFoils
      ang = foilsRaw(foilOrder(fIdx));
      isInferred = fIdx == inferredFoil;
      h = drawSemicircleRay(ax, ang, [0 0 0], '--', 1.4, rayLen * 0.92);
      if fIdx == 1
         hFoil = h;
      end
      labelFoilNumber(ax, ang, fIdx, rayLen, isInferred, labelFs);
   end

   if k == 1
      legendHandles(1) = plot(ax, NaN, NaN, '-', 'Color', pantone.ClassicBlue, 'LineWidth', 2.8);
      legendHandles(2) = hFoil;
      legendHandles(3) = plot(ax, NaN, NaN, '-', 'Color', pantone.AuroraRed, 'LineWidth', 2.8);
   end
end

for k = (n + 1):(nRows * nCols)
   nexttile(tlo);
   axis off;
end

tlo.InnerPosition(1) = 0.02;
tlo.InnerPosition(3) = 0.96;
tlo.InnerPosition(2) = 0.04;
tlo.InnerPosition(4) = 0.90;

lgd = legend(legendHandles, {'Target', 'Foil', 'Response'}, ...
   'Orientation', 'horizontal', 'Box', 'off', 'FontSize', legendFs);
lgd.Units = 'normalized';
lgd.Position = [0.24 0.945 0.52 0.045];
end

function foilOrder = sortFoilsByTargetDistance(targ, foils, piVal)
dists = arrayfun(@(f) min(abs(targ - f), piVal - abs(targ - f)), foils);
[~, foilOrder] = sort(dists, 'ascend');
end

function labelFoilNumber(ax, angle, foilNum, rayLen, circleInferred, fontSize)
labelR = rayLen * 1.07;
x = labelR * cos(angle);
y = labelR * sin(angle);
text(ax, x, y, sprintf('%d', foilNum), ...
   'HorizontalAlignment', 'center', ...
   'VerticalAlignment', 'middle', ...
   'FontSize', fontSize, ...
   'Color', [0 0 0], ...
   'Clipping', 'off');
if circleInferred
   r = 0.050 + 0.0035 * fontSize;
   pos = [x - r, y - r, 2 * r, 2 * r];
   rectangle(ax, 'Position', pos, ...
      'Curvature', [1 1], ...
      'EdgeColor', [0 0 0], ...
      'LineWidth', 1.1, ...
      'FaceColor', 'none');
end
end

function setupSemicircleAxes(ax, xLim, yLim, pantone)
set(ax, ...
   'xlim', xLim, ...
   'ylim', yLim, ...
   'xtick', [], ...
   'ytick', [], ...
   'box', 'off', ...
   'clipping', 'off', ...
   'layer', 'top');
axis(ax, 'equal');
axis(ax, 'off');

th = linspace(0, pi, 120);
plot(ax, cos(th), sin(th), '-', 'Color', pantone.GlacierGray, 'LineWidth', 0.9);
plot(ax, xLim, [0 0], '-', 'Color', pantone.GlacierGray, 'LineWidth', 0.9);
plot(ax, 0, 0, 'o', ...
   'markerfacecolor', pantone.GlacierGray, ...
   'markeredgecolor', 'w', ...
   'markersize', 2.5);
end

function h = drawSemicircleRay(ax, angle, col, ls, lw, len, alpha)
if nargin < 7
   alpha = 1;
end
angle = mod(angle, pi);
if numel(col) == 3
   col = [col, alpha];
end
h = plot(ax, [0 len * cos(angle)], [0 len * sin(angle)], ...
   'Color', col, 'LineStyle', ls, 'LineWidth', lw);
end

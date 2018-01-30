krWinString = '';%'_krWin100';

behav_pattern = [0, 0, 1, 1; 0, 1, 0, 1; ...
    0, 1, 1, 0; 0, 1, 1, 1; 1, 0, 0, 1; ...
    1, 0, 1, 0; 1, 0, 1, 1; 1, 1, 0, 0; ...
    1, 1, 0, 1; 1, 1, 1, 0; 1, 1, 1, 1];

computationToPlot = 'FirstCorrect-R4';

if size(behav_pattern, 1) > 1
    if size(behav_pattern, 1) > 3
        behav_str = 'all';
    else
        behav_str = num2str(reshape(behav_pattern', 1, []));
    end
else
    behav_str = num2str(behav_pattern);
end
behav_str(isspace(behav_str)) = [];

repSubjects = {'JJJ', 'III', 'KKK', 'BBB', 'GGG', 'HHH', 'AAA', 'CCC', 'DDD', 'EEE', 'FFF', 'MM', ...
    'OO', 'PP', 'QQ', 'RR', 'SS', 'TT', 'WW',...
    'YY'};

origSubjects = {'AA', 'BB', 'DD', 'F', 'EE', 'GG', 'HH', 'JJ', ...
    'K', 'M', 'O', 'R', 'S', 'T', 'U', 'V', 'X', 'Y', 'Z', 'N'};
 
numRS = length(repSubjects);
numOS = length(origSubjects);

dataRootR = '/Users/nrafidi/Documents/MATLAB/compEEG-data-rep/';
dataRootO = '/Users/nrafidi/Documents/MATLAB/compEEG-data/';


fnameString = '%sresults/%s/KRanalysis_TGM_Vis%s.mat';


%%
scoreMat = [];
perfMat = [];
firstRoundCorr = [];
IndividualSubjectFirstCorr = cell(numOS + numRS, 1);
for i = 1:numRS
    load(sprintf(fnameString, dataRootR, repSubjects{i}, krWinString));
    
    goodItems = false(size(responseTraj, 1), 1);
    for i_patt = 1:size(behav_pattern, 1)
        goodItems = goodItems | ismember(responseTraj, behav_pattern(i_patt, :), 'rows');
    end
    
    for j = 1:size(responseTraj, 1)
        if goodItems(j)
            firstRoundCorr = cat(1, firstRoundCorr, find(responseTraj(j,:), 1, 'first'));
            IndividualSubjectFirstCorr{i} = cat(1, IndividualSubjectFirstCorr{i}, find(responseTraj(j,:), 1, 'first'));
        end
    end
    
    perfMat = cat(1, perfMat, mean(krLabels(goodItems), 1));
    krTrajToAdd = krTraj(goodItems, :, :, :);
    [num_samp, ~, numK, numC] = size(krTrajToAdd);
    
    krScoreToAdd = nan(num_samp, numK, numC);
    
    for iSamp = 1:num_samp
        krScoreToAdd(iSamp,:,:) = squeeze(krTrajToAdd(iSamp,IndividualSubjectFirstCorr{i}(iSamp),:,:) - ...
            krTrajToAdd(iSamp,4,:,:));
    end
    
    krScore = mean(krScoreToAdd, 1); 
    
    scoreMat = cat(1, scoreMat, krScore);
end

for i = 1:numOS
    load(sprintf(fnameString, dataRootO, origSubjects{i}, krWinString));
    
    goodItems = false(size(responseTraj, 1), 1);
    for i_patt = 1:size(behav_pattern, 1)
        goodItems = goodItems | ismember(responseTraj, behav_pattern(i_patt, :), 'rows');
    end
    
    for j = 1:size(responseTraj, 1)
        if goodItems(j)
            firstRoundCorr = cat(1, firstRoundCorr, find(responseTraj(j,:), 1, 'first'));
            IndividualSubjectFirstCorr{i+numRS} = cat(1, IndividualSubjectFirstCorr{i+numRS}, find(responseTraj(j,:), 1, 'first'));
        end
    end
    
    perfMat = cat(1, perfMat, mean(krLabels(goodItems), 1));
    krTrajToAdd = krTraj(goodItems,:,:,:);
    [num_samp, ~, numK, numC] = size(krTrajToAdd);
    
    krScoreToAdd = nan(num_samp, numK, numC);
    
    for iSamp = 1:num_samp
        krScoreToAdd(iSamp,:,:) = squeeze(krTrajToAdd(iSamp,IndividualSubjectFirstCorr{i+numRS}(iSamp),:,:) - ...
            krTrajToAdd(iSamp,4,:,:));
    end
    
    krScore = mean(krScoreToAdd, 1); 
    
    scoreMat = cat(1, scoreMat, krScore);
end

ctWinToUse = 9:38;
ctWinTime = -100:20:950;
ctWinTime = ctWinTime(ctWinToUse);
krWinToUse = 6:54;
krWinTime = winTime(krWinToUse);
if isempty(krWinString)
    scoreMat = scoreMat(:,krWinToUse, ctWinToUse);
end

%% Plot desired value
[~, numKT, numCT] = size(scoreMat);
score_corr = nan(numKT, numCT);
score_p_R = nan(numKT, numCT);
score_p_L = nan(numKT, numCT);
for kt = 1:numKT
    for ct = 1:numCT
        [score_corr(kt, ct), score_p_R(kt, ct)] = corr(perfMat, scoreMat(:, kt, ct), 'tail', 'right');
        [~, score_p_L(kt, ct)] = corr(perfMat, scoreMat(:, kt, ct), 'tail', 'left');
    end
end

pvals = zeros(numKT, numCT);
pvals(score_p_L < 0.05) = -1;
pvals(score_p_L < 0.01) = -2;
pvals(score_p_L < 0.001) = -5;

pvals(score_p_R < 0.05) = 1;
pvals(score_p_R < 0.01) = 2;
pvals(score_p_R < 0.001) = 5;

%%
f1 = figure;
imagesc(krWinTime, ctWinTime, score_corr');%, [0, trueMax]);
ylabel('Competition Training Time (ms)');
xlabel('KR Test Time (ms)');
colorbar
title(sprintf('Correlation Between Average Competition Drop\nand Subject Performance'));
set(gca, 'FontSize', 18);
set(f1, 'Color', 'w');
export_fig(f1, sprintf('%s/results/figures/KR_TGM_%s_v_abs_pooled%s_v2_behav%s_subCorr.png', dataRootR, computationToPlot, krWinString, behav_str));
export_fig(f1, sprintf('%s/results/figures/KR_TGM_%s_v_abs_pooled%s_v2_behav%s_subCorr.fig', dataRootR, computationToPlot, krWinString, behav_str));
export_fig(f1, sprintf('%s/results/figures/KR_TGM_%s_v_abs_pooled%s_v2_behav%s_subCorr.pdf', dataRootR, computationToPlot, krWinString, behav_str));

f2 = figure;
imagesc(krWinTime, ctWinTime, pvals', [-5, 5]);
ylabel('Competition Training Time (ms)');
xlabel('KR Test Time (ms)');
colorbar
title(sprintf('P value of Correlation Between\nAverage Competition Drop and Subject Performance'));
set(gca, 'FontSize', 18);
set(f2, 'Color', 'w');
% h = suptitle(computationToPlot);
% set(h, 'FontSize', 22, 'FontWeight', 'bold');

% set (f1, 'Units', 'normalized', 'Position', [0,0,1,1]);
export_fig(f2, sprintf('%s/results/figures/KR_TGM_%s_v_PVal_pooled%s_v2_behav%s_subCorr.png', dataRootR, computationToPlot, krWinString, behav_str));
export_fig(f2, sprintf('%s/results/figures/KR_TGM_%s_v_PVal_pooled%s_v2_behav%s_subCorr.fig', dataRootR, computationToPlot, krWinString, behav_str));
export_fig(f2, sprintf('%s/results/figures/KR_TGM_%s_v_PVal_pooled%s_v2_behav%s_subCorr.pdf', dataRootR, computationToPlot, krWinString, behav_str));

%%

cluster_fname = sprintf('%s/results/clusters_pVals_histograms_KRTGM_behav%s_subCorr.mat', dataRootR, behav_str);

pValString = '05';

if ~exist(cluster_fname, 'file')
    options.minClusterSize = 1;
    options.pValThresh = 0.05;
    options.maxPermutations = 1000;
    [clusters, pVals, permutationClusters, permutationHist, permutationSize] = clusterPermTestPooledSub_fullTime_Corr(perfMat, scoreMat, options);
    save(cluster_fname, ...
        'clusters', 'pVals', 'permutationHist', 'permutationSize', ...
        'permutationClusters');
else
    load(cluster_fname);
end

[~, uniClustInd] = unique(cellfun(@num2str, ...
    cellfun(@(x) reshape(x, 1, []), clusters, 'UniformOutput', false), ...
    'UniformOutput', false));
clusters = clusters(uniClustInd);
pVals = pVals(uniClustInd,:);


for i = 1:length(permutationClusters)
    [~, uniClustInd] = unique(cellfun(@num2str, ...
    cellfun(@(x) reshape(x, 1, []), permutationClusters{i}, 'UniformOutput', false), ...
    'UniformOutput', false));
    permutationClusters{i} = permutationClusters{i}(uniClustInd);
end
%%
[~, bestClust] = min(pVals(:, 2));
sizeBestClust = size(clusters{bestClust},1);
f = figure;
histogram(permutationSize, 'FaceAlpha', 1);
hold on
line([sizeBestClust, sizeBestClust], [0, 300]);
legend({'Permuted Cluster Sizes', 'True Max Cluster Size'});
title('Histogram of Permuted Cluster Sizes')
xlabel('Cluster size')
ylabel('Number of permuted clusters')
set(gca, 'FontSize', 14);
set(f, 'Color', 'w');
export_fig(f, sprintf('%s/results/figures/KR_TGM_%s_v_PVal-Cluster_pooled%s_v2_sizeHist_noSizeThresh_p%s_rightTail_behav%s_subCorr.pdf', dataRootR, computationToPlot, krWinString, pValString, behav_str));

%%
f3 = figure;
subplot(1,2,1)
histogram(permutationHist);
title('Summed t stat');
subplot(1,2,2)
histogram(permutationSize);
title('Number of points');
set(f3, 'Color', 'w');
export_fig(f3, sprintf('%s/results/figures/KR_TGM_%s_v_PVal-Cluster_pooled%s_v2_permHist_noSizeThresh_p%s_rightTail_behav%s_subCorr.png', dataRootR, computationToPlot, krWinString,pValString, behav_str));
export_fig(f3, sprintf('%s/results/figures/KR_TGM_%s_v_PVal-Cluster_pooled%s_v2_permHist_noSizeThresh_p%s_rightTail_behav%s_subCorr.fig', dataRootR, computationToPlot, krWinString,pValString, behav_str));
export_fig(f3, sprintf('%s/results/figures/KR_TGM_%s_v_PVal-Cluster_pooled%s_v2_permHist_noSizeThresh_p%s_rightTail_behav%s_subCorr.pdf', dataRootR, computationToPlot, krWinString,pValString, behav_str));
%%
f2 = figure;
imagesc(krWinTime, ctWinTime, pvals', [-5, 5]);
colorbar
hold on
% colors = 'kkm';
% i_color=1;
% for i = 1:length(clusters)
%     if pVals(i,2) <= 0.05
        scatter(krWinTime(clusters{bestClust}(:,1)), ctWinTime(clusters{bestClust}(:,2)), 'r');
%         i_color = i_color+1;
%     end
% %     scatter(krWinTime(clusters{i}(:,1)), ctWinTime(clusters{i}(:,2)), colors(i));
%     
% end
ylabel('Competition Training Time (ms)');
xlabel('KR Test Time (ms)');
title(sprintf('P value Difference\nRemembered vs Forgotten'));
set(gca, 'FontSize', 14);
set(f2, 'Color', 'w');

export_fig(f2, sprintf('%s/results/figures/KR_TGM_%s_v_PVal-Cluster_pooled%s_v2_noSizeThresh_p%s_rightTail_behav%s_subCorr.png', dataRootR, computationToPlot, krWinString,pValString, behav_str));
export_fig(f2, sprintf('%s/results/figures/KR_TGM_%s_v_PVal-Cluster_pooled%s_v2_noSizeThresh_p%s_rightTail_behav%s_subCorr.fig', dataRootR, computationToPlot, krWinString,pValString, behav_str));
export_fig(f2, sprintf('%s/results/figures/KR_TGM_%s_v_PVal-Cluster_pooled%s_v2_noSizeThresh_p%s_rightTail_behav%s_subCorr.pdf', dataRootR, computationToPlot, krWinString,pValString, behav_str));
save(sprintf('%s/results/clusters_pVals_KRTGM_behav%s_subCorr.mat', dataRootR, behav_str), ...
    'clusters', 'pVals', 'IndividualSubjectDataR', 'IndividualSubjectDataF', ...
    'krWinToUse', 'ctWinToUse', 'ctWinTime', 'krWinTime', 'IndividualSubjectFirstCorrF', ...
    'IndividualSubjectFirstCorrR');
%%

% bootstrapCluster_KRTGM_behav(bestClust, computationToPlot, behav_str);

krWinString = '';%'_krWin100';

behav_pattern = [0, 1, 1, 1];
% [0, 0, 1, 1; 0, 1, 0, 1; ...
%     0, 1, 1, 0; 0, 1, 1, 1; 1, 0, 0, 1; ...
%     1, 0, 1, 0; 1, 0, 1, 1; 1, 1, 0, 0; ...
%     1, 1, 0, 1; 1, 1, 1, 0; 1, 1, 1, 1];

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
diffMatR = [];
diffMatF = [];
firstRoundCorr_R = [];
firstRoundCorr_F = [];
IndividualSubjectFirstCorrR = cell(numOS + numRS, 1);
IndividualSubjectFirstCorrF = cell(numOS + numRS, 1);
for i = 1:numRS
    load(sprintf(fnameString, dataRootR, repSubjects{i}, krWinString));
    
    goodItems = false(size(responseTraj, 1), 1);
    for i_patt = 1:size(behav_pattern, 1)
        goodItems = goodItems | ismember(responseTraj, behav_pattern(i_patt, :), 'rows');
    end
    
    for j = 1:size(responseTraj, 1)
        if goodItems(j)
            if krLabels(j) == 1
                firstRoundCorr_R = cat(1, firstRoundCorr_R, find(responseTraj(j,:), 1, 'first'));
                IndividualSubjectFirstCorrR{i} = cat(1, IndividualSubjectFirstCorrR{i}, find(responseTraj(j,:), 1, 'first'));
            else
                firstRoundCorr_F = cat(1, firstRoundCorr_F, find(responseTraj(j,:), 1, 'first'));
                IndividualSubjectFirstCorrF{i} = cat(1, IndividualSubjectFirstCorrF{i}, find(responseTraj(j,:), 1, 'first'));
            end
        end
    end
    
    krTrajToAdd_R = krTraj(krLabels == 1 & goodItems,:,:,:);
    krTrajToAdd_F = krTraj(krLabels == 0 & goodItems,:,:,:);
    [num_samp_R, ~, numK, numC] = size(krTrajToAdd_R);
    num_samp_F = size(krTrajToAdd_F, 1);
    
    krScoreToAdd_R = nan(num_samp_R, numK, numC);
    krScoreToAdd_F = nan(num_samp_F, numK, numC);
    
    for iSamp = 1:num_samp_R
        krScoreToAdd_R(iSamp,:,:) = squeeze(krTrajToAdd_R(iSamp,IndividualSubjectFirstCorrR{i}(iSamp),:,:) - ...
            krTrajToAdd_R(iSamp,4,:,:));
    end
    for iSamp = 1:num_samp_F
        krScoreToAdd_F(iSamp,:,:) = squeeze(krTrajToAdd_F(iSamp,IndividualSubjectFirstCorrF{i}(iSamp),:,:) - ...
            krTrajToAdd_F(iSamp,4,:,:));
    end
    
    krScoreTotal = zscore([krScoreToAdd_R; krScoreToAdd_F], 1);
    
    
    diffMatR = cat(1, diffMatR, krScoreTotal(1:num_samp_R, :, :));
    diffMatF = cat(1, diffMatF, krScoreTotal((num_samp_R + 1):end, :, :));
end

for i = 1:numOS
    load(sprintf(fnameString, dataRootO, origSubjects{i}, krWinString));
    
    goodItems = false(size(responseTraj, 1), 1);
    for i_patt = 1:size(behav_pattern, 1)
        goodItems = goodItems | ismember(responseTraj, behav_pattern(i_patt, :), 'rows');
    end
    
    for j = 1:size(responseTraj, 1)
        if goodItems(j)
            if krLabels(j) == 1
                firstRoundCorr_R = cat(1, firstRoundCorr_R, find(responseTraj(j,:), 1, 'first'));
                IndividualSubjectFirstCorrR{i+numRS} = cat(1, IndividualSubjectFirstCorrR{i+numRS}, find(responseTraj(j,:), 1, 'first'));
            else
                firstRoundCorr_F = cat(1, firstRoundCorr_F, find(responseTraj(j,:), 1, 'first'));
                IndividualSubjectFirstCorrF{i+numRS} = cat(1, IndividualSubjectFirstCorrF{i+numRS}, find(responseTraj(j,:), 1, 'first'));
            end
        end
    end
    
    krTrajToAdd_R = krTraj(krLabels == 1 & goodItems,:,:,:);
    krTrajToAdd_F = krTraj(krLabels == 0 & goodItems,:,:,:);
    [num_samp_R, ~, numK, numC] = size(krTrajToAdd_R);
    num_samp_F = size(krTrajToAdd_F, 1);
    
    krScoreToAdd_R = nan(num_samp_R, numK, numC);
    krScoreToAdd_F = nan(num_samp_F, numK, numC);
    
    for iSamp = 1:num_samp_R
        krScoreToAdd_R(iSamp,:,:) = squeeze(krTrajToAdd_R(iSamp,IndividualSubjectFirstCorrR{i+numRS}(iSamp),:,:) - ...
            krTrajToAdd_R(iSamp,4,:,:));
    end
    for iSamp = 1:num_samp_F
        krScoreToAdd_F(iSamp,:,:) = squeeze(krTrajToAdd_F(iSamp,IndividualSubjectFirstCorrF{i+numRS}(iSamp),:,:) - ...
            krTrajToAdd_F(iSamp,4,:,:));
    end
    
    krScoreTotal = zscore([krScoreToAdd_R; krScoreToAdd_F], 1);
    
    
    diffMatR = cat(1, diffMatR, krScoreTotal(1:num_samp_R, :, :));
    diffMatF = cat(1, diffMatF, krScoreTotal((num_samp_R + 1):end, :, :));
end

% [~, numRounds, numKRT, numCT] = size(krTGM_R);

% % Checking that everything is gaussian
% ksTestOutput = nan(2, numKRT, numCT, 2);
%
%
% for t1 = 1:numKRT
%     for t2 = 1:numCT
%         ksTestOutput(1,t1,t2, 1) = kstest(squeeze(krTGM_R(:,1:2,t1,t2)));
%         ksTestOutput(1,t1,t2, 2) = kstest(squeeze(krTGM_F(:,1:2,t1,t2)));
%
%         ksTestOutput(2,t1,t2, 1) = kstest(squeeze(krTGM_R(:,3:4,t1,t2)));
%         ksTestOutput(2,t1,t2, 2) = kstest(squeeze(krTGM_F(:,3:4,t1,t2)));
%     end
% end

ctWinToUse = 9:38;
ctWinTime = -100:20:950;
ctWinTime = ctWinTime(ctWinToUse);
krWinToUse = 6:54;
krWinTime = winTime(krWinToUse);
if isempty(krWinString)
    diffMatF = diffMatF(:,krWinToUse, ctWinToUse);
    diffMatR = diffMatR(:,krWinToUse, ctWinToUse);
end
numKRT = length(krWinToUse);

%% Plot desired value

numR = size(diffMatR, 1);
numF = size(diffMatF, 1);
numI = min([numR, numF]);
diffMatR = diffMatR(1:numI, :, :);
diffMatF = diffMatF(1:numI, :, :);


[~, rep_corr] = ttest(diffMatR, diffMatF, 'Tail', 'right');
rep_corr = squeeze(rep_corr);

[~, rep_opp] = ttest(diffMatR, diffMatF, 'Tail', 'left');
rep_opp = squeeze(rep_opp);

rep_diff = squeeze(mean(diffMatR - diffMatF, 1));

rep = zeros(size(rep_opp));
rep(rep_opp < 0.05) = -1;
rep(rep_opp < 0.01) = -2;
rep(rep_opp < 0.001) = -5;

rep(rep_corr < 0.05) = 1;
rep(rep_corr < 0.01) = 2;
rep(rep_corr < 0.001) = 5;

figure
meow = rep_corr';
imagesc(krWinTime, ctWinTime, meow, [0, 0.05])

[minval, minind] = min(meow(:));
[row, column] = ind2sub(size(meow), minind);
disp(row)
disp(column)
disp(ctWinTime(row))
disp(krWinTime(column))


%%
f1 = figure;
imagesc(krWinTime, ctWinTime, rep_diff');%, [0, trueMax]);
ylabel('Competition Training Time (ms)');
xlabel('KR Test Time (ms)');
colorbar
title(sprintf('Absolute Difference\nRemembered vs Forgotten'));
set(gca, 'FontSize', 18);
set(f1, 'Color', 'w');
export_fig(f1, sprintf('%s/results/figures/KR_TGM_%s_v_abs_pooled%s_v2_behav%s_Z.png', dataRootR, computationToPlot, krWinString, behav_str));
export_fig(f1, sprintf('%s/results/figures/KR_TGM_%s_v_abs_pooled%s_v2_behav%s_Z.fig', dataRootR, computationToPlot, krWinString, behav_str));
export_fig(f1, sprintf('%s/results/figures/KR_TGM_%s_v_abs_pooled%s_v2_behav%s_Z.pdf', dataRootR, computationToPlot, krWinString, behav_str));

f2 = figure;
imagesc(krWinTime, ctWinTime, rep', [-5, 5]);
ylabel('Competition Training Time (ms)');
xlabel('KR Test Time (ms)');
colorbar
title(sprintf('P value Difference\n%s Trials', behav_str));
set(gca, 'FontSize', 18);
set(f2, 'Color', 'w');
% h = suptitle(computationToPlot);
% set(h, 'FontSize', 22, 'FontWeight', 'bold');

% set (f1, 'Units', 'normalized', 'Position', [0,0,1,1]);
export_fig(f2, sprintf('%s/results/figures/KR_TGM_%s_v_PVal_pooled%s_v2_behav%s_Z.png', dataRootR, computationToPlot, krWinString, behav_str));
export_fig(f2, sprintf('%s/results/figures/KR_TGM_%s_v_PVal_pooled%s_v2_behav%s_Z.fig', dataRootR, computationToPlot, krWinString, behav_str));
export_fig(f2, sprintf('%s/results/figures/KR_TGM_%s_v_PVal_pooled%s_v2_behav%s_Z.pdf', dataRootR, computationToPlot, krWinString, behav_str));
keyboard
%%

cluster_fname = sprintf('%s/results/clusters_pVals_histograms_KRTGM_behav%s_Z.mat', dataRootR, behav_str);

dataToCluster = cat(4, diffMatR, diffMatF);
pValString = '05';

if ~exist(cluster_fname, 'file')
    options.minClusterSize = 1;
    options.pValThresh = 0.05;
    options.maxPermutations = 1000;
    [clusters, pVals, permutationClusters, permutationHist, permutationSize] = clusterPermTestPooledSub_fullTime(dataToCluster, options);
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
export_fig(f, sprintf('%s/results/figures/KR_TGM_%s_v_PVal-Cluster_pooled%s_v2_sizeHist_noSizeThresh_p%s_rightTail_behav%s_Z.pdf', dataRootR, computationToPlot, krWinString, pValString, behav_str));

%%
f3 = figure;
subplot(1,2,1)
histogram(permutationHist);
title('Summed t stat');
subplot(1,2,2)
histogram(permutationSize);
title('Number of points');
set(f3, 'Color', 'w');
export_fig(f3, sprintf('%s/results/figures/KR_TGM_%s_v_PVal-Cluster_pooled%s_v2_permHist_noSizeThresh_p%s_rightTail_behav%s_Z.png', dataRootR, computationToPlot, krWinString,pValString, behav_str));
export_fig(f3, sprintf('%s/results/figures/KR_TGM_%s_v_PVal-Cluster_pooled%s_v2_permHist_noSizeThresh_p%s_rightTail_behav%s_Z.fig', dataRootR, computationToPlot, krWinString,pValString, behav_str));
export_fig(f3, sprintf('%s/results/figures/KR_TGM_%s_v_PVal-Cluster_pooled%s_v2_permHist_noSizeThresh_p%s_rightTail_behav%s_Z.pdf', dataRootR, computationToPlot, krWinString,pValString, behav_str));
%%
f2 = figure;
imagesc(krWinTime, ctWinTime, rep', [-5, 5]);
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

export_fig(f2, sprintf('%s/results/figures/KR_TGM_%s_v_PVal-Cluster_pooled%s_v2_noSizeThresh_p%s_rightTail_behav%s_Z.png', dataRootR, computationToPlot, krWinString,pValString, behav_str));
export_fig(f2, sprintf('%s/results/figures/KR_TGM_%s_v_PVal-Cluster_pooled%s_v2_noSizeThresh_p%s_rightTail_behav%s_Z.fig', dataRootR, computationToPlot, krWinString,pValString, behav_str));
export_fig(f2, sprintf('%s/results/figures/KR_TGM_%s_v_PVal-Cluster_pooled%s_v2_noSizeThresh_p%s_rightTail_behav%s_Z.pdf', dataRootR, computationToPlot, krWinString,pValString, behav_str));
save(sprintf('%s/results/clusters_pVals_KRTGM_behav%s_Z.mat', dataRootR, behav_str), ...
    'clusters', 'pVals', 'IndividualSubjectDataR', 'IndividualSubjectDataF', ...
    'krWinToUse', 'ctWinToUse', 'ctWinTime', 'krWinTime', 'IndividualSubjectFirstCorrF', ...
    'IndividualSubjectFirstCorrR');
%%

bootstrapCluster_KRTGM_behav(bestClust, computationToPlot, behav_str);

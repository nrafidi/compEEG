krWinString = '';%'_krWin100';

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
krTGM_R = [];
krTGM_F = [];
firstRoundCorr_R = [];
firstRoundCorr_F = [];
IndividualSubjectDataR = cell(numOS + numRS, 1);
IndividualSubjectDataF = cell(numOS + numRS, 1);
IndividualSubjectFirstCorrR = cell(numOS + numRS, 1);
IndividualSubjectFirstCorrF = cell(numOS + numRS, 1);
for i = 1:numRS
    load(sprintf(fnameString, dataRootR, repSubjects{i}, krWinString));
    
    goodItems = sum(responseTraj, 2) > 1 | sum(responseTraj(:,1:3),2) > 0;
    
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
    
    IndividualSubjectDataR{i} = krTrajToAdd_R;
    IndividualSubjectDataF{i} = krTrajToAdd_F;
%     numSamp = min([sum(krLabels == 1 & goodItems), sum(krLabels == 0 & goodItems)]);
    
    krTGM_R = cat(1, krTGM_R, krTrajToAdd_R);%(1:numSamp,:,:,:));
    krTGM_F = cat(1, krTGM_F, krTrajToAdd_F);%(1:numSamp,:,:,:));
end

for i = 1:numOS
    load(sprintf(fnameString, dataRootO, origSubjects{i}, krWinString));
    
    goodItems = sum(responseTraj, 2) > 1 | sum(responseTraj(:,1:3),2) > 0;
    
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
%     numSamp = min([sum(krLabels == 1 & goodItems), sum(krLabels == 0 & goodItems)]);
    
    IndividualSubjectDataR{i+numRS} = krTrajToAdd_R;
    IndividualSubjectDataF{i+numRS} = krTrajToAdd_F;

    krTGM_R = cat(1, krTGM_R, krTrajToAdd_R);%(1:numSamp,:,:,:));
    krTGM_F = cat(1, krTGM_F, krTrajToAdd_F);%(1:numSamp,:,:,:));
end

[~, numRounds, numKRT, numCT] = size(krTGM_R);

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
krWinToUse = 6:54;%6:40;
krWinTime = winTime(krWinToUse);
if isempty(krWinString)
    krTGM_F = krTGM_F(:,:,krWinToUse, ctWinToUse);
    krTGM_R = krTGM_R(:,:,krWinToUse, ctWinToUse);
end
numKRT = length(krWinToUse);

%% Plot desired value

computationToPlot = 'FirstCorrect-R4';

numF = size(krTGM_F, 1);
numR = size(krTGM_R, 1);
numCT_short = size(krTGM_F, 4);
numSamp = min([numF, numR]);

switch computationToPlot
    case 'Slope'
        
        diffMatR = nan(numSamp, numKRT, numCT_short);
        diffMatF = nan(numSamp, numKRT, numCT_short);
        
        for i = 1:numSamp
            for j = 1:numKRT
                for k = 1:numCT_short
                    weights = polyfit(1:4, squeeze(krTGM_R(i,:,j,k)),1);
                    diffMatR(i,j,k) = weights(1);
                    weights = polyfit(1:4, squeeze(krTGM_F(i,:,j,k)),1);
                    diffMatF(i,j,k) = weights(1);
                end
            end
        end
        
    case 'Mean(R1R2-R3R4)'
        diffMatR = squeeze(mean(krTGM_R(1:numSamp,1:2,:,:), 2) - mean(krTGM_R(1:numSamp,3:4,:,:), 2));
        diffMatF = squeeze(mean(krTGM_F(1:numSamp,1:2,:,:), 2) - mean(krTGM_F(1:numSamp,3:4,:,:), 2));
    case 'Mean(R3R4)'
        diffMatR = squeeze(mean(krTGM_R(1:numSamp,3:4,:,:), 2));
        diffMatF = squeeze(mean(krTGM_F(1:numSamp,3:4,:,:), 2));
    case 'R2-R3'
        diffMatR = squeeze(krTGM_R(1:numSamp,2,:,:) - krTGM_R(1:numSamp,3,:,:));
        diffMatF = squeeze(krTGM_F(1:numSamp,2,:,:) - krTGM_F(1:numSamp,3,:,:));
    case 'R3-R4'
        diffMatR = squeeze(krTGM_R(1:numSamp,3,:,:) - krTGM_R(1:numSamp,4,:,:));
        diffMatF = squeeze(krTGM_F(1:numSamp,3,:,:) - krTGM_F(1:numSamp,4,:,:));
    case 'R1-R4'
        diffMatR = squeeze(krTGM_R(1:numSamp,1,:,:) - krTGM_R(1:numSamp,4,:,:));
        diffMatF = squeeze(krTGM_F(1:numSamp,1,:,:) - krTGM_F(1:numSamp,4,:,:));
    case 'R2-R4'
        diffMatR = squeeze(krTGM_R(1:numSamp,2,:,:) - krTGM_R(1:numSamp,4,:,:));
        diffMatF = squeeze(krTGM_F(1:numSamp,2,:,:) - krTGM_F(1:numSamp,4,:,:));
    case 'R2-Mean(R3R4)'
        diffMatR = squeeze(krTGM_R(1:numSamp,2,:,:) - mean(krTGM_R(1:numSamp,3:4,:,:), 2));
        diffMatF = squeeze(krTGM_F(1:numSamp,2,:,:) - mean(krTGM_F(1:numSamp,3:4,:,:), 2));
    case 'R1-Mean(R2R3R4)'
        diffMatR = squeeze(krTGM_R(1:numSamp,1,:,:) - mean(krTGM_R(1:numSamp,2:4,:,:), 2));
        diffMatF = squeeze(krTGM_F(1:numSamp,1,:,:) - mean(krTGM_F(1:numSamp,2:4,:,:), 2));
    case 'FirstCorrect-R4'
        diffMatR = nan(numSamp, numKRT, numCT_short);
        diffMatF = nan(numSamp, numKRT, numCT_short);
        for i = 1:numSamp
            diffMatR(i,:,:) = squeeze(krTGM_R(i,firstRoundCorr_R(i),:,:) - krTGM_R(i,4,:,:));
            diffMatF(i,:,:) = squeeze(krTGM_F(i,firstRoundCorr_F(i),:,:) - krTGM_F(i,4,:,:));
        end
    otherwise
        error('Not Implemented')
end

%%

diffMatAll = cat(1, diffMatR, diffMatF);

[~, rep_corr] = ttest(diffMatAll, 0, 'Tail', 'right');
rep_corr = squeeze(rep_corr);

[~, rep_opp] = ttest(diffMatAll, 0, 'Tail', 'left');
rep_opp = squeeze(rep_opp);

rep_diff = squeeze(mean(diffMatAll, 1));

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

disp(ctWinTime(row))
disp(krWinTime(column))
%%
f1 = figure;
subplot(1, 2, 1)
imagesc(krWinTime, ctWinTime, rep_diff');%, [0, trueMax]);
ylabel('Competition Training Time (ms)');
xlabel('KR Test Time (ms)');
colorbar
title(sprintf('Absolute Difference\nRemembered vs Forgotten'));
set(gca, 'FontSize', 18);

subplot(1, 2, 2)
imagesc(krWinTime, ctWinTime, rep', [-5, 5]);
colorbar
title(sprintf('P value Difference\nRemembered vs Forgotten'));
set(gca, 'FontSize', 18);
% h = suptitle(computationToPlot);
% set(h, 'FontSize', 22, 'FontWeight', 'bold');
set(f1, 'Color', 'w');
set (f1, 'Units', 'normalized', 'Position', [0,0,1,1]);
export_fig(f1, sprintf('%s/results/figures/KR_TGM_%s_v_PVal_pooled%s_v2_droponly.png', dataRootR, computationToPlot, krWinString));
export_fig(f1, sprintf('%s/results/figures/KR_TGM_%s_v_PVal_pooled%s_v2_droponly.fig', dataRootR, computationToPlot, krWinString));
export_fig(f1, sprintf('%s/results/figures/KR_TGM_%s_v_PVal_pooled%s_v2_droponly.pdf', dataRootR, computationToPlot, krWinString));

%%
% dataToCluster = cat(4, diffMatAll, diffMatF);
% pValString = '05';
% 
% clusterFname = sprintf('%s/results/clusters_pVals_histograms_KRTGM_droponly.mat', dataRootR);
% 
% if ~exist(clusterFname, 'file')
%     options.minClusterSize = 1;
%     options.pValThresh = 0.05;
%     options.maxPermutations = 1000;
%     [clusters, pVals, permutationClusters, permutationHist, permutationSize] = clusterPermTestPooledSub_fullTime(dataToCluster, options);
%     save(sprintf('%s/results/clusters_pVals_histograms_KRTGM.mat', dataRootR), ...
%         'clusters', 'pVals', 'permutationHist', 'permutationSize', ...
%         'permutationClusters');
% else
%     load(clusterFname);
% end
% 
% [~, uniClustInd] = unique(cellfun(@num2str, ...
%     cellfun(@(x) reshape(x, 1, []), clusters, 'UniformOutput', false), ...
%     'UniformOutput', false));
% clusters = clusters(uniClustInd);
% pVals = pVals(uniClustInd,:);
% 
% 
% for i = 1:length(permutationClusters)
%     [~, uniClustInd] = unique(cellfun(@num2str, ...
%     cellfun(@(x) reshape(x, 1, []), permutationClusters{i}, 'UniformOutput', false), ...
%     'UniformOutput', false));
%     permutationClusters{i} = permutationClusters{i}(uniClustInd);
% end
% 
% sigClust = find(pVals(:,2) < 0.05);
% sizeSigClust = size(clusters{sigClust},1);
% f = figure;
% histogram(permutationSize, 'FaceAlpha', 1);
% hold on
% line([sizeSigClust, sizeSigClust], [0, 300]);
% legend({'Permuted Cluster Sizes', 'True Max Cluster Size'});
% title('Histogram of Permuted Cluster Sizes')
% xlabel('Cluster size')
% ylabel('Number of permuted clusters')
% set(gca, 'FontSize', 14);
% set(f, 'Color', 'w');
% export_fig(f, sprintf('%s/results/figures/KR_TGM_%s_v_PVal-Cluster_pooled%s_v2_sizeHist_noSizeThresh_p%s_rightTail.pdf', dataRootR, computationToPlot, krWinString,pValString));
% 
% %%
% f3 = figure;
% subplot(1,2,1)
% histogram(permutationHist);
% title('Summed t stat');
% subplot(1,2,2)
% histogram(permutationSize);
% title('Number of points');
% set(f3, 'Color', 'w');
% export_fig(f3, sprintf('%s/results/figures/KR_TGM_%s_v_PVal-Cluster_pooled%s_v2_permHist_noSizeThresh_p%s_rightTail.png', dataRootR, computationToPlot, krWinString,pValString));
% export_fig(f3, sprintf('%s/results/figures/KR_TGM_%s_v_PVal-Cluster_pooled%s_v2_permHist_noSizeThresh_p%s_rightTail.fig', dataRootR, computationToPlot, krWinString,pValString));
% export_fig(f3, sprintf('%s/results/figures/KR_TGM_%s_v_PVal-Cluster_pooled%s_v2_permHist_noSizeThresh_p%s_rightTail.pdf', dataRootR, computationToPlot, krWinString,pValString));
% %%
% f2 = figure;
% imagesc(krWinTime, ctWinTime, rep', [-5, 5]);
% colorbar
% hold on
% % colors = 'abcdrbkm';
% for i = 1:length(clusters)
%     if pVals(i,2) <= 0.01
%         scatter(krWinTime(clusters{i}(:,1)), ctWinTime(clusters{i}(:,2)), 'r');
%     end
% %     scatter(krWinTime(clusters{i}(:,1)), ctWinTime(clusters{i}(:,2)), colors(i));
%     
% end
% ylabel('Competition Training Time (ms)');
% xlabel('KR Test Time (ms)');
% title(sprintf('P value Difference\nRemembered vs Forgotten'));
% set(gca, 'FontSize', 14);
% set(f2, 'Color', 'w');
% 
% export_fig(f2, sprintf('%s/results/figures/KR_TGM_%s_v_PVal-Cluster_pooled%s_v2_noSizeThresh_p%s_rightTail.png', dataRootR, computationToPlot, krWinString,pValString));
% export_fig(f2, sprintf('%s/results/figures/KR_TGM_%s_v_PVal-Cluster_pooled%s_v2_noSizeThresh_p%s_rightTail.fig', dataRootR, computationToPlot, krWinString,pValString));
% export_fig(f2, sprintf('%s/results/figures/KR_TGM_%s_v_PVal-Cluster_pooled%s_v2_noSizeThresh_p%s_rightTail.pdf', dataRootR, computationToPlot, krWinString,pValString));
% save(sprintf('%s/results/clusters_pVals_KRTGM.mat', dataRootR), ...
%     'clusters', 'pVals', 'IndividualSubjectDataR', 'IndividualSubjectDataF', ...
%     'krWinToUse', 'ctWinToUse', 'ctWinTime', 'krWinTime', 'IndividualSubjectFirstCorrF', ...
%     'IndividualSubjectFirstCorrR');
% % numPerms = 100;
% % repForPerm = rep';
% %
% % [numRows, numCols] = size(repForPerm);
% %
% % maxBlockSize = nan(numRows, 1);
% % maxBlockSize_perm = nan(numRows, numPerms);
% % percGreater = nan(numRows, 1);
% % for i = 1:numRows
% %     blocks = findContiguousNZ(repForPerm(i,:));
% %     if ~isempty(blocks)
% %         maxBlockSize(i) = max(cellfun(@numel, blocks));
% %         for p = 1:numPerms
% %             rng(p);
% %             maxBlockSize_perm(i,p) = max(cellfun(@numel, ...
% %                 findContiguousNZ(repForPerm(i,randperm(numCols)))));
% %         end
% %         percGreater(i) = sum(maxBlockSize_perm(i,:) > maxBlockSize(i));
% %     end
% % end
% 

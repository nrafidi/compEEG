repSubjects = {'JJJ', 'III', 'KKK', 'BBB', 'GGG', 'HHH', 'AAA', 'CCC', 'DDD', 'EEE', 'FFF', 'MM', ...
    'OO', 'PP', 'QQ', 'RR', 'SS', 'TT', 'WW',...
    'YY'};

origSubjects = {'AA', 'BB', 'DD', 'F', 'EE', 'GG', 'HH', 'JJ', ...
    'K', 'M', 'O', 'R', 'S', 'T', 'U', 'V', 'X', 'Y', 'Z', 'N'};

numRS = length(repSubjects);
numOS = length(origSubjects);

dataRootR = '/Users/nrafidi/Documents/MATLAB/compEEG-data-rep/';
dataRootO = '/Users/nrafidi/Documents/MATLAB/compEEG-data/';


fnameString = '%sresults/%s/KRanalysis_TGM_Vis.mat';

repTGM_R = [];
repTGM_F = [];
origTGM_R = [];
origTGM_F = [];

for i = 1:numRS
    load(sprintf(fnameString, dataRootR, repSubjects{i}));
    repTGM_R = cat(1, repTGM_R, krTraj(krLabels == 1,:,:,:));
    repTGM_F = cat(1, repTGM_F, krTraj(krLabels == 0,:,:,:));
end

for i = 1:numOS
    load(sprintf(fnameString, dataRootO, origSubjects{i}));
    origTGM_R = cat(1, origTGM_R, krTraj(krLabels == 1,:,:,:));
    origTGM_F = cat(1, origTGM_F, krTraj(krLabels == 0,:,:,:));
end

[~, numRounds, numKRT, numCT] = size(repTGM_R);

% Checking that everything is gaussian
% ksTestOutput = nan(2, numKRT, numCT, 4);
% 
% % for r = 1:numRounds
%     for t1 = 1:numKRT
%         for t2 = 1:numCT
%             ksTestOutput(1,t1,t2, 1) = kstest(squeeze(repTGM_R(:,1:2,t1,t2)));
%             ksTestOutput(1,t1,t2, 2) = kstest(squeeze(repTGM_F(:,1:2,t1,t2)));
%             ksTestOutput(1,t1,t2, 3) = kstest(squeeze(origTGM_R(:,1:2,t1,t2)));
%             ksTestOutput(1,t1,t2, 4) = kstest(squeeze(origTGM_F(:,1:2,t1,t2)));
%             
%             ksTestOutput(2,t1,t2, 1) = kstest(squeeze(repTGM_R(:,3:4,t1,t2)));
%             ksTestOutput(2,t1,t2, 2) = kstest(squeeze(repTGM_F(:,3:4,t1,t2)));
%             ksTestOutput(2,t1,t2, 3) = kstest(squeeze(origTGM_R(:,3:4,t1,t2)));
%             ksTestOutput(2,t1,t2, 4) = kstest(squeeze(origTGM_F(:,3:4,t1,t2)));
%         end
%     end
% end




% repTGM_R = squeeze(mean(repTGM_R, 1));
% repTGM_F = squeeze(mean(repTGM_F, 1));
% origTGM_R = squeeze(mean(origTGM_R, 1));
% origTGM_F = squeeze(mean(origTGM_F, 1));



%% Plot Each Round Separately
numOrigF = size(origTGM_F, 1);
numOrigR = size(origTGM_R, 1);
numRepF = size(repTGM_F, 1);
numRepR = size(repTGM_R, 1);

numOrig = min([numOrigF, numOrigR]);
numRep = min([numRepF, numRepR]);

for r = 1:numRounds
    
    [~, orig] = ttest(origTGM_R(1:numOrig,r,:,:), ...
        origTGM_F(1:numOrig,r,:,:));
    orig = squeeze(orig);
    maxO = max(max(abs(orig)));
    [~, rep] = ttest(repTGM_R(1:numRep,r,:,:), ...
        repTGM_F(1:numRep,r,:,:));
    rep = squeeze(rep);
    maxR = max(max(abs(rep)));
    
    
    trueMax = max([maxO, maxR]);
    
    f1 = figure;
    subplot(1,2,1);
    imagesc(winTime(1:numKRT), winTime(1:numCT), orig, [0, trueMax]);
    title('Original');
    subplot(1,2,2);
    imagesc(winTime(1:numKRT), winTime(1:numCT), rep, [0, trueMax]);
    title('Replication');
    suptitle(sprintf('Round %d\nP Values', r));
%     
%     f2 = figure;
%     diffMatrix = abs(orig-rep)/trueMax;
%     imagesc(winTime(1:numKRT), winTime(1:numCT), diffMatrix);
%     title(sprintf('Difference between Original and Replication\nRound %d', r));
    
    f3 = figure;
    plot(winTime(1:numKRT), orig(:,17));
    hold on
    plot(winTime(1:numKRT), rep(:,18), 'r');
    line([winTime(1), winTime(numKRT)], [0.05, 0.05], 'Color', 'k', 'LineStyle', '--');
end
%% Plot meanR1R2

orig = squeeze(mean(origTGM_R(1:2,:,:),1) - mean(origTGM_F(1:2,:,:),1));
maxO = max(max(abs(orig)));
rep = squeeze(mean(repTGM_R(1:2,:,:),1) - mean(repTGM_F(1:2,:,:),1));
maxR = max(max(abs(rep)));

trueMax = max([maxO, maxR]);

f1 = figure;
subplot(1,2,1);
imagesc(winTime(1:numKRT), winTime(1:numCT), orig, [-trueMax, trueMax]);
title('Original');
subplot(1,2,2);
imagesc(winTime(1:numKRT), winTime(1:numCT), rep, [-trueMax, trueMax]);
title('Replication');
suptitle('Mean R1R2');

f2 = figure;
imagesc(winTime(1:numKRT), winTime(1:numCT), abs(orig-rep));
title(sprintf('Difference between Original and Replication\nMean R1R2'));

% f3 = figure;
% plot(winTime(1:numKRT), orig(:,17));
% hold on
% plot(winTime(1:numKRT), rep(:,18), 'r');

%% Plot meanR3R4
numOrigF = size(origTGM_F, 1);
numOrigR = size(origTGM_R, 1);
numRepF = size(repTGM_F, 1);
numRepR = size(repTGM_R, 1);

numOrig = min([numOrigF, numOrigR]);
numRep = min([numRepF, numRepR]);
    
% orig = squeeze(mean(origTGM_R(3:4,:,:),1) - mean(origTGM_F(3:4,:,:),1));
[~, orig] = ttest(mean(origTGM_R(1:numOrig,3:4,:,:), 2), ...
        mean(origTGM_F(1:numOrig,3:4,:,:), 2));
    orig = squeeze(orig);
orig_diff = squeeze(mean(mean(origTGM_R(1:numOrig,3:4,:,:), 2) - mean(origTGM_F(1:numOrig,3:4,:,:), 2), 1));

maxO = max(max(abs(orig)));
% rep = squeeze(mean(repTGM_R(3:4,:,:),1) - mean(repTGM_F(3:4,:,:),1));
[~, rep] = ttest(mean(repTGM_R(1:numRep,3:4,:,:), 2), ...
        mean(repTGM_F(1:numRep,3:4,:,:), 2));
    rep = squeeze(rep);

rep_diff = squeeze(mean(mean(repTGM_R(1:numRep,3:4,:,:), 2) - mean(repTGM_F(1:numRep,3:4,:,:), 2), 1));
    
maxR = max(max(abs(rep)));

trueMax = max([maxO, maxR]);

f1 = figure;
subplot(1,2,1);
imagesc(winTime(1:numKRT), winTime(1:numCT), orig, [0, trueMax]);
hold on
% line([-100, 960], [500, 500], 'Color', 'k');
% line([-100, 960], [700, 700], 'Color', 'k');
hold off
title('Original');
subplot(1,2,2);
imagesc(winTime(1:numKRT), winTime(1:numCT), rep, [0, trueMax]);
hold on
% line([-100, 960], [500, 500], 'Color', 'k');
% line([-100, 960], [700, 700], 'Color', 'k');
hold off
title('Replication');
suptitle('Mean R3R4');
export_fig(f1, sprintf('%s/results/figures/KR_TGM_OvR_R3R4_PVal.png', dataRootR));

[UR, SR, VR] = svd(rep);
[UO, SO, VO] = svd(orig);
rankToKeep = 5;
repEst = UR(:,1:rankToKeep)*SR(1:rankToKeep, 1:rankToKeep)*VR(:,1:rankToKeep)';
origEst = UO(:,1:rankToKeep)*SO(1:rankToKeep, 1:rankToKeep)*VO(:,1:rankToKeep)';

f15 = figure;
subplot(1,2,1);
imagesc(winTime(1:numKRT), winTime(1:numCT), origEst, [0, trueMax]);
hold on
hold off
title('Original');
subplot(1,2,2);
imagesc(winTime(1:numKRT), winTime(1:numCT), repEst, [0, trueMax]);
hold on
hold off
title('Replication');
suptitle('Mean R3R4');
export_fig(f1, sprintf('%s/results/figures/KR_TGM_OvR_R3R4_PVal_Est.png', dataRootR));


f2 = figure;
diffMatrix = abs(orig-rep);
imagesc(winTime(1:numKRT), winTime(1:numCT), diffMatrix);
title(sprintf('Difference between Original and Replication\nMean R3R4'));
export_fig(f2, sprintf('%s/results/figures/KR_TGM_Diff_R3R4_PVal.png', dataRootR));
% 
% [U, S, V] = svd(diffMatrix);
% lowRankDiff = U(:,1:2)*S(1:2,1:2)*V(:,1:2)';
% f25 = figure;
% imagesc(winTime(1:numKRT), winTime(1:numCT), lowRankDiff);

% f3 = figure;
% corrOR = corr(orig, rep);
% imagesc(winTime(1:numCT), winTime(1:numCT), corrOR, [0, 1]);
% xlabel('Comp Time');
% ylabel('Comp Time');
% title(sprintf('Correlation between Original and Replication\nMean R3R4'));
% export_fig(f3, sprintf('%s/results/figures/KR_TGM_Corr_R3R4_PVal.png', dataRootR));
% 
% [maxCorr, indMax] = max(corrOR(:));
% [i, j] = ind2sub([numCT, numCT], indMax);
% while abs(i-j) > 2
%     corrOR(indMax) = 0;
%     [maxCorr, indMax] = max(corrOR(:));
%     [i, j] = ind2sub([numCT, numCT], indMax);
% end

diffMatrix(orig > 0.05 | rep > 0.05) = 1e20;
diffMatrix(winTime(1:numCT) < 200, :) = 1e20;
diffMatrix(winTime(1:numCT) > 530, :) = 1e20;
[minDiff, indMin] = min(diffMatrix(:));
[i, j] = ind2sub([numCT, numKRT], indMin);
% while winTime(i) < 200 || winTime(i) > 350 || winTime(j) < 200
%     diffMatrix(indMin) = 1e20;
%     [minDiff, indMin] = min(diffMatrix(:));
%     [i, j] = ind2sub([numCT, numKRT], indMin);
% end

heatMap = orig < 0.05 & rep < 0.05;
f35 = figure;
imagesc(heatMap)
numNZ = sum(heatMap, 1);

% [~, i] = max(numNZ(6:45));
for i = 29
f4 = figure;
plot(winTime(1:numKRT), orig(:,i));
hold on
plot(winTime(1:numKRT), rep(:,i), 'r');
% line([winTime(1), winTime(numKRT)], [0.05, 0.05], 'Color', 'k', 'LineStyle', '--');
xlabel('KR Time');
ylabel(sprintf('P value Difference between Remembered and Forgotten'));
title(sprintf('R - F, Mean R3R4\nTrained on Comp Time (%0.0f,%0.0f)', winTime(i), winTime(i)));
legend({'Original', 'Replication'});
set(f4, 'Color', 'w');
export_fig(f4, sprintf('%s/results/figures/KR_TGM_MaxCorr_R3R4_PVal_CWin%d.pdf', dataRootR, i));

f5 = figure;
plot(winTime(1:numKRT), orig_diff(:,i));
hold on
plot(winTime(1:numKRT), rep_diff(:,i), 'r');
% line([winTime(1), winTime(numKRT)], [0.05, 0.05], 'Color', 'k', 'LineStyle', '--');
xlabel('KR Time');
ylabel(sprintf('Difference between Remembered and Forgotten'));
title(sprintf('R - F, Mean R3R4\nTrained on Comp Time (%0.0f,%0.0f)', winTime(i), winTime(i)));
legend({'Original', 'Replication'});
set(f5, 'Color', 'w');
export_fig(f5, sprintf('%s/results/figures/KR_TGM_MaxCorr_R3R4_PVal_CWin%d.pdf', dataRootR, i));

end

%% Plot meanR2R2 - meanR3R4

orig = squeeze((mean(origTGM_R(1:2,:,:),1)-mean(origTGM_R(3:4,:,:),1)) - ...
    (mean(origTGM_F(1:2,:,:),1) - mean(origTGM_F(3:4,:,:),1)));
maxO = max(max(abs(orig)));
rep = squeeze((mean(repTGM_R(1:2,:,:),1)-mean(repTGM_R(3:4,:,:),1)) - ...
    (mean(repTGM_F(1:2,:,:),1) - mean(repTGM_F(3:4,:,:),1)));
maxO = max(max(abs(rep)));

trueMax = max([maxO, maxR]);

f1 = figure;
subplot(1,2,1);
imagesc(winTime(1:numKRT), winTime(1:numCT), orig, [-trueMax, trueMax]);
title('Original');
subplot(1,2,2);
imagesc(winTime(1:numKRT), winTime(1:numCT), rep, [-trueMax, trueMax]);
title('Replication');
suptitle('Mean R1R2-R3R4');
export_fig(f1, sprintf('%s/results/figures/KR_TGM_OvR_MeanDiff.png', dataRootR));


f2 = figure;
diffMatrix = abs(orig-rep);
imagesc(winTime(1:numKRT), winTime(1:numCT), diffMatrix);
title(sprintf('Difference between Original and Replication\nMean R1R2-R3R4'));
export_fig(f2, sprintf('%s/results/figures/KR_TGM_Diff_MeanDiff.png', dataRootR));

[U, S, V] = svd(diffMatrix);
lowRankDiff = U(:,1:2)*S(1:2,1:2)*V(:,1:2)';
f25 = figure;
imagesc(winTime(1:numKRT), winTime(1:numCT), lowRankDiff);

f3 = figure;
corrOR = corr(orig, rep);
imagesc(winTime(1:numCT), winTime(1:numCT), corrOR, [0, 1]);
xlabel('Comp Time');
ylabel('Comp Time');
title(sprintf('Correlation between Original and Replication\nMean R1R2-R3R4'));
export_fig(f3, sprintf('%s/results/figures/KR_TGM_Corr_MeanDiff.png', dataRootR));

[maxCorr, indMax] = max(corrOR(:));
[i, j] = ind2sub([numCT, numCT], indMax);
while abs(i-j) > 3
    corrOR(indMax) = 0;
    [maxCorr, indMax] = max(corrOR(:));
    [i, j] = ind2sub([numCT, numCT], indMax);
end


f4 = figure;
plot(winTime(1:numKRT), orig(:,i));
hold on
plot(winTime(1:numKRT), rep(:,j), 'r');
% line([500, 500], [-0.08, 0.06], 'Color', 'k');
% line([700, 700], [-0.08, 0.06], 'Color', 'k');
xlabel('KR Time');
ylabel(sprintf('Difference between Remembered and Forgotten'));
title(sprintf('R - F, Mean R1R2-R3R4\nTrained on Comp Time (%0.0f,%0.0f)', winTime(i), winTime(j)));
legend({'Original', 'Replication'});
export_fig(f4, sprintf('%s/results/figures/KR_TGM_MaxCorr_MeanDiff.pdf', dataRootR));
%% Plot R2 - R3

orig = squeeze((mean(origTGM_R(2,:,:),1)-mean(origTGM_R(3,:,:),1)) - ...
    (mean(origTGM_F(2,:,:),1) - mean(origTGM_F(3,:,:),1)));
rep = squeeze((mean(repTGM_R(2,:,:),1)-mean(repTGM_R(3,:,:),1)) - ...
    (mean(repTGM_F(2,:,:),1) - mean(repTGM_F(3,:,:),1)));

f1 = figure;
subplot(1,2,1);
imagesc(winTime(1:numKRT), winTime(1:numCT), orig);
title('Original');
subplot(1,2,2);
imagesc(winTime(1:numKRT), winTime(1:numCT), rep);
title('Replication');
suptitle('R2-R3');

f2 = figure;
imagesc(winTime(1:numKRT), winTime(1:numCT), abs(orig-rep));
title(sprintf('Difference between Original and Replication\nR2-R3'));

f3 = figure;
plot(winTime(1:numKRT), orig(:,17));
hold on
plot(winTime(1:numKRT), rep(:,18), 'r');
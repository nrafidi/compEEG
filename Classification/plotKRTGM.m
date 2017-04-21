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

repTGM_R = squeeze(mean(repTGM_R, 1));
repTGM_F = squeeze(mean(repTGM_F, 1));
origTGM_R = squeeze(mean(origTGM_R, 1));
origTGM_F = squeeze(mean(origTGM_F, 1));

[numRounds, numKRT, numCT] = size(repTGM_R);

%% Plot Each Round Separately

for r = 1:numRounds
    
    orig = squeeze(origTGM_R(r,:,:) - origTGM_F(r,:,:));
    maxO = max(max(abs(orig)));
    rep = squeeze(repTGM_R(r,:,:) - repTGM_F(r,:,:));
    maxR = max(max(abs(rep)));
    
    trueMax = max([maxO, maxR]);
    
    f1 = figure;
    subplot(1,2,1);
    imagesc(winTime(1:numKRT), winTime(1:numCT), orig, [-trueMax, trueMax]);
    title('Original');
    subplot(1,2,2);
    imagesc(winTime(1:numKRT), winTime(1:numCT), rep, [-trueMax, trueMax]);
    title('Replication');
    suptitle(sprintf('Round %d', r));
    
    f2 = figure;
    diffMatrix = abs(orig-rep)/trueMax;
    imagesc(winTime(1:numKRT), winTime(1:numCT), diffMatrix);
    title(sprintf('Difference between Original and Replication\nRound %d', r));
    
    f3 = figure;
    plot(winTime(1:numKRT), orig(:,17));
    hold on
    plot(winTime(1:numKRT), rep(:,18), 'r');
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

f3 = figure;
plot(winTime(1:numKRT), orig(:,17));
hold on
plot(winTime(1:numKRT), rep(:,18), 'r');

%% Plot meanR3R4

orig = squeeze(mean(origTGM_R(3:4,:,:),1) - mean(origTGM_F(3:4,:,:),1));
maxO = max(max(abs(orig)));
rep = squeeze(mean(repTGM_R(3:4,:,:),1) - mean(repTGM_F(3:4,:,:),1));
maxR = max(max(abs(rep)));

trueMax = max([maxO, maxR]);

f1 = figure;
subplot(1,2,1);
imagesc(winTime(1:numKRT), winTime(1:numCT), orig, [-trueMax, trueMax]);
hold on
line([-100, 960], [500, 500], 'Color', 'k');
line([-100, 960], [700, 700], 'Color', 'k');
hold off
title('Original');
subplot(1,2,2);
imagesc(winTime(1:numKRT), winTime(1:numCT), rep, [-trueMax, trueMax]);
hold on
line([-100, 960], [500, 500], 'Color', 'k');
line([-100, 960], [700, 700], 'Color', 'k');
hold off
title('Replication');
suptitle('Mean R3R4');
export_fig(f1, sprintf('%s/results/figures/KR_TGM_OvR_R3R4.png', dataRootR));


f2 = figure;
diffMatrix = abs(orig-rep);
newMatrix = zeros(size(diffMatrix));
newMatrix(diffMatrix <= 0.1 & abs(orig) > 0.019 & abs(rep) > 0.013) = 1;
imagesc(winTime(1:numKRT), winTime(1:numCT), diffMatrix);
title(sprintf('Difference between Original and Replication\nMean R3R4'));
export_fig(f2, sprintf('%s/results/figures/KR_TGM_Diff_R3R4.png', dataRootR));

[U, S, V] = svd(diffMatrix);
lowRankDiff = U(:,1:2)*S(1:2,1:2)*V(:,1:2)';
f25 = figure;
imagesc(winTime(1:numKRT), winTime(1:numCT), lowRankDiff);

f3 = figure;
corrOR = corr(orig, rep);
imagesc(winTime(1:numCT), winTime(1:numCT), corrOR, [0, 1]);
xlabel('Comp Time');
ylabel('Comp Time');
title(sprintf('Correlation between Original and Replication\nMean R3R4'));
export_fig(f3, sprintf('%s/results/figures/KR_TGM_Corr_R3R4.png', dataRootR));

[maxCorr, indMax] = max(corrOR(:));
[i, j] = ind2sub([numCT, numCT], indMax);

f4 = figure;
plot(winTime(1:numKRT), orig(:,i));
hold on
plot(winTime(1:numKRT), rep(:,j), 'r');
% line([500, 500], [-0.08, 0.06], 'Color', 'k');
% line([700, 700], [-0.08, 0.06], 'Color', 'k');
xlabel('KR Time');
ylabel(sprintf('Difference between Remembered and Forgotten'));
title(sprintf('R - F, Mean R3R4\nTrained on Comp Time (%0.0f,%0.0f)', winTime(i), winTime(j)));
legend({'Original', 'Replication'});
export_fig(f4, sprintf('%s/results/figures/KR_TGM_MaxCorr_R3R4.pdf', dataRootR));


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
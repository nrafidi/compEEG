% plot KR Slide
clear;
subjects = {'AA', 'BB', 'DD', 'EE', 'F', 'GG', 'HH', 'JJ', ...
    'K', 'M', 'O', 'R', 'S', 'T', 'U', 'V', 'X', 'Y', 'Z', 'N'};
numSub = length(subjects);
compWinToUse = 17;%17, 24:28, 38

timePoints = -100:20:960;
numTimePoints = length(timePoints);
numRounds = 4;

trajAll = [];%nan(numSub, numRounds, numTimePoints);
trajCorr = cell(numRounds, 1);%nan(numSub, numRounds, numTimePoints);
trajInc = cell(numRounds, 1);%nan(numSub, numRounds, numTimePoints);
trajSubCorr = [];%nan(numSub, numRounds, numTimePoints);
trajSubInc = [];%nan(numSub, numRounds, numTimePoints);
for s = 1:numSub
    
    load(['../../compEEG-data/results/' subjects{s} ...
        '/KRanalysis_SlidingFeat_CWin' num2str(compWinToUse) '.mat']);
    trajInds = logical(responseTraj);
    labels = logical(krLabels);
    trajAll = cat(1, trajAll, krTraj);
    trajSubCorr = cat(1, trajSubCorr, krTraj(labels,:,:));
    trajSubInc = cat(1, trajSubInc, krTraj(~labels,:,:));
    for r = 1:numRounds
        trajCorr{r} = cat(1, trajCorr{r}, krTraj(trajInds(:,r), r, :));
        trajInc{r} = cat(1, trajInc{r}, krTraj(~trajInds(:,r), r, :));
    end
end

%% Corr/inc on same plot
for r = 1:numRounds
    f = figure;
    meanTrajCorr = squeeze(nanmean(trajCorr{r}, 1));
    stdTrajCorr = squeeze(nanstd(trajCorr{r},0,1))./sqrt(size(trajCorr{r}, 1));
    meanTrajInc = squeeze(nanmean(trajInc{r}, 1));
    stdTrajInc = squeeze(nanstd(trajInc{r},0,1))./sqrt(size(trajInc{r}, 1));
    
    plot(timePoints, meanTrajCorr, 'b');
    hold on;
    plot(timePoints, meanTrajInc, 'r');
    legend({'Correct', 'Incorrect'});
    errorbar(timePoints, meanTrajCorr, stdTrajCorr, 'Color', 'b');
    errorbar(timePoints, meanTrajInc, stdTrajInc, 'Color', 'r');
    ylim([0.3 0.7]);
    xlim([min(timePoints), max(timePoints)]);
    ylabel('P(C)');
    xlabel('Start of Window Relative to Stimulus Onset (ms)');
    title(sprintf('Correct vs Incorrect Trials\nRound %d\nTrained on Comp Time %d ms', r, timePoints(compWinToUse)));
    set(gcf, 'Color', 'w');
    export_fig(f, sprintf('../../compEEG-data/results/figures/krSlide_IC_round%d_CWin%d.pdf', r, compWinToUse));
end

%Subseqeuently rem/for on same plot
for r = 1:numRounds
    f = figure;
    meanTrajCorr = squeeze(nanmean(trajSubCorr(:,r,:), 1));
    stdTrajCorr = squeeze(nanstd(trajSubCorr(:,r,:),0,1))./sqrt(size(trajSubCorr, 1));
    meanTrajInc = squeeze(nanmean(trajSubInc(:,r,:), 1));
    stdTrajInc = squeeze(nanstd(trajSubInc(:,r,:),0,1))./sqrt(size(trajSubInc, 1));
    
    plot(timePoints, meanTrajCorr, 'b');
    hold on;
    plot(timePoints, meanTrajInc, 'r');
    legend({'Remembered', 'Forgotten'});
    errorbar(timePoints, meanTrajCorr, stdTrajCorr, 'Color', 'b');
    errorbar(timePoints, meanTrajInc, stdTrajInc, 'Color', 'r');
    ylim([0.3 0.7]);
    xlim([min(timePoints), max(timePoints)]);
    ylabel('P(C)');
    xlabel('Start of Window Relative to Stimulus Onset (ms)');
    title(sprintf('Remembered vs Forgotten Items\nRound %d\nTrained on Comp Time %d ms', r, timePoints(compWinToUse)));
    set(gcf, 'Color', 'w');
    export_fig(f, sprintf('../../compEEG-data/results/figures/krSlide_RF_round%d_CWin%d.pdf', r, compWinToUse));
end

%Rounds
f = figure;
colors = 'rgbm';
hold on;
for r = 1:numRounds
    meanTrajCorr = squeeze(nanmean(trajSubCorr(:,r,:), 1));
    stdTrajCorr = squeeze(nanstd(trajSubCorr(:,r,:),0,1))./sqrt(size(trajSubCorr, 1));
    errorbar(timePoints, meanTrajCorr, stdTrajCorr, 'Color', colors(r));
end
ylim([0.3 0.7]);
xlim([min(timePoints), max(timePoints)]);
ylabel('P(C)');
xlabel('Start of Window Relative to Stimulus Onset (ms)');
legend({'Round 1', 'Round 2', 'Round 3', 'Round 4'});
title(sprintf('Subsequently Remembered Items\nTrained on Comp Time %d ms', timePoints(compWinToUse)));
set(gcf, 'Color', 'w');
export_fig(f, sprintf('../../compEEG-data/results/figures/krSlide_Remembered_CWin%d.pdf', compWinToUse));

f = figure;
colors = 'rgbm';
hold on;
for r = 1:numRounds
    meanTrajInc = squeeze(nanmean(trajSubInc(:,r,:), 1));
    stdTrajInc = squeeze(nanstd(trajSubInc(:,r,:),0,1))./sqrt(size(trajSubInc, 1));
    errorbar(timePoints, meanTrajInc, stdTrajInc, 'Color', colors(r));
    
end
ylim([0.3 0.7]);
xlim([min(timePoints), max(timePoints)]);
ylabel('P(C)');
xlabel('Start of Window Relative to Stimulus Onset (ms)');
legend({'Round 1', 'Round 2', 'Round 3', 'Round 4'});

title(sprintf('Subsequently Forgotten Items\nTrained on Comp Time %d ms', timePoints(compWinToUse)));
set(gcf, 'Color', 'w');
export_fig(f, sprintf('../../compEEG-data/results/figures/krSlide_Forgotten_CWin%d.pdf', compWinToUse));

%Mean R1R2 - R3R4
f = figure;
colors = 'rb';
hold on;
for r = 1:2
    rInd = (r-1)*2 + 1;
    meanTrajCorr = squeeze(nanmean(nanmean(trajSubCorr(:,rInd:(rInd+1),:), 2), 1));
    stdTrajCorr = squeeze(nanstd(nanmean(trajSubCorr(:,rInd:(rInd+1),:), 2),0,1))./sqrt(size(trajSubCorr, 1));
    errorbar(timePoints, meanTrajCorr, stdTrajCorr, 'Color', colors(r));
end
ylim([0.3 0.7]);
xlim([min(timePoints), max(timePoints)]);
ylabel('P(C)');
xlabel('Start of Window Relative to Stimulus Onset (ms)');
legend({'Mean Round 1&2', 'Mean Round 3&4'});
title(sprintf('Subsequently Remembered Items\nTrained on Comp Time %d ms', timePoints(compWinToUse)));
set(gcf, 'Color', 'w');
export_fig(f, sprintf('../../compEEG-data/results/figures/krSlide_Remembered_Mean_CWin%d.pdf', compWinToUse));

f = figure;
colors = 'rb';
hold on;
for r = 1:2
    rInd = (r-1)*2 + 1;
    meanTrajInc = squeeze(nanmean(nanmean(trajSubInc(:,rInd:(rInd+1),:), 2), 1));
    stdTrajIc = squeeze(nanstd(nanmean(trajSubInc(:,rInd:(rInd+1),:), 2),0,1))./sqrt(size(trajSubInc, 1));
    errorbar(timePoints, meanTrajInc, stdTrajInc, 'Color', colors(r));
    
end
ylim([0.3 0.7]);
xlim([min(timePoints), max(timePoints)]);
ylabel('P(C)');
xlabel('Start of Window Relative to Stimulus Onset (ms)');
legend({'Mean Round 1&2', 'Mean Round 3&4'});
title(sprintf('Subsequently Forgotten Items\nTrained on Comp Time %d ms', timePoints(compWinToUse)));
set(gcf, 'Color', 'w');
export_fig(f, sprintf('../../compEEG-data/results/figures/krSlide_Forgotten_Mean_CWin%d.pdf', compWinToUse));

%Mean R3R4
f = figure;
colors = 'rb';
hold on;
meanTrajCorr = squeeze(nanmean(nanmean(trajSubCorr(:,3:4,:), 2), 1));
stdTrajCorr = squeeze(nanstd(nanmean(trajSubCorr(:,3:4,:), 2),0,1))./sqrt(size(trajSubCorr, 1));
errorbar(timePoints, meanTrajCorr, stdTrajCorr, 'Color', colors(1));
meanTrajInc = squeeze(nanmean(nanmean(trajSubInc(:,3:4,:), 2), 1));
stdTrajInc= squeeze(nanstd(nanmean(trajSubInc(:,3:4,:), 2),0,1))./sqrt(size(trajSubInc, 1));
errorbar(timePoints, meanTrajInc, stdTrajInc, 'Color', colors(2));


ylim([0.3 0.7]);
xlim([min(timePoints), max(timePoints)]);
ylabel('P(C)');
xlabel('Start of Window Relative to Stimulus Onset (ms)');
legend({'Remembered', 'Forgotten'});
title(sprintf('Mean R3 and R4 Competition Level over Time\nTrained on Comp Time %d ms', timePoints(compWinToUse)));
set(gcf, 'Color', 'w');
export_fig(f, sprintf('../../compEEG-data/results/figures/krSlide_RvF_MeanR3R4_CWin%d.pdf', compWinToUse));

%% MeanDiff
f = figure;
colors = 'rb';
hold on;
meanTrajCorr12 = squeeze(nanmean(trajSubCorr(:,1:2,:), 2));
meanTrajCorr34 = squeeze(nanmean(trajSubCorr(:,3:4,:), 2));
meanTrajCorr = nanmean(meanTrajCorr12 - meanTrajCorr34);
stdTrajCorr = nanstd(meanTrajCorr12 - meanTrajCorr34,0,1)./sqrt(size(trajSubCorr, 1));
errorbar(timePoints, meanTrajCorr, stdTrajCorr, 'Color', colors(1));
meanTrajInc12 = squeeze(nanmean(trajSubInc(:,1:2,:), 2));
meanTrajInc34 = squeeze(nanmean(trajSubInc(:,3:4,:), 2));
meanTrajInc = nanmean(meanTrajInc12 - meanTrajInc34);
stdTrajInc = nanstd(meanTrajInc12 - meanTrajInc34,0,1)./sqrt(size(trajSubInc, 1));
errorbar(timePoints, meanTrajInc, stdTrajInc, 'Color', colors(2));


ylim([-0.1 0.1]);
xlim([min(timePoints), max(timePoints)]);
line([min(timePoints), max(timePoints)], [0, 0], 'LineStyle', '--', 'Color', 'k');
ylabel('P(C)');
xlabel('Start of Window Relative to Stimulus Onset (ms)');
legend({'Remembered', 'Forgotten'});
title(sprintf('MeanR1R2-MeanR3R4 Competition Level over Time\nTrained on Comp Time %d ms', timePoints(compWinToUse)));
set(gcf, 'Color', 'w');
export_fig(f, sprintf('../../compEEG-data/results/figures/krSlide_RvF_MeanDiff_CWin%d.pdf', compWinToUse));
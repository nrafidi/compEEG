% plot KR Slide
clear;
subjects = {'AA', 'BB', 'DD', 'EE', 'F', 'GG', 'HH', 'JJ', ...
    'K', 'M', 'O', 'R', 'S', 'T', 'U', 'V', 'X', 'Y', 'Z'};
numSub = length(subjects);


timePoints = -280:20:780;
preStimTime = timePoints < 0;
numTimePoints = length(timePoints);
winWidths = 50;
numWinWidths = length(winWidths);
numRounds = 4;

trajAll = nan(numSub, numRounds, numWinWidths, numTimePoints);
trajCorr = nan(numSub, numRounds, numWinWidths, numTimePoints);
trajInc = nan(numSub, numRounds, numWinWidths, numTimePoints);
trajSubCorr = nan(numSub, numRounds, numWinWidths, numTimePoints);
trajSubInc = nan(numSub, numRounds, numWinWidths, numTimePoints);
for s = 1:numSub
    
    load(['../../compEEG-data/results/' subjects{s} ...
        '/KRanalysis_SlidingFeat_lateComp_preStim.mat']);
    trajInds = logical(responseTraj);
    for r = 1:numRounds
        for w = 1:numWinWidths
            baseLineSig = squeeze(mean(krTraj{w}(:,r,preStimTime), 3));
            krTraj{w}(:,r,:) = krTraj{w}(:,r,:) - repmat(baseLineSig, 1,1,numTimePoints);
            labels = logical(krLabels{w});
            trajCorr(s, r, w, :) = squeeze(mean(krTraj{w}(trajInds(:,r), r, :)));
            trajInc(s, r, w, :) = squeeze(mean(krTraj{w}(~trajInds(:,r), r, :)));
            trajAll(s, r, w, :) = (trajCorr(s, r, w, :)+trajInc(s, r, w, :))./2;
            trajSubCorr(s, r, w, :) = squeeze(mean(krTraj{w}(labels, r, :)));
            trajSubInc(s, r, w, :) = squeeze(mean(krTraj{w}(~labels, r, :)));
        end
    end
end

%Corr/inc on same plot
for r = 1:numRounds
    f = figure;
    for w = 1:numWinWidths
        meanTrajCorr = squeeze(nanmean(trajCorr(:,r,w,:), 1));
        stdTrajCorr = squeeze(nanstd(trajCorr(:,r,w,:),0,1))./sqrt(numSub);
        meanTrajInc = squeeze(nanmean(trajInc(:,r,w,:), 1));
        stdTrajInc = squeeze(nanstd(trajInc(:,r,w,:),0,1))./sqrt(numSub);
        if numWinWidths > 1
            subplot(2,1,w);
        end
        plot(timePoints, meanTrajCorr, 'b');
        hold on;
        plot(timePoints, meanTrajInc, 'r');
        legend({'Correct', 'Incorrect'});
        errorbar(timePoints, meanTrajCorr, stdTrajCorr, 'Color', 'b');
        errorbar(timePoints, meanTrajInc, stdTrajInc, 'Color', 'r');
        ylim([-0.2 0.2]);
        xlim([min(timePoints), max(timePoints)]);
        ylabel('P(C)');
        xlabel('Start of Window (ms)');
        title(sprintf('Window Size %d', winWidths(w)));
    end
    suptitle(sprintf('Round %d', r));
    saveas(f, sprintf('../../compEEG-data/results/figures/krSlide_IC_round%d_preStim_Base.png', r));
end

%Subseqeuently corr/inc on same plot
for r = 1:numRounds
    f = figure;
    for w = 1:numWinWidths
        meanTrajCorr = squeeze(nanmean(trajSubCorr(:,r,w,:), 1));
        stdTrajCorr = squeeze(nanstd(trajSubCorr(:,r,w,:),0,1))./sqrt(numSub);
        meanTrajInc = squeeze(nanmean(trajSubInc(:,r,w,:), 1));
        stdTrajInc = squeeze(nanstd(trajSubInc(:,r,w,:),0,1))./sqrt(numSub);
        if numWinWidths > 1
            subplot(2,1,w);
        end
        plot(timePoints, meanTrajCorr, 'b');
        hold on;
        plot(timePoints, meanTrajInc, 'r');
        legend({'Remembered', 'Forgotten'});
        errorbar(timePoints, meanTrajCorr, stdTrajCorr, 'Color', 'b');
        errorbar(timePoints, meanTrajInc, stdTrajInc, 'Color', 'r');
        ylim([-0.2 0.2]);
        xlim([min(timePoints), max(timePoints)]);
        ylabel('P(C)');
        xlabel('Start of Window (ms)');
        title(sprintf('Window Size %d', winWidths(w)));
    end
    suptitle(sprintf('Round %d', r));
    saveas(f, sprintf('../../compEEG-data/results/figures/krSlide_RF_round%d_preStim_Base.png', r));
end

%All
f = figure;

colors = 'rgbm';
for w = 1:numWinWidths
    if numWinWidths > 1
        subplot(2,1,w);
    end
    hold on;
    meanTrajAll = cell(numRounds, 1);
    stdTrajAll = cell(numRounds, 1);
    for r = 1:numRounds
        meanTrajAll{r} = squeeze(nanmean(trajAll(:,r,w,:), 1));
        stdTrajAll{r} = squeeze(nanstd(trajAll(:,r,w,:),0,1))./sqrt(numSub);
        plot(timePoints, meanTrajAll{r}, colors(r));
        
    end
    legend({'R1', 'R2', 'R3', 'R4'}, 'Location', 'northeast');
    for r = 1:numRounds
        errorbar(timePoints, meanTrajAll{r}, stdTrajAll{r}, 'Color', colors(r));
    end
    ylim([-0.2 0.2]);
    xlim([min(timePoints), max(timePoints)]);
    ylabel('P(C)');
    xlabel('Start of Window (ms)');
    title(sprintf('Window Size %d', winWidths(w)));
    
end
suptitle('All Samples');
saveas(f, sprintf('../../compEEG-data/results/figures/krSlide_all_preStim_Base.png'));

f = figure;
colors = 'rgbm';
for w = 1:numWinWidths
    if numWinWidths > 1
        subplot(2,1,w);
    end
    hold on;
    for r = 1:numRounds
        meanTrajCorr = squeeze(nanmean(trajSubCorr(:,r,w,:), 1));
        stdTrajCorr = squeeze(nanstd(trajSubCorr(:,r,w,:),0,1))./sqrt(numSub);
        
        errorbar(timePoints, meanTrajCorr, stdTrajCorr, 'Color', colors(r));
        
    end
    ylim([-0.2 0.2]);
    xlim([min(timePoints), max(timePoints)]);
    ylabel('P(C)');
    xlabel('Start of Window (ms)');
    title(sprintf('Window Size %d', winWidths(w)));
    legend({'Round 1', 'Round 2', 'Round 3', 'Round 4'});
end
suptitle('Subsequently Remembered');
saveas(f, sprintf('../../compEEG-data/results/figures/krSlide_Remembered_preStim_Base.png'));

f = figure;
colors = 'rgbm';
for w = 1:numWinWidths
    if numWinWidths > 1
        subplot(2,1,w);
    end
    hold on;
    for r = 1:numRounds
        meanTrajInc = squeeze(nanmean(trajSubInc(:,r,w,:), 1));
        stdTrajInc = squeeze(nanstd(trajSubInc(:,r,w,:),0,1))./sqrt(numSub);
        
        errorbar(timePoints, meanTrajInc, stdTrajInc, 'Color', colors(r));
        
    end
    ylim([-0.2 0.2]);
    xlim([min(timePoints), max(timePoints)]);
    ylabel('P(C)');
    xlabel('Start of Window (ms)');
    title(sprintf('Window Size %d', winWidths(w)));
    legend({'Round 1', 'Round 2', 'Round 3', 'Round 4'});
end
suptitle('Subsequently Forgotten');
saveas(f, sprintf('../../compEEG-data/results/figures/krSlide_Forgotten_preStim_Base.png'));
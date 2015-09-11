% plot KR Slide
subjects = {'AA', 'BB', 'DD', 'EE', 'F', 'GG', 'HH', 'JJ', ...
    'K', 'M', 'N', 'O', 'R', 'S', 'T', 'U', 'V', 'X', 'Y', 'Z'};
numSub = length(subjects);


timePoints = 0:20:780;
numTimePoints = length(timePoints);
winWidths = 50:50:200;
numWinWidths = length(winWidths);
numRounds = 4;

trajAll = nan(numSub, numRounds, numWinWidths, numTimePoints);
trajCorr = nan(numSub, numRounds, numWinWidths, numTimePoints);
trajInc = nan(numSub, numRounds, numWinWidths, numTimePoints);
trajSubCorr = nan(numSub, numRounds, numWinWidths, numTimePoints);
trajSubInc = nan(numSub, numRounds, numWinWidths, numTimePoints);
for s = 1:numSub
    
    load(['../../compEEG-data/results/' subjects{s} ...
        '/KRanalysis_SlidingFeat.mat']);
    trajInds = logical(responseTraj);
    for r = 1:numRounds
        for w = 1:numWinWidths
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
        subplot(2,2,w);
        plot(timePoints, meanTrajCorr, 'b');
        hold on;
        plot(timePoints, meanTrajInc, 'r');
        legend({'Correct', 'Incorrect'});
        errorbar(timePoints, meanTrajCorr, stdTrajCorr, 'Color', 'b');
        errorbar(timePoints, meanTrajInc, stdTrajInc, 'Color', 'r');
        ylim([0.3 0.7]);
        xlim([min(timePoints), max(timePoints)]);
        ylabel('P(C)');
        xlabel('Start of Window (ms)');
        title(sprintf('Window Size %d', winWidths(w)));
    end
    suptitle(sprintf('Round %d', r));
    saveas(f, sprintf('../../compEEG-data/results/figures/krSlide_IC_round%d.png', r));
end

%Subseqeuently corr/inc on same plot
for r = 1:numRounds
    f = figure;
    for w = 1:numWinWidths
        meanTrajCorr = squeeze(nanmean(trajSubCorr(:,r,w,:), 1));
        stdTrajCorr = squeeze(nanstd(trajSubCorr(:,r,w,:),0,1))./sqrt(numSub);
        meanTrajInc = squeeze(nanmean(trajSubInc(:,r,w,:), 1));
        stdTrajInc = squeeze(nanstd(trajSubInc(:,r,w,:),0,1))./sqrt(numSub);
        subplot(2,2,w);
        plot(timePoints, meanTrajCorr, 'b');
        hold on;
        plot(timePoints, meanTrajInc, 'r');
        legend({'Remembered', 'Forgotten'});
        errorbar(timePoints, meanTrajCorr, stdTrajCorr, 'Color', 'b');
        errorbar(timePoints, meanTrajInc, stdTrajInc, 'Color', 'r');
        ylim([0.3 0.7]);
        xlim([min(timePoints), max(timePoints)]);
        ylabel('P(C)');
        xlabel('Start of Window (ms)');
        title(sprintf('Window Size %d', winWidths(w)));
    end
    suptitle(sprintf('Round %d', r));
    saveas(f, sprintf('../../compEEG-data/results/figures/krSlide_RF_round%d.png', r));
end

%All
for r = 1:numRounds
    f = figure;
    for w = 1:numWinWidths
        meanTrajAll = squeeze(nanmean(trajAll(:,r,w,:), 1));
        stdTrajAll = squeeze(nanstd(trajAll(:,r,w,:),0,1))./sqrt(numSub);
        subplot(2,2,w);
        plot(timePoints, meanTrajAll);
        errorbar(timePoints, meanTrajAll, stdTrajAll);
        ylim([0.3 0.7]);
        xlim([min(timePoints), max(timePoints)]);
        ylabel('P(C)');
        xlabel('Start of Window (ms)');
        title(sprintf('Window Size %d', winWidths(w)));
    end
    suptitle(sprintf('Round %d', r));
    saveas(f, sprintf('../../compEEG-data/results/figures/krSlide_all_round%d.png', r));
end
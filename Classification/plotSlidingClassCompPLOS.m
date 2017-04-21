% Plot sliding classification data
load ../../compEEG-data/results/CompEEG_5FCV_Slide_Accs.mat
subAccs = squeeze(mean(subAccs, 3));
[numSub, numT] = size(subAccs);
numT = numT-1;
winTime = winTime(1:numT);
subAccs = subAccs(:,1:numT);

erpSubAvg = nanmean(subAccs);
erpSubStd = nanstd(subAccs)./sqrt(numSub);
subUpp = erpSubAvg + erpSubStd;
subLow = erpSubAvg - erpSubStd;
X = [winTime, fliplr(winTime)];
Y = [subLow, fliplr(subUpp)];
f = figure;
hold on;
fh = fill(X, Y, 'c');
set(fh, 'EdgeAlpha', 0);
plot(winTime, erpSubAvg, 'b');
line([min(winTime) - 20, max(winTime)+20], [0.5, 0.5], 'Color', 'k');
line([min(winTime) - 20, max(winTime)+20], [0.55, 0.55], 'Color', 'k', 'LineStyle', '--');
line([0, 0], [0.46, 0.66], 'Color', 'k');
% scatter(220, erpSubAvg(17), 500, 'ro', 'LineWidth', 5);
% errorbar(winTime, erpSubAvg, erpSubStd, 'b');

xlabel('Start Time of Window Relative to Onset (ms)');
xlim([min(winTime) - 20, max(winTime)+20]);
ylim([0.46, 0.66]);

ylabel('Classification Accuracy');
title(sprintf('Competition accuracy over time\n%s', ...
    'using 50ms window averages'));
set(gcf, 'color', 'w');
set(gca, 'fontsize', 16);
export_fig(f, '../../compEEG-data/results/figures/compCV_allSub_50ms_PLOS.pdf');
% 
% freqSubAvg = squeeze(mean(collectFreq, 4));
% freqSubStd = (squeeze(std(collectFreq, 0, 4)))./sqrt(numSub);
% 
% f = figure;
% freqTitles = {'\alpha', '\alpha', '\alpha'};
% for freq = 1:numFreqTypes
%     subplot(numFreqTypes, 1, freq);
%     hold on;
%     for wWd = 1:numWinWidths
%         plot(timeForFreq, freqSubAvg(freq,:, wWd), colors(wWd));
%     end
%     legend({'50ms', '100ms', '150ms', '200ms'});
%     line([-20, 800], [0.5, 0.5], 'Color',  'k');
%     for wWd = 1:numWinWidths
%         errorbar(timeForFreq, freqSubAvg(freq,:, wWd), freqSubStd(freq,:, wWd), colors(wWd));
%     end
%     xlabel('Start Time of Window Relative to Onset (ms)');
%     xlim([-20, 800]);
%     ylim([0.45, 0.7]);
%     
%     ylabel('Classification Accuracy');
%     title(freqTitles{freq});
% end
% 
% suptitle(sprintf('Competition accuracy over time\n%s', ...
%     'using different sized window averages of freq. power'));
% saveas(f, '../../compEEG-data/results/figures/freq_allSub_allWidths.png');
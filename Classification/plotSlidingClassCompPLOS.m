% Plot sliding classification data
load ../../compEEG-data/results/CompEEG_CV_Slide_Accs_PDTW.mat
[numSub, numT] = size(subAccs);
winTime = winTime(1:numT);

erpSubAvg = mean(subAccs);
erpSubStd = std(subAccs)./sqrt(numSub);
f = figure;
hold on;
plot(winTime, erpSubAvg, 'b');
line([min(winTime) - 20, max(winTime)+20], [0.5, 0.5], 'Color', 'k');
line([min(winTime) - 20, max(winTime)+20], [0.55, 0.55], 'Color', 'k', 'LineStyle', '--');
errorbar(winTime, erpSubAvg, erpSubStd, 'b');

xlabel('Start Time of Window Relative to Onset (ms)');
xlim([min(winTime) - 20, max(winTime)+20]);
ylim([0.48, 0.68]);

ylabel('Classification Accuracy');
title(sprintf('Competition accuracy over time\n%s', ...
    'using 50ms window averages'));
set(gcf, 'color', 'w');
set(gca, 'fontsize', 16);
export_fig(f, '../../compEEG-data/results/figures/compCV_allSub_50ms_PDTW_PLOS.pdf');
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
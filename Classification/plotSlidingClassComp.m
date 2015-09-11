% Plot sliding classification data
numWinWidths = 4;
numERPwin = 40;
numFreqTypes = 3;
numFreqWin = 40;

timeForERP = 0:20:780;
timeForFreq = 0:20:780;

colors = 'rgbm';

filesToPlot = dir('../../compEEG-data/results/CompEEG_CV_*_Def_Slide_All.mat');
numSub = length(filesToPlot);

collectERP = nan(numERPwin, numWinWidths, numSub);
collectFreq = nan(numFreqTypes, numFreqWin, numWinWidths, numSub);

% subjects = {'L', 'P', 'W', 'CC', 'FF', 'H', 'I', 'J', ...
%     'AA', 'BB', 'DD', 'EE', 'GG', 'HH', 'JJ', 'M', 'N',...
%     'F', 'K', 'O', 'R', 'S', 'T', 'U', 'V', 'X', 'Y', 'Z'};
for s = 1:numSub
    %     subFig = figure;
    load(['../../compEEG-data/results/' filesToPlot(s).name]);
    
    collectERP(:, :, s) = subAcc(1:numERPwin, :);%mean(subAcc(1:numERPwin, :), 2);
    
    for f = 1:numFreqTypes
        collectFreq(f, :, :, s) = ...
            subAcc(((f-1)*numFreqWin + 1 + numERPwin):(f*numFreqWin + numERPwin), :);
        % mean(subAcc((numERPwin+f):numFreqTypes:end,:), 2);
        %         subFStd = cat(2, subFStd, std(subAcc((numERPwin+f):numFreqTypes:end,:), 0, 2));
        %         plot(timeForFreq, squeeze(collectFreq(f,:,s)), colors(f));
    end
    
    %     subStd = std(subAcc(1:numERPwin, :),0, 2);
    %
    %     subplot(2, 1, 1);
    %     plot(timeForERP, collectERP(:,s));
    %     errorbar(timeForERP, collectERP(:,s), subStd);
    %     xlabel('Start Time of Window Relative to Onset (ms)');
    %     xlim([-50, 850]);
    %     ylim([0.45, 0.7]);
    %     line([-50, 850], [0.5, 0.5], 'Color', 'k');
    %     ylabel('Classification Accuracy');
    %     title(sprintf('Competition accuracy over time\n%s', ...
    %         'using 50ms window averages of voltage'));
    %
    %     subplot(2,1,2);
    %     hold on;
    %     subFStd = [];
    
    
    %     legend({'\theta', '\alpha', '\beta'});
    %
    %     for freq = 1:numFreqTypes
    %         errorbar(timeForFreq, squeeze(collectFreq(freq,:,s)), subFStd(:,freq), colors(freq));
    %     end
    %
    %     xlabel('Start Time of Window Relative to Onset (ms)');
    %     ylabel('Classification Accuracy');
    %     xlim([-50, 650]);
    %     ylim([0.45, 0.7]);
    %     line([-50, 650], [0.5, 0.5], 'Color',  'k');
    %     title(sprintf('Competition accuracy over time\n%s', ...
    %         'using 200ms window averages of freq. power'));
    
    
    %     suptitle(subjects{s});
    %     saveas(subFig, ['../../compEEG-data/results/figures/' ...
    %         subjects{s} '_comp_slide.png']);
end

erpSubAvg = squeeze(mean(collectERP, 3));
erpSubStd = (squeeze(std(collectERP, 0, 3)))./sqrt(numSub);
f = figure;
hold on;
for wWd = 1:numWinWidths
    plot(timeForERP, erpSubAvg(:, wWd), colors(wWd));
end
line([-20, 800], [0.5, 0.5], 'Color', 'k');
for wWd = 1:numWinWidths
    errorbar(timeForERP, erpSubAvg(:, wWd), erpSubStd(:,wWd), colors(wWd));
end
xlabel('Start Time of Window Relative to Onset (ms)');
xlim([-20, 800]);
ylim([0.45, 0.7]);

ylabel('Classification Accuracy');
title(sprintf('Competition accuracy over time\n%s', ...
    'using different sized window averages of voltage'));
saveas(f, '../../compEEG-data/results/figures/voltage_allSub_allWidths.png');

freqSubAvg = squeeze(mean(collectFreq, 4));
freqSubStd = (squeeze(std(collectFreq, 0, 4)))./sqrt(numSub);

f = figure;
freqTitles = {'\theta', '\alpha', '\beta'};
for freq = 1:numFreqTypes
    subplot(numFreqTypes, 1, freq);
    hold on;
    for wWd = 1:numWinWidths
        plot(timeForFreq, freqSubAvg(freq,:, wWd), colors(wWd));
    end
    legend({'50ms', '100ms', '150ms', '200ms'});
    line([-20, 800], [0.5, 0.5], 'Color',  'k');
    for wWd = 1:numWinWidths
        errorbar(timeForFreq, freqSubAvg(freq,:, wWd), freqSubStd(freq,:, wWd), colors(wWd));
    end
    xlabel('Start Time of Window Relative to Onset (ms)');
    xlim([-20, 800]);
    ylim([0.45, 0.7]);
    
    ylabel('Classification Accuracy');
    title(freqTitles{freq});
end

suptitle(sprintf('Competition accuracy over time\n%s', ...
    'using different sized window averages of freq. power'));
saveas(f, '../../compEEG-data/results/figures/freq_allSub_allWidths.png');
% Plot sliding classification data

numERPwin = 16;
numFreqTypes = 3;
numFreqWin = 4;

timeForERP = 0:50:750;
timeForFreq = 0:200:600;

colors = 'rgb';

filesToPlot = dir('../../compEEG-data/results/CompEEG_CV_*_Slide.mat');
numSub = length(filesToPlot);

collectERP = nan(numERPwin, numSub);
collectFreq = nan(numFreqTypes, numFreqWin, numSub);

subjects = {'CC', 'FF', 'H', 'I', 'J', 'L', 'P', 'W'};
for s = 1:numSub
    subFig = figure;
    load(['../../compEEG-data/results/' filesToPlot(s).name]);
    
    collectERP(:, s) = mean(subAcc(1:numERPwin, :), 2);
    subStd = std(subAcc(1:numERPwin, :),0, 2);
    
    subplot(2, 1, 1);
    plot(timeForERP, collectERP(:,s));
    errorbar(timeForERP, collectERP(:,s), subStd);
    xlabel('Start Time of Window Relative to Onset (ms)');
    xlim([-50, 850]);
    ylim([0.45, 0.7]);
    line([-50, 850], [0.5, 0.5], 'Color', 'k');
    ylabel('Classification Accuracy');
    title(sprintf('Competition accuracy over time\n%s', ...
        'using 50ms window averages of voltage'));
    
    subplot(2,1,2);
    hold on;
    subFStd = [];
    for f = 1:numFreqTypes
        collectFreq(f, :, s) = ...
            mean(subAcc((numERPwin+f):numFreqTypes:end,:), 2);
        subFStd = cat(2, subFStd, std(subAcc((numERPwin+f):numFreqTypes:end,:), 0, 2));
        plot(timeForFreq, squeeze(collectFreq(f,:,s)), colors(f));
    end
    
    legend({'\theta', '\alpha', '\beta'});

for freq = 1:numFreqTypes
    errorbar(timeForFreq, squeeze(collectFreq(freq,:,s)), subFStd(:,freq), colors(freq));
end

xlabel('Start Time of Window Relative to Onset (ms)');
ylabel('Classification Accuracy');
xlim([-50, 650]);
ylim([0.45, 0.7]);
line([-50, 650], [0.5, 0.5], 'Color',  'k');
title(sprintf('Competition accuracy over time\n%s', ...
    'using 200ms window averages of freq. power'));


suptitle(subjects{s});
saveas(subFig, ['../../compEEG-data/results/figures/' ...
    subjects{s} '_comp_slide.png']);
end

erpSubAvg = mean(collectERP, 2);
erpSubStd = std(collectERP, 0, 2);
f = figure;
plot(timeForERP, erpSubAvg);
errorbar(timeForERP, erpSubAvg, erpSubStd);
xlabel('Start Time of Window Relative to Onset (ms)');
xlim([-50, 850]);
ylim([0.45, 0.7]);
line([-50, 850], [0.5, 0.5], 'Color', 'k');
ylabel('Classification Accuracy');
title(sprintf('Competition accuracy over time\n%s', ...
    'using 50ms window averages of voltage'));
saveas(f, '../../compEEG-data/results/figures/voltage_pilot.png');

freqSubAvg = squeeze(mean(collectFreq, 3));
freqSubStd = squeeze(std(collectFreq, 0, 3));
f = figure;
hold on;

for freq = 1:numFreqTypes
    plot(timeForFreq, freqSubAvg(freq,:), colors(freq));
end
legend({'\theta', '\alpha', '\beta'});

for freq = 1:numFreqTypes
    errorbar(timeForFreq, freqSubAvg(freq,:), freqSubStd(freq,:), colors(freq));
end

xlabel('Start Time of Window Relative to Onset (ms)');
ylabel('Classification Accuracy');
xlim([-50, 650]);
ylim([0.45, 0.7]);
line([-50, 650], [0.5, 0.5], 'Color',  'k');
title(sprintf('Competition accuracy over time\n%s', ...
    'using 200ms window averages of freq. power'));
saveas(f, '../../compEEG-data/results/figures/freq_pilot.png');
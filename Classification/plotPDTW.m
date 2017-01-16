itemAvg = '';%'_ItemAvg';
fileString = ...
    '/Users/nrafidi/Documents/MATLAB/compEEG-data/results/CompEEG_CV_Slide_Accs%s%s.mat';
figString = ...
    '/Users/nrafidi/Documents/MATLAB/compEEG-data/results/figures/CompEEG_CV_Slide_PDTW%s%s.mat';
time= -100:20:1000;

load(sprintf(fileString, '', itemAvg));
oldData = subAccs;
load(sprintf(fileString, '_PDTW', itemAvg));
newData = subAccs;

[numSub, numTime] = size(newData);
numSub = 8;
time = time(1:numTime);

numPlots = ceil(numSub/4);

f = figure;
for s = 1:numSub
    subplot(4, numPlots, s);
    
    plot(time, oldData(s,:), 'b');
    hold on
    plot(time, newData(s,:), 'r');
    xlim([min(time), max(time)]);
    line([min(time), max(time)], [0.5, 0.5], 'Color', 'k', 'LineStyle', '--');
    line([0, 0], [0.4, 0.7], 'Color', 'k');
    
    if s == numSub
        legend('Old', 'New', 'Chance', 'Location', 'NorthEast');
    end
end
set(f, 'Color', 'w');
export_fig(f, sprintf(figString, itemAvg, '_singleSub.png'));

f = figure;

plot(time, mean(oldData), 'b');
hold on
plot(time, mean(newData), 'r');
xlim([min(time), max(time)]);
line([min(time), max(time)], [0.5, 0.5], 'Color', 'k', 'LineStyle', '--');
line([0, 0], [0.4, 0.7], 'Color', 'k');

if s == numSub
    legend('Old', 'New', 'Chance', 'Location', 'NorthEast');
end
title('Distinguishing Retrieval States from EEG Data');
set(f, 'Color', 'w');
export_fig(f, sprintf(figString, itemAvg, '_subAvg.png'));
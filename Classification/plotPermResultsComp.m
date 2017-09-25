%%
threshold = 0.001/55;
resultDir = '/Users/nrafidi/Documents/MATLAB/compEEG-data%s/results/';
addpath /Users/nrafidi/Documents/MATLAB/Toolboxes/export_fig/


resultDirRep = sprintf(resultDir, '-rep');

resultDir = sprintf(resultDir, '');


fileString = 'CompEEG_5FCV_win*_permAccs.mat';


filesRep = dir([resultDirRep fileString]);
files = dir([resultDir fileString]);


%%
numFiles = length(files);
numFilesRep = length(filesRep);

assert(numFiles == numFilesRep);

trueAccs = nan(numFiles, 1);
trueSubAccsPooled = [];
pVals = nan(numFiles, 1);
timeToPlot = nan(numFiles, 1);
for f = 1:numFiles
    load([resultDir files(f).name]);
    permSubAccsOrig = permSubAccs;
    trueSubAccsOrig = trueSubAccs;
    load([resultDirRep filesRep(f).name]);
    permSubAccs = cat(1, permSubAccs, permSubAccsOrig);
    trueSubAccs = cat(1, trueSubAccs, trueSubAccsOrig);
    timeToPlot(f) = winTime(str2double(files(f).name(17:18)));
    trueAccs(f) = mean(mean(trueSubAccs));
    trueSubAccs = mean(trueSubAccs, 2);
    trueSubAccsPooled = cat(2, trueSubAccsPooled, trueSubAccs);
    permAccs = squeeze(mean(permSubAccs, 3));
    subPvals = sum(permAccs > repmat(trueSubAccs, 1, 100), 2)/100;
    subProbs = norminv(subPvals+eps, 0, 1);
    [~, pVals(f)] = ttest(subProbs); 
end

numSub = size(trueSubAccs, 1);

erpSubAvg = nanmean(trueSubAccsPooled);
erpSubStd = nanstd(trueSubAccsPooled)./sqrt(numSub);
subUpp = erpSubAvg + erpSubStd;
subLow = erpSubAvg - erpSubStd;
X = [timeToPlot', fliplr(timeToPlot')];
Y = [subLow, fliplr(subUpp)];

h = figure;
hold on
fh = fill(X, Y, 'c');
set(fh, 'EdgeAlpha', 0);
plot(timeToPlot, erpSubAvg, 'b');
for f = 1:(numFiles-1)
    if pVals(f) <= threshold
        scatter(timeToPlot(f), trueAccs(f)+0.005, 'r*');
    end
end
legend({'Standard Deviation', 'Accuracy', 'Bonferroni Significant'});
line([min(timeToPlot), max(timeToPlot)], [0.55, 0.55], 'LineStyle', '--', 'Color', 'k');
xlim([min(timeToPlot), max(timeToPlot)]);
ylim([min(trueAccs)-0.005, max(trueAccs)+0.01]);
xlabel('Time relative to stimulus onset (ms)');
ylabel('Accuracy');

title(sprintf('Competition Accuracy Over Time'));

set(h, 'Color', 'w');
export_fig(h, [resultDirRep 'figures/compPerm_07-40.pdf']);
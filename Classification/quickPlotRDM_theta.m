%%
resultDir = '/Users/nrafidi/Documents/MATLAB/compEEG-data%s/results/';
addpath /Users/nrafidi/Documents/MATLAB/Toolboxes/export_fig/


repSubjects = {'III', 'JJJ', 'KKK', 'BBB', 'GGG', 'HHH', 'AAA', 'CCC', 'DDD', 'EEE', 'FFF', 'MM', ...
    'OO', 'PP', 'QQ', 'RR', 'SS', 'TT', 'WW',...
    'YY'};
subjects =  {'L', 'P', 'W', 'CC', 'FF', 'H', 'I', 'J', ...
    'AA', 'BB', 'DD', 'EE', 'GG', 'HH', 'JJ', 'M', 'N',...
    'F', 'K', 'O', 'R', 'S', 'T', 'U', 'V', 'X', 'Y', 'Z'};



resultDirRep = sprintf(resultDir, '-rep');

resultDir = sprintf(resultDir, '');

loadFileString = 'CompEEG_%s_CV_Slide_Models_PLOS_theta.mat';
saveFileString = 'CompEEG_CV_Slide_Models_theta.mat';


subModelsOrig = cell(length(subjects), 1);
for s = 1:length(subjects)
    sub = subjects{s};
    loadFile = sprintf(loadFileString, sub);
    load([resultDir loadFile]);
    subModelsOrig{s} = cell2mat(models(:, 1)')';
end
save([resultDir saveFileString], 'subModelsOrig');

subModels = cell(length(repSubjects), 1);
for s = 1:length(repSubjects)
    sub = repSubjects{s};
    loadFile = sprintf(loadFileString, sub);
    load([resultDirRep loadFile]);
    subModels{s} = cell2mat(models(:, 1)')';
end
save([resultDirRep saveFileString], 'subModels');
%%
subModels = cat(1, subModels, subModelsOrig);
numSub = size(subModels, 1);
avgRDM = nan(54);
timeVec = -100:20:960;
minTime = min(timeVec);
maxTime = max(timeVec);

for s = 1:numSub
    subModel = subModels{s};
    subDist = pdist(subModel, 'cosine');
    subDist = squareform(subDist);
    subDist = subDist./max(max(subDist));
    avgRDM = avgRDM + subDist;
end
avgRDM = avgRDM./numSub;
f = figure;
imagesc(timeVec, timeVec, subDist);
hold on
% styles = {'-', '--', '-.'};
% handles = nan(3, 1);
% for mT = 1:length(myTimes)
%     myTime = myTimes(mT);
%     handles(mT) = line([minTime maxTime], [myTime myTime], 'Color', 'k', 'LineStyle', styles{mT});
%     line([myTime myTime], [minTime maxTime], 'Color', 'k', 'LineStyle', styles{mT});
% end
title(sprintf('Average RDM of Weight Vectors over time\nTheta'));
xlabel('Time Relative to Onset (ms)');
ylabel('Time Relative to Onset (ms)');
colorbar
% legend(handles, {'Peak 1', 'Peak 2', 'Peak 3'});

set(gca, 'fontsize', 18);
set(gcf, 'color', 'w');
% set (f, 'Units', 'normalized', 'Position', [0,0,1,1]);
export_fig(f, '../../compEEG-data-rep/results/figures/avgRDM_50ms_theta_PLOS.pdf');

%%
timeToPlot = 7:40;
f = figure;
imagesc(timeVec(timeToPlot), timeVec(timeToPlot), subDist(timeToPlot, timeToPlot));
hold on
styles = {'-', '--', '-.'};
handles = nan(3, 1);
for mT = 1:length(myTimes)
    myTime = myTimes(mT);
    handles(mT) = line([minTime maxTime], [myTime myTime], 'Color', 'k', 'LineStyle', styles{mT});
    line([myTime myTime], [minTime maxTime], 'Color', 'k', 'LineStyle', styles{mT});
end
% title(sprintf('Average RDM of Weight Vectors over time\n Using 50ms window averages'));
% xlabel('Start Time of Window Relative to Onset (ms)');
% ylabel('Start Time of Window Relative to Onset Onset (ms)');
colorbar
% legend(handles, {'Peak 1', 'Peak 2', 'Peak 3'});

set(gca, 'fontsize', 16);
set(gcf, 'color', 'w');
export_fig(f, '../../compEEG-data-rep/results/figures/avgRDM_50ms_PLOS_sig.pdf');
% %%
% for t = 1:length(myTimes)
%     f = figure;
%     timeToPlot = find(myTimes(t) == timeVec);
%     plot(timeVec, avgRDM(:, timeToPlot), 'b');
%     hold on
%     styles = {'-', '--', '-.'};
%     handles = nan(3, 1);
%     for mT = 1:length(myTimes)
%         myTime = myTimes(mT);
%         handles(mT) = line([myTime myTime], [0 1], 'Color', 'k', 'LineStyle', styles{mT});
%     end
%     legend(handles, {'Peak 1', 'Peak 2', 'Peak 3'});
%     xlim([minTime, maxTime]);
%     xlabel('Time Relative to Stimulus Onset (ms)');
%     ylabel('Cosine Distance');
%     title(sprintf('Distance of all weight vectors to %d ms weight vector', myTimes(t)));
%     set(gcf, 'Color', 'w');
%     set(gca, 'fontsize', 16);
% %     export_fig(f, sprintf('../../compEEG-data/results/figures/avgRDMSlice_50ms_%d_PLOS.pdf', myTimes(t)));
% end
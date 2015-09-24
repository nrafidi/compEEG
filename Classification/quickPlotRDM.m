% Quickly plot the RSM for all subs
load ../../compEEG-data/results/CompEEG_CV_Slide_Models.mat

numSub = size(subModels, 1);
avgRDM = zeros(size(subModels, 2));
timeVec = 0:20:780;
myTime = 240;
minTime = min(timeVec);
maxTime = max(timeVec);

for s = 1:numSub
subModel = subModels(s,:);
subModel = cell2mat(subModel);
subModel = subModel';
subDist = pdist(subModel, 'cosine');
subDist = squareform(subDist);
subDist = subDist./max(max(subDist));
avgRDM = avgRDM + subDist;
end
avgRDM = avgRDM./numSub;
f = figure;
imagesc(timeVec, timeVec, subDist);
hold on
line([minTime maxTime], [myTime myTime], 'Color', 'k');
line([myTime myTime], [minTime maxTime], 'Color', 'k');
title(sprintf('Average RDM of Weight Vectors over time\n Using 50ms window averages'));
xlabel('Time Post Onset (ms)');
ylabel('Time Post Onset (ms)');
colorbar
legend({'Time of Peak Decoding Accuracy'});
saveas(f, '../../compEEG-data/results/figures/avgRDM_50ms.png');
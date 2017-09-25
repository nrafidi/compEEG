%Bootstrap cluster value
function [trueClusterT, permClusterT, bootGrid] = bootstrapCluster_KRTGM(clusterToUse, computationToPlot)

dataRootR = '/Users/nrafidi/Documents/MATLAB/compEEG-data-rep/';
load(sprintf('%s/results/clusters_pVals_KRTGM.mat', dataRootR));
numPerms = 1000;
numSubjects = size(IndividualSubjectDataR, 1);

% True Cluster Size
[krTGM_R, krTGM_F, firstRoundCorr_R, firstRoundCorr_F] = collectData(IndividualSubjectDataR, IndividualSubjectDataF, ...
    IndividualSubjectFirstCorrR, IndividualSubjectFirstCorrF, krWinToUse, ctWinToUse);

trueClusterT = scoreCluster(krTGM_R, krTGM_F, firstRoundCorr_R, firstRoundCorr_F, clusters{clusterToUse}, computationToPlot);

permClusterT = nan(numPerms, 1);
numTKR = length(krWinToUse);
numTC = length(ctWinToUse);
bootGrid = zeros(numTKR, numTC);
subjectList = 1:numSubjects;
for p = 1:numPerms
    tic
    subjectSample = datasample(subjectList, numSubjects);
    
    [krTGM_R, krTGM_F, firstRoundCorr_R, firstRoundCorr_F] = collectData( ...
        IndividualSubjectDataR(subjectSample), IndividualSubjectDataF(subjectSample), ...
    IndividualSubjectFirstCorrR(subjectSample), IndividualSubjectFirstCorrF(subjectSample), krWinToUse, ctWinToUse);

    [permClusterT(p), diffMatR, diffMatF] = scoreCluster(krTGM_R, krTGM_F, firstRoundCorr_R, firstRoundCorr_F, clusters{clusterToUse}, computationToPlot);
    
    options.minClusterSize = 1;
    options.pValThresh = 0.05;
    options.clusterTestStatistic = 'summed_t_values';
    dataToCluster = cat(4, diffMatR, diffMatF);
    clustersBoot = findClusters(dataToCluster, options);
    [~, uniClustInd] = unique(cellfun(@num2str, ...
    cellfun(@(x) reshape(x, 1, []), clustersBoot, 'UniformOutput', false), ...
    'UniformOutput', false));
    clustersBoot = clustersBoot(uniClustInd);
    
    clustersBootInds = cellfun(@(x) sub2ind([numTKR, numTC], x(:,1), x(:,2)), ...
        clustersBoot, 'UniformOutput', false);
    
    for iClust = 1:length(clustersBootInds)
        bootGrid(clustersBootInds{iClust}) = bootGrid(clustersBootInds{iClust}) + 1;
    end
    
    toc
    disp(p)
end

bootGrid = bootGrid./numPerms;
f = figure;
imagesc(krWinTime, ctWinTime, bootGrid');
colorbar
xlabel('KR Time (ms)');
ylabel('Competition Time (ms)');
title(sprintf('Proportion of Bootstrap Draws\nCluster Participation'));
export_fig(f, sprintf('%s/results/figures/clusterBootstrapGrid.png', dataRootR));

f = figure;
histogram(permClusterT);
counts = histcounts(permClusterT);
hold on
line([trueClusterT, trueClusterT], [0, max(counts)+5]);
legend({'Bootstrap Distribution', 'True Value'});
xlabel('Cluster Summed T Stat');
ylabel('Number of Bootstrap Draws');

export_fig(f, sprintf('%s/results/figures/clusterBootstrap.png', dataRootR));
export_fig(f, sprintf('%s/results/figures/clusterBootstrap.pdf', dataRootR));
export_fig(f, sprintf('%s/results/figures/clusterBootstrap.fig', dataRootR));
end

function [krTGM_R, krTGM_F, firstRoundCorr_R, firstRoundCorr_F] = collectData(IndividualSubjectDataR, IndividualSubjectDataF, ...
    IndividualSubjectFirstCorrR, IndividualSubjectFirstCorrF, krWinToUse, ctWinToUse)

numSubjects = length(IndividualSubjectDataR);

krTGM_R = [];
krTGM_F = [];
firstRoundCorr_R = [];
firstRoundCorr_F = [];
for s = 1:numSubjects
    krTGM_R = cat(1, krTGM_R, IndividualSubjectDataR{s});
    krTGM_F = cat(1, krTGM_F, IndividualSubjectDataF{s});
    firstRoundCorr_R = cat(1, firstRoundCorr_R, IndividualSubjectFirstCorrR{s});
    firstRoundCorr_F = cat(1, firstRoundCorr_F, IndividualSubjectFirstCorrF{s});
end

krTGM_F = krTGM_F(:,:,krWinToUse, ctWinToUse);
krTGM_R = krTGM_R(:,:,krWinToUse, ctWinToUse);

end

function [clusterT, diffMatR,diffMatF] = scoreCluster(krTGM_R, krTGM_F, firstRoundCorr_R, firstRoundCorr_F, cluster, computationToPlot)

numF = size(krTGM_F, 1);
[numR, ~, numKRT, numCT_short] = size(krTGM_R);

numSamp = min([numF, numR]);

switch computationToPlot
    case 'Slope'
        
        diffMatR = nan(numSamp, numKRT, numCT_short);
        diffMatF = nan(numSamp, numKRT, numCT_short);
        
        for i = 1:numSamp
            for j = 1:numKRT
                for k = 1:numCT_short
                    weights = polyfit(1:4, squeeze(krTGM_R(i,:,j,k)),1);
                    diffMatR(i,j,k) = weights(1);
                    weights = polyfit(1:4, squeeze(krTGM_F(i,:,j,k)),1);
                    diffMatF(i,j,k) = weights(1);
                end
            end
        end
        
    case 'Mean(R1R2-R3R4)'
        diffMatR = squeeze(mean(krTGM_R(1:numSamp,1:2,:,:), 2) - mean(krTGM_R(1:numSamp,3:4,:,:), 2));
        diffMatF = squeeze(mean(krTGM_F(1:numSamp,1:2,:,:), 2) - mean(krTGM_F(1:numSamp,3:4,:,:), 2));
    case 'Mean(R3R4)'
        diffMatR = squeeze(mean(krTGM_R(1:numSamp,3:4,:,:), 2));
        diffMatF = squeeze(mean(krTGM_F(1:numSamp,3:4,:,:), 2));
    case 'R2-R3'
        diffMatR = squeeze(krTGM_R(1:numSamp,2,:,:) - krTGM_R(1:numSamp,3,:,:));
        diffMatF = squeeze(krTGM_F(1:numSamp,2,:,:) - krTGM_F(1:numSamp,3,:,:));
    case 'R3-R4'
        diffMatR = squeeze(krTGM_R(1:numSamp,3,:,:) - krTGM_R(1:numSamp,4,:,:));
        diffMatF = squeeze(krTGM_F(1:numSamp,3,:,:) - krTGM_F(1:numSamp,4,:,:));
    case 'R1-R4'
        diffMatR = squeeze(krTGM_R(1:numSamp,1,:,:) - krTGM_R(1:numSamp,4,:,:));
        diffMatF = squeeze(krTGM_F(1:numSamp,1,:,:) - krTGM_F(1:numSamp,4,:,:));
    case 'R2-R4'
        diffMatR = squeeze(krTGM_R(1:numSamp,2,:,:) - krTGM_R(1:numSamp,4,:,:));
        diffMatF = squeeze(krTGM_F(1:numSamp,2,:,:) - krTGM_F(1:numSamp,4,:,:));
    case 'R2-Mean(R3R4)'
        diffMatR = squeeze(krTGM_R(1:numSamp,2,:,:) - mean(krTGM_R(1:numSamp,3:4,:,:), 2));
        diffMatF = squeeze(krTGM_F(1:numSamp,2,:,:) - mean(krTGM_F(1:numSamp,3:4,:,:), 2));
    case 'R1-Mean(R2R3R4)'
        diffMatR = squeeze(krTGM_R(1:numSamp,1,:,:) - mean(krTGM_R(1:numSamp,2:4,:,:), 2));
        diffMatF = squeeze(krTGM_F(1:numSamp,1,:,:) - mean(krTGM_F(1:numSamp,2:4,:,:), 2));
    case 'FirstCorrect-R4'
        diffMatR = nan(numSamp, numKRT, numCT_short);
        diffMatF = nan(numSamp, numKRT, numCT_short);
        for i = 1:numSamp
            diffMatR(i,:,:) = squeeze(krTGM_R(i,firstRoundCorr_R(i),:,:) - krTGM_R(i,4,:,:));
            diffMatF(i,:,:) = squeeze(krTGM_F(i,firstRoundCorr_F(i),:,:) - krTGM_F(i,4,:,:));
        end
    otherwise
        error('Not Implemented')
end

[~, ~, ~, stats] = ttest(diffMatR, diffMatF, 'Tail', 'right');

clusterT = sum(stats.tstat(sub2ind([numKRT, numCT_short], cluster(:,1), cluster(:,2))));
end

function [clusters, clusterTestStat] = findClusters(data, options)
% Are we dealing with 1D series data, or 2D series?
numDim = length(size(data));

if numDim == 3
    [numSamp, numD1Points, numCond] = size(data);
elseif numDim == 4
    [numSamp, numD1Points, numD2Points, numCond] = size(data);
else
    error('Unsupported number of dimensions');
end

% Collect differences for calculating t values (for cluster
% calculation) and p values (for thresholding)
if numDim == 3
    pooledDifferences = diff(data, 1, 3);
    acrossSubjectTValues = nan(numD1Points, 1);
    acrossSubjectPValues = nan(numD1Points, 1);
else
    pooledDifferences = diff(data, 1, 4);
    acrossSubjectTValues = nan(numD1Points, numD2Points);
    acrossSubjectPValues = nan(numD1Points, numD2Points);
end


% Calculate t values and p values for each time or time-sensor pair across
% subjects (uncorrected)
for t = 1:numD1Points
    if numDim == 3
        [~, acrossSubjectPValues(t), ~, stats] = ttest(data(:,t, 1), data(:,t, 2),'Tail', 'right');
        acrossSubjectTValues(t) = stats.tstat;
    else
        for se = 1:numD2Points
            [~, acrossSubjectPValues(t, se), ~, stats] = ttest(data(:, t, se,1), data(:, t, se,2),'Tail', 'right');
            acrossSubjectTValues(t, se) = stats.tstat;
        end
    end
end

% Find the samples that survive the p value thresholding (linear indexing)
if numDim == 3
    thresholdedPoints = find(acrossSubjectPValues <= options.pValThresh);
else
    [thresholdedPointsI, thresholdedPointsJ] = find(acrossSubjectPValues <= options.pValThresh);
    [thresholdedPointsI, sortInd] = sort(thresholdedPointsI);
    thresholdedPointsJ = thresholdedPointsJ(sortInd);
    thresholdedPoints = [thresholdedPointsI, thresholdedPointsJ];
end

% Create clusters - first select contiguous D1 points
c = 1;
clustersByD1 = cell(1);
if ~isempty(thresholdedPoints)
    clustersByD1{1} = thresholdedPoints(1, :);
    for t = 2:size(thresholdedPoints,1)
        if (thresholdedPoints(t, 1) - thresholdedPoints(t-1, 1))  > 1
            clustersByD1 = cat(1, clustersByD1, thresholdedPoints(t, :));
            c = c+1;
        else
            clustersByD1{c} = cat(1, clustersByD1{c}, thresholdedPoints(t, :));
        end
    end
end

if numDim == 4
    [thresholdedPointsJ, sortInd] = sort(thresholdedPointsJ);
    thresholdedPointsI = thresholdedPointsI(sortInd);
    thresholdedPoints = [thresholdedPointsI, thresholdedPointsJ];
    c = 1;
    clustersByD2 = cell(1);
    if ~isempty(thresholdedPoints)
        clustersByD2{1} = thresholdedPoints(1, :);
        for t = 2:size(thresholdedPoints,1)
            if (thresholdedPoints(t, 2) - thresholdedPoints(t-1, 2))  > 1
                clustersByD2 = cat(1, clustersByD2, thresholdedPoints(t, :));
                c = c+1;
            else
                clustersByD2{c} = cat(1, clustersByD2{c}, thresholdedPoints(t, :));
            end
        end
    end
end

if numDim == 3 || isempty(clustersByD1{1})
    clusters = clustersByD1;
else
    clusters = cell(1);
    cc = 1;
    for c = 1:length(clustersByD1)
        newCluster = sortrows(clustersByD1{c}, [2 1]);
        clusters{cc,1} = newCluster(1,:);
        for t = 2:size(newCluster,1)
            
            piece1 = ((newCluster(t, 2) - newCluster(t-1, 2))  > 1) || ...
                (abs(newCluster(t, 1) - newCluster(t-1, 1))  > 1);
            piece2 = (abs(newCluster(t, 2) - newCluster(t-1, 2))  > 0) && ...
                (abs(newCluster(t, 1) - newCluster(t-1, 1))  > 0);
            
            if piece1  || piece2
                %                 if piece2
                %                     keyboard
                %                 end
                clusters = cat(1, clusters, newCluster(t, :));
                cc = cc+1;
            else
                clusters{cc,1} = cat(1, clusters{cc,1}, newCluster(t, :));
            end
        end
        cc = cc+1;
    end
    
    for c = 1:length(clustersByD2)
        newCluster = sortrows(clustersByD2{c}, [1 2]);
        clusters{cc,1} = newCluster(1,:);
        for t = 2:size(newCluster,1)
            if ((newCluster(t, 1) - newCluster(t-1, 1))  > 1 || ...
                    abs(newCluster(t, 2) - newCluster(t-1, 2))  > 1)  || ...
                    (abs(newCluster(t, 2) - newCluster(t-1, 2))  > 0 && ...
                    abs(newCluster(t, 1) - newCluster(t-1, 1)) > 0)
                clusters = cat(1, clusters, newCluster(t, :));
                cc = cc+1;
            else
                clusters{cc,1} = cat(1, clusters{cc,1}, newCluster(t, :));
            end
        end
        cc = cc+1;
    end
    origClusters = clusters;
    
    % Discard duplicates
    uniClusters = {};
    for ic = 1:length(clusters)
        isDuplicate = false;
        for jc = (ic+1):length(clusters)
            cluster1 = sortrows(clusters{ic},1);
            cluster2 = sortrows(clusters{jc},1);
            if all(size(cluster1) == size(cluster2))
                if all(cluster1 == cluster2)
                    isDuplicate = true;
                end
            end
        end
        if ~isDuplicate
            uniClusters = cat(1, uniClusters, clusters{ic});
        end
    end
    
    % Find if any clusters are touching each other
    clusters = {};
    absorbedClusters = [];
    for ic = 1:length(uniClusters)
        if ~ismember(ic, absorbedClusters)
            cluster1 = uniClusters{ic};
            for jc = (ic+1):length(uniClusters)
                cluster2 = uniClusters{jc};
                for ic1 = 1:size(cluster1, 1)
                    d1Touch = abs(cluster1(ic1, 1) - cluster2(:,1)) <=1;
                    d1_oneaway = abs(cluster1(ic1, 1) - cluster2(:,1)) ==1;
                    d2Touch = abs(cluster1(ic1,2) - cluster2(:,2)) <=1;
                    d2_oneaway = abs(cluster1(ic1,2) - cluster2(:,2)) ==1;
                    
                    isTouching = (d1Touch & d2Touch) & ~(d1_oneaway & d2_oneaway);
                    
                    if any(isTouching)
                        absorbedClusters = cat(1, absorbedClusters, jc);
                        cluster1 = cat(1, cluster1, cluster2);
                        
                        break
                    end
                end
            end
            clusters = cat(1, clusters, unique(cluster1, 'rows'));
        end
    end
end

% Discard clusters that are too small
numClusters = length(clusters);
thresholdedClusters = false(numClusters, 1);
for c = 1:numClusters
    if length(clusters{c}) >= options.minClusterSize
        thresholdedClusters(c) = true;
    end
end

clusters = clusters(thresholdedClusters);
numClusters = sum(thresholdedClusters);

% Compute the cluster test statistic
clusterTestStat = nan(numClusters, 1);
for c = 1:numClusters
    % Currently only one option
    switch options.clusterTestStatistic
        case 'summed_t_values'
            if numDim == 2
                clusterTestStat(c) = sum(acrossSubjectTValues(clusters{c}));
            else
                indForData = sub2ind([numD1Points, numD2Points], clusters{c}(:,1), clusters{c}(:,2));
                clusterTestStat(c) = sum(acrossSubjectTValues(indForData));
            end
    end
end

end

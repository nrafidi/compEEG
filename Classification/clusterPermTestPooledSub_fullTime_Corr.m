function [clusters, monteCarloPvals, permutationClusters, ...
    permutationHist, sizePermClusters] = ...
    clusterPermTestPooledSub_fullTime_Corr(dataX, dataY, varargin)
% CLUSTERPERMTEST: runs the cluster permutation test described in Maris &
% Oostenveld 2007 for two conditions in a within-subjects MEG study
%
%
% Optional Inputs:
% options: a struct with any of the following fields:
%   pValThresh - uncorrected pvalue used to threshold initial pvalues to get
%       clusters (default 0.05)
%   minClusterSize - clusters greater than this size will be considered.
%       Note that this number should be very different depending on whether
%       your samples are time points or sensor-time points. (default 10)
%   clusterTestStatistic - summary test statistic for a cluster. Only one
%       supported option: 'summed_t_values' the sum of the t values of the
%       samples within the cluster (default)
%   maxPermutations - maximum number of condition assignment permutations to
%       try of the 2^numSubjects possible (default 1000)
%
% Outputs:
% significantClusters: a cell array where each entry is an array containing
%   the indices of a significant cluster, e.g. the locations in the
%   corresponding time vector of a significant cluster. Will be empty if
%   none are found. For space-time data, the first column of the cell array
%   entry will be the dim1 index and the second column will be the matching
%   dim2 index.
% monteCarloPvals: the Monte Carlo p values for each of the significant
%   clusters. Note that p < 0.025 is the threshold for significance.

% Checking Arguments
if nargin < 3
    options = struct();
elseif nargin == 3
    options = varargin{1};
else
    error('Too many input arguments');
end

% Are we dealing with 1D series data, or 2D series?
numDim = length(size(dataY));

if numDim == 2
    [numSamp, numD1Points] = size(dataY);
elseif numDim == 3
    [numSamp, numD1Points, numD2Points] = size(dataY);
else
    error('Unsupported number of dimensions');
end

% Setting options to defaults if they are not provided
if ~isfield(options, 'pValThresh')
    options.pValThresh = 0.05;
end
if ~isfield(options, 'minClusterSize')
    options.minClusterSize = 10;
end
if ~isfield(options, 'clusterTestStatistic')
    options.clusterTestStatistic = 'summed_t_values';
end
if ~isfield(options, 'maxPermutations')
    options.maxPermutations = 1000;
end

[clusters, clusterTestStat] = findClusters(dataX, dataY, options);
numClusters = size(clusters, 1);

% Get the permutations of condition assignments
numPermutations = options.maxPermutations;
permutationClusters = cell(numPermutations, 1);
permutationTestStat = [];
numPermClusters = [];
sizePermClusters = [];
permutationHist = [];
percFlipped = nan(numPermutations, 1);

% Compute the cluster statistic on the max cluster data for each of the
% permuted condition assignments
for p = 1:numPermutations
    disp(p)
    permuted = randperm(numSamp);
    
    if numDim == 2
        permutedData = dataY(permuted, :);
    else
        permutedData = dataY(permuted, :, :);
    end   
    
    [permutationClusters{p}, permutationTestStat_single] = ...
        findClusters(dataX, permutedData, options);
    
    if isempty(permutationClusters{p})
        permutationTestStat_single = 0;
        sizePermClusters = cat(1, sizePermClusters, 0);
    else
        sizePermClusters = cat(1, sizePermClusters, max(cellfun(@(x) size(x, 1), permutationClusters{p})));
        
%         if max(permutationTestStat_single) >= 20 && ~haveSaved
%             meow = zeros(54, 21);
%             for c = 1:length(permutationClusters{p})
%                 linInd = sub2ind([54, 21], permutationClusters{p}{c}(:,1), permutationClusters{p}{c}(:,2));
%                 meow(linInd) = c;
%             end
%             f = figure;
%             imagesc(meow');
%             title(permutationTestStat_single)
%             
%             export_fig(f, '/Users/nrafidi/Documents/MATLAB/compEEG-data-rep/results/figures/sampleClust.png');
%             haveSaved = true;
%         end
        %         keyboard
    end
    
    [maxAbs, indMax] = max(abs(permutationTestStat_single));
    permutationHist = cat(1, permutationHist, permutationTestStat_single(indMax));
    permutationTestStat = cat(1, permutationTestStat, maxAbs);
    numPermClusters = cat(1, numPermClusters, length(permutationClusters{p}));
end

% For each cluster, its MC p value is the fraction of permutation values
% that are greater than its own. This is a two sided test, so we compare
% absolute values
% keyboard

monteCarloPvals = nan(numClusters, 2);
for c = 1:numClusters
    monteCarloPvals(c,1) = sum(abs(permutationTestStat) > abs(clusterTestStat(c)))/length(permutationTestStat);
    monteCarloPvals(c,2) = sum(sizePermClusters > size(clusters{c},1))/length(sizePermClusters);
end
disp(nnz(numPermClusters))
disp(sum(numPermClusters > numClusters));

end

function [clusters, clusterTestStat] = findClusters(dataX, dataY, options)
% Are we dealing with 1D series data, or 2D series?
numDim = length(size(dataY));

if numDim == 2
    [numSamp, numD1Points] = size(dataY);
elseif numDim == 3
    [numSamp, numD1Points, numD2Points] = size(dataY);
else
    error('Unsupported number of dimensions');
end

% Collect differences for calculating t values (for cluster
% calculation) and p values (for thresholding)
if numDim == 2
    pooledCorr = nan(numD1Points, 1);
    acrossSubjectPValues = nan(numD1Points, 1);
else
    pooledCorr = nan(numD1Points, numD2Points);
    acrossSubjectPValues = nan(numD1Points, numD2Points);
end


% Calculate t values and p values for each time or time-sensor pair across
% subjects (uncorrected)
for t = 1:numD1Points
    if numDim == 2
        [pooledCorr(t), acrossSubjectPValues(t)] = corr(dataX, dataY(:,t),'tail', 'right');
    else
        for se = 1:numD2Points
            [pooledCorr(t, se), acrossSubjectPValues(t, se)] = corr(dataX, dataY(:,t, se),'tail', 'right');
        end
    end
end

% Find the samples that survive the p value thresholding (linear indexing)
if numDim == 2
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

if numDim == 3
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

if numDim == 2 || isempty(clustersByD1{1})
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
    
%     meow = zeros(49, 30);
%     for c = 1:length(clusters)
%         linInd = sub2ind([49, 30], clusters{c}(:,1), clusters{c}(:,2));
%         if all(meow(linInd) == 0)
%             meow(linInd) = c;
%         else
%             meow(linInd) = 20;
%         end
%     end
%     figure
%     imagesc(meow);
%     colorbar
%     keyboard;
    
    
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
    %
    %     meow = zeros(54, 21);
    %     for c = 1:length(uniClusters)
    %         linInd = sub2ind([54, 21], uniClusters{c}(:,1), uniClusters{c}(:,2));
    %         if all(meow(linInd) == 0)
    %             meow(linInd) = c;
    %         else
    %             meow(linInd) = 20;
    %         end
    %     end
    %     figure
    %     imagesc(meow);
    %     colorbar
    %     keyboard;
    
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
                        
                        %                         meow = zeros(54, 21);
                        %                         linInd = sub2ind([54, 21], cluster1(:,1), cluster1(:,2));
                        %                         meow(linInd) = 1;
                        %                         linInd = sub2ind([54, 21], cluster2(:,1), cluster2(:,2));
                        %                         meow(linInd) = 2;
                        %
                        %                         figure
                        %                         imagesc(meow);
                        %
                        %                         keyboard
                        
                        cluster1 = cat(1, cluster1, cluster2);
                        
                        break
                    end
                end
            end
            clusters = cat(1, clusters, unique(cluster1, 'rows'));
        end
    end
end
%
% if ~isempty(clusters{1})
%     meow = zeros(54, 21);
%     for c = 1:length(clusters)
%         linInd = sub2ind([54, 21], clusters{c}(:,1), clusters{c}(:,2));
%         if all(meow(linInd) == 0)
%             meow(linInd) = c;
%         else
%             meow(linInd) = 20;
%         end
%     end
%     figure
%     imagesc(meow);
%     keyboard
% end

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
    clusterTestStat(c) = length(clusters{c});
end

end
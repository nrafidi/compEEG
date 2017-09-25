function [significantClusters, monteCarloPvals] = clusterPermTestMultiSamp(IndividualSubjectData, varargin)
% CLUSTERPERMTEST: runs the cluster permutation test described in Maris &
% Oostenveld 2007 for two conditions in a within-subjects MEG study
%
% Inputs:
% IndividualSubjectData: a cell array where each entry is the data for one
%   subject in one of the following two forms:
%       numSamp x dim1 x condition
%       numSamp x dim1 x dim2 x condition
%
% dim1 and dim2 could be time, space, whatever
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
if nargin < 2
    options = struct();
elseif nargin == 2
    options = varargin{1};
else
    error('Too many input arguments');
end

% Getting basic parameters about the input data
numSubjects = length(IndividualSubjectData);

% Are we dealing with 1D series data, or 2D series?
numDim = length(size(IndividualSubjectData{1}));

if numDim == 3
    [~, numD1Points, numCond] = size(IndividualSubjectData{1});
elseif numDim == 4
    [~, numD1Points, numD2Points, numCond] = size(IndividualSubjectData{1});
else
    error('Unsupported number of dimensions');
end

% Checking the number of conditions (right now can only handle two)
if numCond ~= 2
    error('Unsupported number of conditions');
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

% Collect differences within subject for calculating t values (for cluster
% calculation) and p values (for thresholding)
if numDim == 3
    withinSubjectDifferences = [];
    acrossSubjectTValues = nan(numD1Points, 1);
    acrossSubjectPValues = nan(numD1Points, 1);
else
    withinSubjectDifferences = [];
    acrossSubjectTValues = nan(numD1Points, numD2Points);
    acrossSubjectPValues = nan(numD1Points, numD2Points);
end

for s = 1:numSubjects
    if numDim == 3
        withinSubjectDifferences = cat(1, withinSubjectDifferences, diff(IndividualSubjectData{s}, 1, 3));
    else
        withinSubjectDifferences = cat(1, withinSubjectDifferences, diff(IndividualSubjectData{s}, 1, 4));
    end
end

% Calculate t values and p values for each time or time-sensor pair across
% subjects (uncorrected)
for t = 1:numD1Points
    if numDim == 3
        [~, acrossSubjectPValues(t), ~, stats] = ttest(withinSubjectDifferences(:,t));
        acrossSubjectTValues(t) = stats.tstat;
    else
        for se = 1:numD2Points
            [~, acrossSubjectPValues(t, se), ~, stats] = ttest(withinSubjectDifferences(:, t, se));
            acrossSubjectTValues(t, se) = stats.tstat;
        end
    end
end

% Find the samples that survive the p value thresholding (linear indexing)
thresholdedPoints = find(acrossSubjectPValues < options.pValThresh);
% If using 2D data, convert these indices to a matrix so that we can
% cluster D1 and D2 separately
if numDim == 4
    [I, J] = ind2sub([numD1Points, numD2Points], thresholdedPoints);
    [I, ind] = sort(I);
    thresholdedPoints = [I J(ind)];
end

% Create clusters - first select contiguous D1 points
c = 1;
clustersByD1 = cell(1);
clustersByD1{1} = thresholdedPoints(1, :);
for t = 2:length(thresholdedPoints)
    if (thresholdedPoints(t, 1) - thresholdedPoints(t-1, 1))  > 1
        clustersByD1 = cat(1, clustersByD1, thresholdedPoints(t, :));
        c = c+1;
    else
        clustersByD1{c} = cat(1, clustersByD1{c}, thresholdedPoints(t, :));
    end
end
clusters = clustersByD1;
% If using sensor data, divide the above clusters into ones that are also
% spatially proximal
% if numDim == 3
%     clusters = clustersByD1;
% else
%     numD1Clusters = length(clustersByD1);
%     clustersByD2 = cell(1);
%     for c = 1:numD1Clusters
%         clusterD2 = clustersByD1{c}(:, 2);
%         numClusterD2 = length(clusterD2);
%         
%         % Compute pairwise distances between the sensors in this cluster
%         clusterSensorCoord = options.sensorCoord(clusterD2, :);
%         distanceSquare = squareform(pdist(clusterSensorCoord));
%         distanceSquare = distanceSquare < options.clusterSpaceConstraint;
%         
%         
%         % Get the indices of the clusters within the appropriate proximity
%         clusterSensorGroups = {find(distanceSquare(:,1))};
%         
%         % Add the appropriate clusters, being careful not to add duplicates
%         for cc = 2:numClusterD2
%             clusterToAdd = find(distanceSquare(:,cc));
%             canAdd = true;
%             for ccc = 1:length(clusterSensorGroups)
%                 if length(clusterToAdd) == length(clusterSensorGroups{ccc})
%                     if sum(clusterToAdd == clusterSensorGroups{ccc}) == length(clusterToAdd)
%                         canAdd = false;
%                     end
%                 end
%             end
%             if canAdd
%                 clusterSensorGroups = cat(1, clusterSensorGroups, clusterToAdd);
%             end
%         end
%         
%         for cc = 2:length(clusterSensorGroups)
%             clustersByD2 = cat(1, clustersByD2, clustersByD1{c}(clusterSensorGroups{cc}, :));
%         end
%     end
%     
%     % Parse again by time to make sure we still have contiguous time
%     % windows
%     c = 1;
%     clusters = cell(1);
%     for cc = 1:length(clustersByD2)
%         for t = 1:length(clustersByD2{cc})
%             if t > 1
%                 if (clustersByD2{cc}(t, 1) - clustersByD2{cc}(t-1, 1))  > 1
%                     clusters = cat(1, clusters, clustersByD2{cc}(t, :));
%                     c = c+1;
%                 else
%                     clusters{c} = cat(1, clusters{c}, clustersByD2{cc}(t, :));
%                 end
%             else
%                 clusters{c} = cat(1, clusters{c}, clustersByD2{cc}(t, :));
%             end
%         end
%     end
% end

% Discard clusters that are too small
numClusters = length(clusters);
thresholdedClusters = false(numClusters, 1);
for c = 1:numClusters
    if length(clusters{c}) > options.minClusterSize
        thresholdedClusters(c) = true;
    end
end
if ~any(thresholdedClusters)
    error('No clusters survived thresholding');
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

% Find the max valued cluster
[~, maxCluster] = max(clusterTestStat);
maxClusterLength = size(clusters{maxCluster}, 1);

% Aggregate data corresponding to the max valued cluster
IndividualSubjectClusters = cell(numSubjects, 1);
for s = 1:numSubjects
    if numDim == 3
        IndividualSubjectClusters{s} = IndividualSubjectData{s}(clusters{maxCluster}, :);
    else
        indForData = sub2ind([numD1Points, numD2Points], clusters{maxCluster}(:,1), clusters{maxCluster}(:,2));
        indexedData = [];
        for n = 1:numCond
            dataForCond = squeeze(IndividualSubjectData{s}(:, :, 1));
            indexedData = cat(2, indexedData, dataForCond(indForData));%logicIndForData));
        end
        IndividualSubjectClusters{s} = indexedData;
    end
end

% Get the permutations of condition assignments
numPermutations = numCond^numSubjects;
numPermutations = min(numPermutations, options.maxPermutations);
IndividualSubjectPermutationDiffs = cell(numPermutations, 1);
acrossSubjectPermutationTValues = nan(numPermutations, maxClusterLength);

% Compute the cluster statistic on the max cluster data for each of the
% permuted condition assignments
for p = 1:numPermutations
    IndividualSubjectPermutationDiffs{p} = nan(numSubjects, maxClusterLength);
    permuted = de2bi((p-1), numSubjects);
    permuted = [permuted + 1; 2 - permuted];
    
    % For each permutation, calculate the cluster statistic of the max
    % cluster
    for s = 1:numSubjects
        IndividualSubjectPermutationDiffs{p}(s,:) = diff(IndividualSubjectClusters{s}(:, permuted(:,s))');
    end
    for t = 1:maxClusterLength
        [~, ~, ~, stats] = ttest(IndividualSubjectPermutationDiffs{p}(:,t));
        acrossSubjectPermutationTValues(p, t) = stats.tstat;
    end
end

% The histogram of permutation statistic values
permutationHist = sum(acrossSubjectPermutationTValues, 2);

% For each cluster, its MC p value is the fraction of permutation values
% that are greater than its own. This is a two sided test, so we compare
% absolute values
monteCarloPvals = nan(numClusters, 1);
for c = 1:numClusters
    monteCarloPvals(c) = sum(abs(permutationHist) >= abs(clusterTestStat(c)))/numPermutations;
end

% Return the signficant clusters and their p values
significantClusters = clusters(monteCarloPvals < 0.025);
monteCarloPvals = monteCarloPvals(monteCarloPvals < 0.025);

end
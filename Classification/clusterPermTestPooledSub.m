function [significantClusters, monteCarloPvals] = clusterPermTestPooledSub(pooledData, varargin)
% CLUSTERPERMTEST: runs the cluster permutation test described in Maris &
% Oostenveld 2007 for two conditions in a within-subjects MEG study
%
% Inputs:
% pooledData: a double array in one of the following two forms:
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

% Are we dealing with 1D series data, or 2D series?
numDim = length(size(pooledData));

if numDim == 3
    [numSamp, numD1Points, numCond] = size(pooledData);
elseif numDim == 4
    [numSamp, numD1Points, numD2Points, numCond] = size(pooledData);
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

% Collect differences for calculating t values (for cluster
% calculation) and p values (for thresholding)
if numDim == 3
    pooledDifferences = diff(pooledData, 1, 3);
    acrossSubjectTValues = nan(numD1Points, 1);
    acrossSubjectPValues = nan(numD1Points, 1);
else
    pooledDifferences = diff(pooledData, 1, 4);
    acrossSubjectTValues = nan(numD1Points, numD2Points);
    acrossSubjectPValues = nan(numD1Points, numD2Points);
end


% Calculate t values and p values for each time or time-sensor pair across
% subjects (uncorrected)
for t = 1:numD1Points
    if numDim == 3
        [~, acrossSubjectPValues(t), ~, stats] = ttest(pooledDifferences(:,t));
        acrossSubjectTValues(t) = stats.tstat;
    else
        for se = 1:numD2Points
            [~, acrossSubjectPValues(t, se), ~, stats] = ttest(pooledDifferences(:, t, se));
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
if numDim == 3
    maxClusterData = pooledData(:, clusters{maxCluster}, :);
else
    indForData = sub2ind([numD1Points, numD2Points], clusters{maxCluster}(:,1), clusters{maxCluster}(:,2));
    maxClusterData = [];
    for n = 1:numCond
        dataForCond = [];
        for s = 1:numSamp
            sampData = squeeze(pooledData(s, :, :, n));
            dataForCond = cat(1, dataForCond, sampData(indForData)');
        end
        maxClusterData = cat(3, maxClusterData, dataForCond);%logicIndForData));
    end
end

% Get the permutations of condition assignments
numPermutations = options.maxPermutations;
permutationDiffs = cell(numPermutations, 1);
acrossSubjectPermutationTValues = nan(numPermutations, maxClusterLength);

% Compute the cluster statistic on the max cluster data for each of the
% permuted condition assignments
for p = 1:numPermutations
    permutationDiffs{p} = nan(numSamp, maxClusterLength);
    permuted = randi([0,1], 1, numSamp);
    permuted = [permuted + 1; 2 - permuted];
    
    % For each permutation, calculate the cluster statistic of the max
    % cluster
    for s = 1:numSamp
        permutationDiffs{p}(s,:) = diff(maxClusterData(s, :, permuted(:,s)), 1, 3);
    end
    for t = 1:maxClusterLength
        [~, ~, ~, stats] = ttest(permutationDiffs{p}(:,t));
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
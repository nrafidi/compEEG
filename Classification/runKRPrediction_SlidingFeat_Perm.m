function [ itemTraj ] = runKRPrediction_SlidingFeat_Perm(sub, krData, krLabels, compWin, krWin)
% Applies a classifier trained on a subject's competition data to their KR
% data

addpath ./logisticRegression/
numChan = 64;
dataRoot = ['/Users/nrafidi/Documents/MATLAB/compEEG-data/preproc-final/' sub '/'];

fCPrefix = [dataRoot 'CompEEG_'];
fResPrefix = ['/Users/nrafidi/Documents/MATLAB/compEEG-data/results/' sub '/'];
if ~exist(fResPrefix, 'dir')
    mkdir(fResPrefix);
end
fSuffix = '_Vis_BP2-200_N60_Ref_Epochs_Base_ICA1-2_Features_Overlap_Time.mat';

rng('shuffle');

compFname = [fCPrefix sub fSuffix];
load(compFname);
numSamp = size(featData, 1);
featData = double(featData);
[compData, mu, sigma] = zscore(featData);
compLabels = labels(randperm(numSamp),1); %#ok<*NODEF>

compWinToUse = (compWin*numChan + 1):((compWin+1)*numChan);

B = logReg(compData(:,compWinToUse), compLabels, [1 10], false);

% old_krLabels = krLabels;
krData = krData(krLabels(:,1) == 1, :);
krData = (krData - repmat(mu, size(krData, 1), 1))./repmat(sigma, size(krData, 1), 1);
krLabels = krLabels(krLabels(:,1) ==1, 2);


winToTry = find(winTime >= krWin(1) & winTime <= krWin(2));
numWinToTry = length(winToTry);
uniqueItems = unique(krLabels);
itemTraj = nan(length(uniqueItems), 4, numWinToTry);
for i = 1:length(uniqueItems)
    
    items = krLabels == (i+1);
    for w = 1:numWinToTry
        win = winToTry(w);
        timeWindow = ((win-1)*numChan + 1):(win*numChan);
        P = [krData(items,timeWindow) ones(sum(items), 1)]*B;
        probs = exp(P)./(1+exp(P));
        if length(probs) == 4
            itemTraj(i, :, w) = probs';
        end
    end
    
end

end


function [ itemTraj ] = runKRPrediction_SlidingFeat(sub, krData, krLabels, compWin, isRepExp)
% Applies a classifier trained on a subject's competition data to their KR
% data

addpath ./logisticRegression/
numChan = 64;

if isRepExp
    dataRoot = ['/Users/nrafidi/Documents/MATLAB/compEEG-data-rep/preproc-final/' sub '/'];
    fResPrefix = ['/Users/nrafidi/Documents/MATLAB/compEEG-data-rep/results/' sub '/'];
else
    dataRoot = ['/Users/nrafidi/Documents/MATLAB/compEEG-data/preproc-final/' sub '/'];
    fResPrefix = ['/Users/nrafidi/Documents/MATLAB/compEEG-data/results/' sub '/'];
end

fBname = [fResPrefix 'CompEEG_B_cWin' num2str(compWin) '.mat'];
fCPrefix = [dataRoot 'CompEEG_'];

if ~exist(fResPrefix, 'dir')
    mkdir(fResPrefix);
end
fSuffixString = '%s_BP2-%d_N60_Ref_Epochs_Base_ICA1-2_Features_Overlap_Time.mat';
if strcmp(sub, 'NN')
    fSuffix = sprintf(fSuffixString, '_Vis', 127);
elseif strcmp(sub, 'HHH')
    fSuffix = sprintf(fSuffixString, '', 200);
else
    fSuffix = sprintf(fSuffixString, '_Vis', 200);
end
rng('shuffle');

compFname = [fCPrefix sub fSuffix];

if exist(fBname, 'file')
    load(fBname);
else
    
    load(compFname);
    
    featData = double(featData);
    [compData, mu, sigma] = zscore(featData);
    compLabels = labels(:,1); %#ok<*NODEF>
    
    compWinToUse = (compWin*numChan + 1):((compWin+1)*numChan);
    
    B = logReg(compData(:,compWinToUse), compLabels, [1 10], false);
    
    save(fBname, 'B', 'mu', 'sigma');
end
% old_krLabels = krLabels;
krData = krData(krLabels(:,1) == 1, :);
krLabels = krLabels(krLabels(:,1) ==1, 2);


numWinToTry = size(krData,2)/numChan;
uniqueItems = unique(krLabels);
if any(uniqueItems == 255)
    uniqueItems = uniqueItems(uniqueItems ~= 255);
end
itemTraj = nan(length(uniqueItems), 4, numWinToTry);
for i = 1:length(uniqueItems)
    
    items = krLabels == (i+1);
    for w = 1:numWinToTry
        timeWindow = ((w-1)*numChan + 1):(w*numChan);
        krData(items, timeWindow) = (krData(items,timeWindow) - repmat(mu(timeWindow), sum(items), 1))./repmat(sigma(timeWindow), sum(items), 1);
        P = [krData(items,timeWindow) ones(sum(items), 1)]*B;
        probs = exp(P)./(1+exp(P));
        if length(probs) == 4
            itemTraj(i, :, w) = probs';
        end
    end
    
end

end


function [ itemTraj ] = runKRPrediction_SlidingFeat_theta(sub, krData, krLabels, compWin, isRepExp)
% Applies a classifier trained on a subject's competition data to their KR
% data

addpath ./logisticRegression/
numChan = 64;

rootDir = '/Users/nrafidi/Documents/MATLAB/compEEG-data%s/';

if isRepExp
    resultDir = [sprintf(rootDir, '-rep') 'results/'];
    dataDir = [sprintf(rootDir, '-rep') 'preproc-final/' sub '/'];
else
    resultDir = [sprintf(rootDir, '') 'results/'];
    dataDir = [sprintf(rootDir, '') 'preproc-final/' sub '/'];
end

fSuffixString = '%s_BP2-%d_N60_Ref_Hilbert-theta_Epochs_Features_Overlap_Time.mat';
if strcmp(sub, 'NN')
    fSuffix = sprintf(fSuffixString, '_Vis', 127);
% elseif strcmp(sub, 'HHH')
%     fSuffix = sprintf(fSuffixString, '', 200);
else
    fSuffix = sprintf(fSuffixString, '_Vis', 200);
end

compFname = [dataDir 'CompEEG_' sub fSuffix];
load(compFname);    
featData = double(featData);
[~, mu, sigma] = zscore(featData);
% disp('meow');

loadFileString = 'CompEEG_%s_CV_Slide_Models_PLOS_theta.mat';

rng('shuffle');

load([resultDir sprintf(loadFileString, sub)]);
B = models{compWin, 1};

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


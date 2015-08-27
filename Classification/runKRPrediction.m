function [ itemTraj ] = runKRPrediction(sub)
% Applies a classifier trained on a subject's competition data to their KR
% data

addpath ./logisticRegression/

dataRoot = ['/Users/nrafidi/Documents/MATLAB/compEEG-data/preproc-final/' sub '/'];

fCPrefix = [dataRoot 'CompEEG_'];
fResPrefix = ['/Users/nrafidi/Documents/MATLAB/compEEG-data/results/' sub '/'];
if ~exist(fResPrefix, 'dir')
    mkdir(fResPrefix);
end
fSuffix = '_Features.mat';

rng('shuffle');

fDir = dir([fCPrefix sub '_*' fSuffix]);
load([dataRoot fDir(end).name]);

[compData, mu, sigma] = zscore(featData);
compLabels = labels(:,1); %#ok<*NODEF>

B = logReg(compData, compLabels, [1 10], false);

fDir = dir([fCPrefix '_KR*' fSuffix]);
load([dataRoot fDir(end).name]);
krData = featData;
krLabels = labels(:,2);

krData = krData(labels(:,1) == 1, :);
krData = (krData - repmat(mu, size(krData, 1), 1))./repmat(sigma, size(krData, 1), 1);
krLabels = krLabels(labels(:,1) ==1, :);

uniqueItems = unique(krLabels);
itemTraj = nan(length(uniqueItems), 4);
for i = 1:length(uniqueItems)
    
    
    items = krLabels == (i+1);
    
    P = [krData(items,:) ones(sum(items), 1)]*B;
%     probs = exp(P)./(1+exp(P));
    probs = P;
    if length(probs) == 4
        itemTraj(i, :) = probs';
    else
%         keyboard;
    end
end

figure;
plot(1:4, itemTraj);

save([fResPrefix 'CompEEG__KR_' sub '_itemTraj_logit.mat'], 'itemTraj');

end


% Applies a classifier trained on a subject's competition data to their KR
% data
addpath ./GNB/

subjects = 'CE';
numSub = length(subjects);
fCPrefix = '/Users/nrafidi/Documents/MATLAB/compEEG-data/preproc-final/CompEEG_';
fKPrefix = [fCPrefix 'KR_'];
fResPrefix = '/Users/nrafidi/Documents/MATLAB/compEEG-data/results/CompEEG_KR_';
fSuffix = '_Preprocessed.mat';

rng('shuffle');

for s = 1:length(subjects)
    
    sub = subjects(s);
    load([fCPrefix sub fSuffix]);
    [compData, mu, sigma] = zscore(squeeze(mean(data, 3)));
    compLabels = labels(:,1) + 1;
    
    B = nbayes_train(compData, compLabels, 1);
    
    load([fKPrefix sub fSuffix]);
    krData = squeeze(mean(data, 3));
    krLabels = labels(:,2);
    
    krData = krData(labels(:,1) == 1, :);
    krData = (krData - repmat(mu, size(krData, 1), 1))./repmat(sigma, size(krData, 1), 1);
    krLabels = krLabels(labels(:,1) ==1, :);
    
    uniqueItems = unique(krLabels);
    itemTraj = nan(length(uniqueItems), 4);
    for i = 1:length(uniqueItems)
        

        items = krLabels == (i+1);
        
        P = nbayes_apply(krData(items,:), B);
        probs = exp(P);
        probs = probs./repmat(sum(probs,2), 1, 2);
        
        if length(probs) == 4
            itemTraj(i, :) = probs(:,2)';
        end
    end
    
    figure;
    plot(1:4, itemTraj);
    
    save([fResPrefix sub '_itemTraj.mat'], 'itemTraj');
end



% Runs cross-validated classification within subject for the competition
% data
addpath ./GNB/

subjects = 'C';
numSub = length(subjects);
% numFolds = 10;
numFeats = 48;
fPrefix = '/Users/nrafidi/Documents/MATLAB/compEEG-data/preproc-final/CompEEG_';
fSuffix = '_Preprocessed.mat';

subAccs = nan(numSub, numFolds);
subModels = cell(numSub, numFolds);
subProbs = cell(numSub, numFolds);
subConf = cell(numSub, numFolds);
subFolds = cell(numSub, 1);

rng('shuffle');

for s = 1:length(subjects)
    
    sub = subjects(s);
    load([fPrefix sub fSuffix]);
    
       
    Y = labels(:,1);
    
    N = size(data, 1);
    numFolds = 10;
    
    
    X = squeeze(mean(data, 3));
%     X = reshape(data(:,:,time > 0.2), N, []);
    folds = crossvalind('Kfold', N, numFolds);
    subFolds{s} = folds;
    for f = 1:numFolds
        testInds = (folds == f);
        
        [trainData, mu, sigma] = zscore(X(~testInds, :));
        testData = (X(testInds,:) - repmat(mu, sum(testInds), 1))./repmat(sigma,sum(testInds), 1);
            
%         [B, stats] = mnrfit(trainData, Y(~testInds) + 1);
%         P = mnrval(B, testData);
%         Yhat = double(P(:,2) > 0.5);

        B = nbayes_train(trainData, Y(~testInds)+1, 1);
        B.sortedFeats = sortGNBfeaturesByMuDifference(B);
        selectedFeats = B.sortedFeats(1:numFeats);
        
        P = nbayes_apply(testData(:, selectedFeats), B);
        [~, predictedLabelIndices] = max(P, [], 2);
        Yhat = B.labelVocab(predictedLabelIndices) - 1;
% 
%         [W, stats] = mnrfit(trainData(:,selectedFeats), Y(~testInds) + 1);
%         P = mnrval(W, testData(:,selectedFeats));
%         Yhat = double(P(:,2) > 0.5);

        subAccs(s, f) = sum(Yhat == Y(testInds))/length(Yhat);
        subModels{s, f} = B;
        probs = exp(P);
        subProbs{s, f} = probs./repmat(sum(probs,2), 1, 2);
        subConf{s, f} = subProbs{s,f}(Yhat ~= Y(testInds), 1);
    end
    
end

mean(subAccs, 2)

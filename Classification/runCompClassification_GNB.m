% Runs cross-validated classification within subject for the competition
% data
addpath ./GNB/

subjects = 'C';
numSub = length(subjects);
numFolds = 10;
% numFeats = 1:10:500;
% doZ = 0;
fPrefix = '/Users/nrafidi/Documents/MATLAB/compEEG-data/preproc-final/CompEEG_';
fSuffix = '_Features.mat';

subAccs = nan(numSub, numFolds);
subModels = cell(numSub, numFolds);
subProbs = cell(numSub, numFolds);
subConf = cell(numSub, numFolds);
subFolds = cell(numSub, 1);

rng('shuffle');

testVals = nan(20, 50, 2);

for reps = 1:20
    nn = 1;
    for numFeats = 1:10:500
        for doZ = 0:1
            for s = 1:length(subjects)
                
                sub = subjects(s);
                load([fPrefix sub fSuffix]);
                
                
                Y = labels(:,1);
                N = size(featData, 1);
                X = featData;
                
                folds = crossvalind('Kfold', N, numFolds);
                subFolds{s} = folds;
                for f = 1:numFolds
                    testInds = (folds == f);
                    
                    if doZ
                        [trainData, mu, sigma] = zscore(X(~testInds, :));
                        testData = (X(testInds,:) - repmat(mu, sum(testInds), 1))./repmat(sigma,sum(testInds), 1);
                    else
                        trainData = X(~testInds, :);
                        testData = X(testInds,:);
                    end
                    
                    B = nbayes_train(trainData, Y(~testInds)+1, 1);
                    B.sortedFeats = sortGNBfeaturesByMuDifference(B);
                    selectedFeats = B.sortedFeats(1:numFeats);
                    
                    P = nbayes_apply(testData(:, selectedFeats), B);
                    [~, predictedLabelIndices] = max(P, [], 2);
                    Yhat = B.labelVocab(predictedLabelIndices) - 1;
                    
                    subAccs(s, f) = sum(Yhat == Y(testInds))/length(Yhat);
                    subModels{s, f} = B;
                    probs = exp(P);
                    subProbs{s, f} = probs./repmat(sum(probs,2), 1, 2);
                    subConf{s, f} = subProbs{s,f}(Yhat ~= Y(testInds), 1);
                end
                
            end
            
            testVals(reps, nn, doZ+1) = mean(subAccs, 2);
        end
        nn = nn+1;
    end
end

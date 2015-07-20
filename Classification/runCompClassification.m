% Runs cross-validated classification within subject for the competition
% data
addpath ./logisticRegression/

subjects = 'CDEHILOP';
numSub = length(subjects);
numFolds = 5;
% numFeats = 1:10:500;
doZ = 1;
fPrefix = '/Users/nrafidi/Documents/MATLAB/compEEG-data/preproc-final/CompEEG_';
fSuffix = '_Features_Less.mat';

subAccs = nan(numSub, numFolds);
subModels = cell(numSub, numFolds);
subProbs = cell(numSub, numFolds);
subConf = cell(numSub, numFolds);
subFolds = cell(numSub, 1);

rng('shuffle');

testVals = nan(2, 1);

% for reps = 1:20
% for doZ = 0:1
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
            trainData = X(~testInds, :); %#ok<*UNRCH>
            testData = X(testInds,:);
        end
        tic
        if f > 1
            [B, lambda] = logReg(trainData, Y(~testInds), [1 10], false, B);
        else
            [B, lambda] = logReg(trainData, Y(~testInds), [1 10], false);
        end
        toc
        disp(lambda);
        P = [testData ones(sum(testInds), 1)]*B;
        P = exp(P)./(1 + exp(P));
        Yhat = double(P > 0.5);
        
        subAccs(s, f) = sum(Yhat == Y(testInds))/length(Yhat);
        disp(subAccs(s, f));
        subModels{s, f} = B;
        subProbs{s, f} = P;
        subConf{s, f} = subProbs{s,f}(Yhat ~= Y(testInds), 1);
    end
    
    subAcc = subAccs(s,:);
    subModel = subModels(s,:);
    subFold = subFolds{s};
    
    save(['/Users/nrafidi/Documents/MATLAB/compEEG-data/results/CompEEG_CV_' sub '_Less.mat'],...
        'subAcc', 'subModel', 'subFold');
    
end

disp(mean(subAccs, 2));

% save(['/Users/nrafidi/Documents/MATLAB/compEEG-data/results/CompEEG_CV_' subjects '.mat']);

%     testVals(doZ+1) = mean(subAccs, 2);
% end
% end

% addpath logisticRegression/
addpath GNB/
dataRoot = '/Users/nrafidi/Documents/MATLAB/compEEG-data/preproc-final/';
resultRoot = '/Users/nrafidi/Documents/MATLAB/compEEG-data/results/';

numT = 44;
numDim = 128;
numFolds = 5;
numFeats = 50;

load(sprintf('%sallSubData_KR_5fold_withPassLabels.mat', dataRoot));

accPerFold = nan(numT, numFolds);
% weightVecs = cell(numT, numFolds);
GNBmodels = cell(numT, numFolds);
for fold = 1:numFolds
    numSub = size(allSubTrainData, 1);
    
%     testLabels = cell2mat(allSubTestPassLabels(:,fold));
%     trainLabels = cell2mat(allSubTrainPassLabels(:,fold));
%     testLabels = cell2mat(allSubTestLabels(:,fold));
%     trainLabels = cell2mat(allSubTrainLabels(:,fold));

    trainLabels = [];
    testLabels = [];
    for s = 1:numSub
        trainLabels = cat(1, trainLabels, ...
            s*ones(length(allSubTrainLabels{s, fold}), 1));
        testLabels = cat(1, testLabels, ...
            s*ones(length(allSubTestLabels{s, fold}), 1));
    end
    
    numTestSamples = length(testLabels);
    h = waitbar(0, sprintf('Fold %d', fold));
    for t = 1:numT
        waitbar(t/numT, h);
        
        trainData = cell2mat(allSubTrainData(:, t, fold));
        testData = cell2mat(allSubTestData(:, t, fold));
        
        [trainData, meanWin, varWin] = zscore(trainData);
        testData = bsxfun(@rdivide, bsxfun(@minus, testData, meanWin), ...
            varWin);
        
%         weightVecs{t, fold} = logReg(trainData, trainLabels, [1 2], 0);
        GNBmodels{t, fold} = nbayes_train(trainData, trainLabels, 1);
        GNBmodels{t, fold}.sortedFeats = sortGNBfeaturesByMuDifference(GNBmodels{t, fold});
        selectedFeats = GNBmodels{t, fold}.sortedFeats(1:numFeats);
        likelihoods = nbayes_apply(testData, GNBmodels{t, fold}, selectedFeats);
        
        [~, predictedLabels] = max(likelihoods, [], 2);
%         P = [testData ones(numTestSamples, 1)]*weightVecs{t, fold};
%         accPerFold(t, fold) = sum(double(P > 0) == testLabels)/numTestSamples;
        accPerFold(t, fold) = sum(predictedLabels == testLabels)/numTestSamples;
        fprintf('Accuracy at %d = %d\n', t, accPerFold(t, fold));
        
    end
    
    close(h);
    
end

save(sprintf('%sKR_Sub_slide.mat', resultRoot), 'accPerFold', 'GNBmodels');%'weightVecs');
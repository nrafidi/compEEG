addpath logisticRegression/
% addpath GNB/
dataRoot = '/Users/nrafidi/Documents/MATLAB/compEEG-data/preproc-final/';
resultRoot = '/Users/nrafidi/Documents/MATLAB/compEEG-data/results/';

numT = 44;
numDim = 128;
numFolds = 5;

load(sprintf('%sallSubData_KR_5fold_withPassLabels.mat', dataRoot));

accPerFold = nan(numT, numFolds);
weightVecs = cell(numT, numFolds);
subIDVecs = cell(numT, numFolds);
%%
for fold = 4:numFolds
    numSub = size(allSubTrainData, 1);
    
    testLabels = cell2mat(allSubTestPassLabels(:,fold));
    trainLabels = cell2mat(allSubTrainPassLabels(:,fold));

    subIDtrainLabels = [];
    subIDtestLabels = [];
    for s = 1:numSub
        subIDtrainLabels = cat(1, subIDtrainLabels, ...
            s*ones(length(allSubTrainLabels{s, fold}), 1));
        subIDtestLabels = cat(1, subIDtestLabels, ...
            s*ones(length(allSubTestLabels{s, fold}), 1));
    end
    
    numTestSamples = length(testLabels);
    h = waitbar(0, sprintf('Fold %d', fold));
    for t = 1:numT
        waitbar(t/numT, h);
        
        trainDataOrig = cell2mat(allSubTrainData(:, t, fold));
        testDataOrig = cell2mat(allSubTestData(:, t, fold));
        
%         [trainDataOrig, meanWin, varWin] = zscore(trainDataOrig);
%         testDataOrig = bsxfun(@rdivide, bsxfun(@minus, testDataOrig, meanWin), ...
%             varWin);
        
        %Regress out sub ID
        subIDVecs{t, fold} = lscov(subIDtrainLabels, trainDataOrig);
        trainData = trainDataOrig - subIDtrainLabels*subIDVecs{t, fold};
        testData = testDataOrig - subIDtestLabels*subIDVecs{t, fold};
        
        
%         [trainData, meanWin, varWin] = zscore(trainData);
%         testData = bsxfun(@rdivide, bsxfun(@minus, testData, meanWin), ...
%             varWin);
        
        weightVecs{t, fold} = logReg(trainData, trainLabels, [1 2], 0);
        
        P = [testData ones(numTestSamples, 1)]*weightVecs{t, fold};
        accPerFold(t, fold) = sum(double(P > 0) == testLabels)/numTestSamples;
        fprintf('Accuracy at %d = %d\n', t, accPerFold(t, fold));
        
    end
    
    close(h);
    
end

save(sprintf('%sKR_gSbS_regS_slide_noZ.mat', resultRoot), 'accPerFold', 'weightVecs', 'subIDVecs');
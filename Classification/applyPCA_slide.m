addpath logisticRegression/
dataRoot = '/Users/nrafidi/Documents/MATLAB/compEEG-data/preproc-final/';
resultRoot = '/Users/nrafidi/Documents/MATLAB/compEEG-data/results/';

numT = 44;
numDim = 128;


load(sprintf('%sallSubData_KR_5fold_withLabels.mat', dataRoot));

accOverComps = nan(numT, 20, 5);
for numCompToUse = 1:20
    for fold = 1:5
        numSub = size(allSubTrainData, 1);
        
        acc_PCA = nan(numT, 1);
        
        testLabels = cell2mat(allSubTestLabels(:,fold));
        trainLabels = cell2mat(allSubTrainLabels(:,fold));
        
        numTestSamples = length(testLabels);
        h = waitbar(0, sprintf('Fold %d', fold));
        for t = 1:numT
            waitbar(t/numT, h);
            allSubTrainData_PCA = cell(numSub, 1);
            allSubTestData_PCA = cell(numSub, 1);
            for s = 1:numSub
                numTrain = size(allSubTrainData{s, t, fold}, 1);
                [coeff, scores] = pca([allSubTrainData{s, t, fold}; allSubTestData{s, t, fold}]);
                allSubTrainData_PCA{s} = scores(1:numTrain,1:numCompToUse);
                allSubTestData_PCA{s} = scores((numTrain+1):end,1:numCompToUse);
                %             allSubTestData_PCA{s} = allSubTestData{s, t, fold}*coeff(:,1:numCompToUse);
            end
            
            
            trainData = cell2mat(allSubTrainData_PCA);
            
            testData = cell2mat(allSubTestData_PCA);
            
            %         [trainData, meanWin, varWin] = zscore(trainData);
            %         testData = bsxfun(@rdivide, bsxfun(@minus, testData, meanWin), ...
            %             varWin);
            
            [w_PCA, lambda] = logReg(trainData, trainLabels, 0, 0);
            %         disp(lambda);
            P_PCA = [testData ones(numTestSamples, 1)]*w_PCA;
            acc_PCA(t) = sum(double(P_PCA > 0) == testLabels)/numTestSamples;
            fprintf('PCA accuracy at %d = %d\n', t, acc_PCA(t));
            %
        end
        
        close(h);
        
        save(sprintf('%sKRClassification_PCAfold%d_slide_numComp%d.mat', ...
            resultRoot, fold, numCompToUse), 'acc_PCA');
        
        accOverComps(:, numCompToUse, fold) = acc_PCA;
    end
    
end
addpath logisticRegression/
% addpath GNB/
dataRoot = '/Users/nrafidi/Documents/MATLAB/compEEG-data/preproc-final/';
resultRoot = '/Users/nrafidi/Documents/MATLAB/compEEG-data/results/';


load(sprintf('%sallSubData_KR.mat', dataRoot));

[numSub, numT] = size(allSubData);

accPerFold = nan(numT, numSub, 3);
regressionFits = cell(numT, numSub, 2);

weightVecs = cell(numT, numSub);
weightVecsGvB = cell(numT, numSub);
weightVecsScoreReg = cell(numT, numSub);
weightVecsGvBAvg = cell(numT, numSub);
weightVecsScoreRegAvg = cell(numT, numSub);


for s = 1:numSub
    trainInd = (1:numSub) ~= s;
    
    trainDataCell = allSubData(trainInd, :);
    testDataCell = allSubData(s, :);
    
    trainDataAvgCell = cellfun(@mean, trainDataCell, 'UniformOutput', false);
    testDataAvgCell = cellfun(@mean, testDataCell, 'UniformOutput', false);
    
    trainLabels = cell2mat(allSubLabels(trainInd,:));
    testLabels = cell2mat(allSubLabels(s,:));
    trainLabelsGvB = cell2mat(allSubPassLabels(trainInd,:));
    testLabelsGvB = cell2mat(allSubPassLabels(s,:));
    trainLabelsScore = cell2mat(allSubScores(trainInd, :));
    testLabelsScore = cell2mat(allSubScores(s, :));
    
    trainLabelsGvBAvg = cell2mat(cellfun(@mean, allSubPassLabels(trainInd, :), 'UniformOutput', false));
    testLabelsGvBAvg = cell2mat(cellfun(@mean, allSubPassLabels(s, :), 'UniformOutput', false));
    trainLabelsScoreAvg = cell2mat(cellfun(@mean, allSubScores(trainInd, :), 'UniformOutput', false));
    testLabelsScoreAvg = cell2mat(cellfun(@mean, allSubScores(s, :), 'UniformOutput', false));
    
    numTestSamples = length(testLabels);
    
    h = waitbar(0, sprintf('Subject %d', s));
    for t = 1:numT
        tic
        waitbar(t/numT, h);
        
        trainData = cell2mat(trainDataCell(:,t));
        testData = cell2mat(testDataCell(:,t));
        
        [trainDataZ, meanWin, varWin] = zscore(trainData);
        testDataZ = bsxfun(@rdivide, bsxfun(@minus, testData, meanWin), ...
            varWin);
        
        trainDataAvg = cell2mat(trainDataAvgCell(:,t));
        testDataAvg = cell2mat(testDataAvgCell(:,t));
        
        [trainDataAvgZ, meanWin, varWin] = zscore(trainDataAvg);
        testDataAvgZ = bsxfun(@rdivide, bsxfun(@minus, testDataAvg, meanWin), ...
            varWin);
        
        %Single Trial classification
        fprintf('1\n');
        weightVecs{t, s} = logReg(trainDataZ, trainLabels, [1 2], 0);
        fprintf('2\n');
        weightVecsGvB{t, s} = logReg(trainDataZ, trainLabelsGvB, [1 2], 0);
        fprintf('3\n');
        weightVecsScoreReg{t, s} = lscov(trainDataZ, trainLabelsScore);
        
        
        P = [testDataZ ones(numTestSamples, 1)]*weightVecs{t, s};
        accPerFold(t, s, 1) = sum(double(P > 0) == testLabels)/numTestSamples;
        
        fprintf('%0.2f\n', accPerFold(t, s, 1));
        
        P = [testDataZ ones(numTestSamples, 1)]*weightVecsGvB{t, s};
        accPerFold(t, s, 2) = sum(double(P > 0) == testLabelsGvB)/numTestSamples;
    
        fprintf('%0.2f\n', accPerFold(t, s, 2));

        regressionFits{t, s, 1} = testDataZ*weightVecsScoreReg{t, s};
        
        %Avg classification
        fprintf('4\n');
        weightVecsGvBAvg{t, s} = logReg(trainDataAvgZ, trainLabelsGvBAvg, [1 2], 0);
        fprintf('5\n');
        weightVecsScoreRegAvg{t, s} = lscov(trainDataAvgZ, trainLabelsScoreAvg);
        
        P = [testDataAvgZ ones(1, 1)]*weightVecsGvBAvg{t, s};
        accPerFold(t, s, 3) = double(double(P > 0) == testLabelsGvBAvg);
        regressionFits{t, s, 2} = testDataAvgZ*weightVecsScoreRegAvg{t, s};
        
        toc
    end
    
    close(h);
    
end

save(sprintf('%sKR_CrossSub_slide.mat', resultRoot), 'accPerFold', 'regressionFits', 'weightVecs', ...
    'weightVecsGvB', 'weightVecsScoreReg', 'weightVecsGvBAvg', 'weightVecsScoreRegAvg');
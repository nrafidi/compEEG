addpath logisticRegression/

subjects = {'AA', 'BB', 'DD', 'F', 'EE', 'GG', 'HH', 'JJ', ...
    'K', 'M', 'O', 'R', 'S', 'T', 'U', 'V', 'X', 'Y', 'Z'};

numSub = length(subjects);
rootDir = '/Users/nrafidi/Documents/MATLAB/compEEG-data/';
dataFile = '%spreproc-final/%s/KR_dataLabels_winSize50-100.mat';

methodToUse = 'mean_concat';%'concat_all', 'mean_concat', 'mean_diff';
runSingleSub = false;

winSizeToUse = 1;%1 or 2
numFolds = 5;
numChans = 64;
numWins = 44;
wordLabels = 1:60;
allSubFoldAssignments = cell(numSub,1);
allSubTrainData = cell(numSub, numWins, numFolds);
allSubTrainLabels = cell(numSub, numFolds);
allSubTrainWordLabels = cell(numSub, numFolds);
allSubTestData = cell(numSub, numWins, numFolds);
allSubTestLabels = cell(numSub, numFolds);
allSubTestWordLabels = cell(numSub, numFolds);
totalNumItems = 0;
withinSubAccs = cell(numSub, 1);
for s = 1:numSub
    sub = subjects{s};
    load(sprintf(dataFile, rootDir, sub));
    subData = krData{winSizeToUse};
    subLabels = krLabels{winSizeToUse};
    
    usableItems = ~cellfun(@isempty, subData);
    numItems = sum(usableItems);
    
    if numItems < 30
        keyboard;
    end
    
    subData = subData(usableItems);
    subWordLabels = wordLabels(usableItems);
    fprintf('Percent Correct = %d\n', (sum(subLabels)/numItems)*100);
    %     subLabels = subLabels(usableItems);
    subFolds = crossvalind('Kfold', numItems, numFolds);
    allSubFoldAssignments{s} = subFolds;
%     h = waitbar(0, sub);
    subAcc = nan(numItems, numWins);
    for f = 1:numFolds
%         waitbar(f/numFolds, h, sub);
        trainInd = subFolds ~= f;
        trainDataCell = subData(trainInd);
        trainLabels = subLabels(trainInd);
        trainWordLabels = subWordLabels(trainInd);
        testDataCell = subData(~trainInd);
        testLabels = subLabels(~trainInd);
        testWordLabels = subWordLabels(~trainInd);
        for w = 1:numWins
            
            startInd = 1 + (w-1)*numChans;
            endInd = startInd + numChans - 1;
            switch methodToUse
                case 'concat_all'
                    trainData = cellfun(@(x) [x(1,startInd:endInd), ...
                        x(2,startInd:endInd), x(3,startInd:endInd), ...
                        x(4,startInd:endInd)], trainDataCell, ...
                        'UniformOutput', false);
                    trainData = cell2mat(trainData');
                    testData = cellfun(@(x) [x(1,startInd:endInd), ...
                        x(2,startInd:endInd), x(3,startInd:endInd), ...
                        x(4,startInd:endInd)], testDataCell, ...
                        'UniformOutput', false);
                    testData = cell2mat(testData');
                case 'mean_concat'
                    trainData = cellfun(@(x) [mean([x(1,startInd:endInd); ...
                        x(2,startInd:endInd)]), mean([x(3,startInd:endInd); ...
                        x(4,startInd:endInd)])], trainDataCell, ...
                        'UniformOutput', false);
                    trainData = cell2mat(trainData');
                    testData = cellfun(@(x) [mean([x(1,startInd:endInd); ...
                        x(2,startInd:endInd)]), mean([x(3,startInd:endInd); ...
                        x(4,startInd:endInd)])], testDataCell, ...
                        'UniformOutput', false);
                    testData = cell2mat(testData');
                case 'mean_diff'
                    trainData = cellfun(@(x) mean([x(1,startInd:endInd); ...
                        x(2,startInd:endInd)]) - mean([x(3,startInd:endInd); ...
                        x(4,startInd:endInd)]), trainDataCell, ...
                        'UniformOutput', false);
                    trainData = cell2mat(trainData');
                    testData = cellfun(@(x) mean([x(1,startInd:endInd); ...
                        x(2,startInd:endInd)]) - mean([x(3,startInd:endInd); ...
                        x(4,startInd:endInd)]), testDataCell, ...
                        'UniformOutput', false);
                    testData = cell2mat(testData');
                case 'round_4'
                    trainData = cellfun(@(x) x(4,startInd:endInd), trainDataCell, ...
                        'UniformOutput', false);
                    trainData = cell2mat(trainData');
                    testData = cellfun(@(x) x(4,startInd:endInd), testDataCell, ...
                        'UniformOutput', false);
                    testData = cell2mat(testData');
            end
            allSubTrainData{s, w, f} = trainData;
            allSubTestData{s, w, f} = testData;
            
            if runSingleSub
                [trainData, mu, sigma] = zscore(trainData);
                testData = bsxfun(@minus, testData, mu);
                testData = bsxfun(@rdivide, testData, sigma);
                
                weights = logReg(trainData, trainLabels, [1,2], 0);
                probs = [testData ones(size(testData,1), 1)]*weights;
                probs = exp(probs)./(1 + exp(probs));
                Yhat = double(probs > 0.5);
                
                subAcc(~trainInd, w) = double(Yhat == testLabels);
            end
        end
        allSubTrainLabels{s, f} = trainLabels;
        allSubTestLabels{s, f} = testLabels;
        allSubTrainWordLabels{s, f} = trainWordLabels';
        allSubTestWordLabels{s, f} = testWordLabels';
    end
    totalNumItems = totalNumItems + numItems;
%     close(h);
    if runSingleSub
        save(sprintf('%sresults/%s/KRClassification_EEG_%s_w%d.mat', ...
            rootDir, sub, methodToUse, winSizeToUse), 'subAcc');
    end
end

save(sprintf('%spreproc-final/allSubData_KR_5fold_withLabels.mat', rootDir), ...
    'allSubTrainData', 'allSubTestData', 'allSubTrainLabels', 'allSubTestLabels', ...
    'allSubFoldAssignments', 'allSubTrainWordLabels', 'allSubTestWordLabels');
% keyboard
%%
% Cross-subject Decoding
crossSubAcc = cell(numWins, 1);
for f = 1:numFolds
    trainLabels = cell2mat(allSubTrainLabels(:,f));
    testLabels = cell2mat(allSubTestLabels(:,f));
    for w = 12%1:numWins
        tic
        trainData = cell2mat(allSubTrainData(:, w, f));
        testData = cell2mat(allSubTestData(:, w, f));
        
        [trainData, mu, sigma] = zscore(trainData);
        testData = bsxfun(@minus, testData, mu);
        testData = bsxfun(@rdivide, testData, sigma);
        
        weights = logReg(trainData, trainLabels, [1,2], 0);
        probs = [testData ones(size(testData,1), 1)]*weights;
        probs = exp(probs)./(1 + exp(probs));
        Yhat = double(probs > 0.5);
        
        mu_z = mean(trainData, 1);
        
        crossSubAcc{w} = cat(1, crossSubAcc{w}, double(Yhat == testLabels));
        save(sprintf('%sresults/weightsMus_krClass_%s_f%d_w%d.mat', ...
            rootDir, methodToUse, f, w), 'weights', 'mu_z');
        toc
    end
end
save(sprintf('%sresults/KRClassification_EEG_CrossSub_%s_w%d_subset.mat', rootDir, ...
    methodToUse, winSizeToUse), 'crossSubAcc');
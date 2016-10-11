addpath logisticRegression/

subjects = {'AA', 'BB', 'DD', 'F', 'EE', 'GG', 'HH', 'JJ', ...
    'K', 'M', 'O', 'R', 'S', 'T', 'U', 'V', 'X', 'Y', 'Z'};

numSub = length(subjects);
rootDir = '/Users/nrafidi/Documents/MATLAB/compEEG-data/';
dataFile = '%spreproc-final/%s/KR_dataLabels_winSize50-100.mat';

methodToUse = 'mean_concat';%'concat_all', 'mean_concat', 'mean_diff';
runSingleSub = false;

winSizeToUse = 1;%1 or 2
numFolds = [];
numChans = 64;
numWins = 44;
wordLabels = 1:60;

if ~isempty(numFolds)
    allSubFoldAssignments = cell(numSub,1);
    allSubTrainData = cell(numSub, numWins, numFolds);
    allSubTrainLabels = cell(numSub, numFolds);
    allSubTrainWordLabels = cell(numSub, numFolds);
    allSubTestData = cell(numSub, numWins, numFolds);
    allSubTestLabels = cell(numSub, numFolds);
    allSubTestWordLabels = cell(numSub, numFolds);
else
    allSubData = cell(numSub, numWins);
    allSubLabels = cell(numSub, 1);
    allSubWordLabels = cell(numSub, 1);
    allSubScores = cell(numSub, 1);
    allSubPassLabels = cell(numSub, 1);
end

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
    allSubWordLabels{s, 1} = wordLabels(usableItems);
    subScore = sum(subLabels)/numItems;
    allSubScores{s, 1} = subScore*ones(size(subLabels, 1), 1);
    fprintf('Percent Correct = %d\n', subScore*100);
    
    if ~isempty(numFolds)
        subFolds = crossvalind('Kfold', numItems, numFolds);
        allSubFoldAssignments{s} = subFolds;
        for f = 1:numFolds %#ok<*BDSCI>
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
                
                
            end
            allSubTrainLabels{s, f} = trainLabels;
            allSubTestLabels{s, f} = testLabels;
            allSubTrainWordLabels{s, f} = trainWordLabels';
            allSubTestWordLabels{s, f} = testWordLabels';
        end
    else
        for w = 1:numWins
            
            startInd = 1 + (w-1)*numChans;
            endInd = startInd + numChans - 1;
            switch methodToUse
                case 'concat_all'
                    subDataWin = cellfun(@(x) [x(1,startInd:endInd), ...
                        x(2,startInd:endInd), x(3,startInd:endInd), ...
                        x(4,startInd:endInd)], subData, ...
                        'UniformOutput', false);
                    subDataWin = cell2mat(subDataWin');
                case 'mean_concat'
                    subDataWin = cellfun(@(x) [mean([x(1,startInd:endInd); ...
                        x(2,startInd:endInd)]), mean([x(3,startInd:endInd); ...
                        x(4,startInd:endInd)])], subData, ...
                        'UniformOutput', false);
                    subDataWin = cell2mat(subDataWin');
                case 'mean_diff'
                    subDataWin = cellfun(@(x) mean([x(1,startInd:endInd); ...
                        x(2,startInd:endInd)]) - mean([x(3,startInd:endInd); ...
                        x(4,startInd:endInd)]), subData, ...
                        'UniformOutput', false);
                    subDataWin = cell2mat(subDataWin');
                case 'round_4'
                    subDataWin = cellfun(@(x) x(4,startInd:endInd), subData, ...
                        'UniformOutput', false);
                    subDataWin = cell2mat(subDataWin');
            end
            allSubData{s, w} = subDataWin;
            
            
        end
        allSubLabels{s, 1} = subLabels;
        if subScore >= 0.5
            allSubPassLabels{s, 1} = ones(length(subLabels), 1);
        else
            allSubPassLabels{s, 1} = zeros(length(subLabels), 1);
        end
    end
    
end

if ~isempty(numFolds)
save(sprintf('%spreproc-final/allSubData_KR_5fold_withLabels.mat', rootDir), ...
    'allSubTrainData', 'allSubTestData', 'allSubTrainLabels', 'allSubTestLabels', ...
    'allSubFoldAssignments', 'allSubTrainWordLabels', 'allSubTestWordLabels');
else
    save(sprintf('%spreproc-final/allSubData_KR.mat', rootDir), ...
    'allSubData', 'allSubLabels', 'allSubScores', ...
    'allSubPassLabels', 'allSubWordLabels');
end
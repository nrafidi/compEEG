addpath ./logisticRegression/

krWinString = '';

repSubjects = {'JJJ', 'III', 'KKK', 'BBB', 'GGG', 'HHH', 'AAA', 'CCC', 'DDD', 'EEE', 'FFF', 'MM', ...
    'OO', 'PP', 'QQ', 'RR', 'SS', 'TT', 'WW',...
    'YY'};

origSubjects = {'AA', 'BB', 'DD', 'F', 'EE', 'GG', 'HH', 'JJ', ...
    'K', 'M', 'O', 'R', 'S', 'T', 'U', 'V', 'X', 'Y', 'Z', 'N'};

numRS = length(repSubjects);
numOS = length(origSubjects);

dataRootR = '/Users/nrafidi/Documents/MATLAB/compEEG-data-rep/';
dataRootO = '/Users/nrafidi/Documents/MATLAB/compEEG-data/';


fnameString = '%sresults/%s/KRanalysis_TGM_Vis%s.mat';

%%
ctWinToUse = 14;
krWinToUse = 37;

%%
IndividualSubjectBrainDataR = cell(numOS + numRS, 1);
IndividualSubjectBrainDataF = cell(numOS + numRS, 1);
IndividualSubjectBehavDataR = cell(numOS + numRS, 1);
IndividualSubjectBehavDataF = cell(numOS + numRS, 1);
IndividualSubjectFirstCorrR = cell(numOS + numRS, 1);
IndividualSubjectFirstCorrF = cell(numOS + numRS, 1);
for i = 1:numRS
    load(sprintf(fnameString, dataRootR, repSubjects{i}, krWinString));
    
    goodItems = sum(responseTraj, 2) > 1 | sum(responseTraj(:,1:3),2) > 0;
    
    for j = 1:size(responseTraj, 1)
        if goodItems(j)
            if krLabels(j) == 1
                IndividualSubjectFirstCorrR{i} = cat(1, IndividualSubjectFirstCorrR{i}, find(responseTraj(j,:), 1, 'first'));
            else
                IndividualSubjectFirstCorrF{i} = cat(1, IndividualSubjectFirstCorrF{i}, find(responseTraj(j,:), 1, 'first'));
            end
        end
    end
    
    krTrajToAdd_R = krTraj(krLabels == 1 & goodItems,:,:,:);
    krTrajToAdd_F = krTraj(krLabels == 0 & goodItems,:,:,:);
    krRespTrajToAdd_R = responseTraj(krLabels == 1 & goodItems,:);
    krRespTrajToAdd_F = responseTraj(krLabels == 0 & goodItems,:);
    
    winData_R = squeeze(krTrajToAdd_R(:, :, krWinToUse, ctWinToUse));
    winData_F = squeeze(krTrajToAdd_F(:, :, krWinToUse,ctWinToUse));
    
    for j = 1:size(winData_R, 1)
        IndividualSubjectBrainDataR{i} = cat(1, IndividualSubjectBrainDataR{i}, ...
            winData_R(j, IndividualSubjectFirstCorrR{i}(j)) - winData_R(j, 4));
    end
    for j = 1:size(winData_F, 1)
        IndividualSubjectBrainDataF{i} = cat(1, IndividualSubjectBrainDataF{i}, ...
            winData_F(j, IndividualSubjectFirstCorrF{i}(j)) - winData_F(j, 4));
    end
    IndividualSubjectBehavDataR{i} = krRespTrajToAdd_R;
    IndividualSubjectBehavDataF{i} = krRespTrajToAdd_F;
end

for i = 1:numOS
    load(sprintf(fnameString, dataRootO, origSubjects{i}, krWinString));
    
    goodItems = sum(responseTraj, 2) > 1 | sum(responseTraj(:,1:3),2) > 0;
    
    for j = 1:size(responseTraj, 1)
        if goodItems(j)
            if krLabels(j) == 1
                IndividualSubjectFirstCorrR{i+numRS} = cat(1, IndividualSubjectFirstCorrR{i+numRS}, find(responseTraj(j,:), 1, 'first'));
            else
                IndividualSubjectFirstCorrF{i+numRS} = cat(1, IndividualSubjectFirstCorrF{i+numRS}, find(responseTraj(j,:), 1, 'first'));
            end
        end
    end
    
    krTrajToAdd_R = krTraj(krLabels == 1 & goodItems,:,:,:);
    krTrajToAdd_F = krTraj(krLabels == 0 & goodItems,:,:,:);
    krRespTrajToAdd_R = responseTraj(krLabels == 1 & goodItems,:);
    krRespTrajToAdd_F = responseTraj(krLabels == 0 & goodItems,:);
    
    winData_R = squeeze(krTrajToAdd_R(:, :, krWinToUse, ctWinToUse));
    winData_F = squeeze(krTrajToAdd_F(:, :, krWinToUse,ctWinToUse));
    
    for j = 1:size(winData_R, 1)
        IndividualSubjectBrainDataR{i+numRS} = cat(1, IndividualSubjectBrainDataR{i+numRS}, ...
            winData_R(j, IndividualSubjectFirstCorrR{i+numRS}(j)) - winData_R(j, 4));
    end
    for j = 1:size(winData_F, 1)
        IndividualSubjectBrainDataF{i+numRS} = cat(1, IndividualSubjectBrainDataF{i+numRS}, ...
            winData_F(j, IndividualSubjectFirstCorrF{i+numRS}(j)) - winData_F(j, 4));
    end
    IndividualSubjectBehavDataR{i+numRS} = krRespTrajToAdd_R;
    IndividualSubjectBehavDataF{i+numRS} = krRespTrajToAdd_F;
end
%% CV to compute bayes factor
rng(12191989)

numSub = numOS + numRS;
numFolds = 5;
numPerms = 1000;

cvInd = crossvalind('Kfold', numSub, numFolds);

foldAccs = cell(numFolds, 1);
foldPermAccs = cell(numFolds, numPerms);
for sub = 1:numFolds
    fprintf('Fold %d out of %d\n', sub, numFolds);
    trainBehavDataR = cell2mat(IndividualSubjectBehavDataR(sub ~= cvInd));
    trainBehavDataF = cell2mat(IndividualSubjectBehavDataF(sub ~= cvInd));
    
    numR = size(trainBehavDataR, 1);
    numF = size(trainBehavDataF, 1);
    
    minNum = min([numR, numF]);
    
    fprintf('Num Samples = %d\n', minNum)
    
    indR = randsample(numR, minNum);
    indF = randsample(numF, minNum);
      
    trainBehavDataR = trainBehavDataR(indR, :);
    trainBehavDataF = trainBehavDataF(indF, :);
    
    trainLabels = [ones(minNum, 1); ...
        zeros(minNum, 1)];
    
    trainDataBehav = [trainBehavDataR; trainBehavDataF];
    
    testDataBehav = [cell2mat(IndividualSubjectBehavDataR(sub == cvInd)); cell2mat(IndividualSubjectBehavDataF(sub == cvInd))];
    
    testLabels = [ones(size(cell2mat(IndividualSubjectBehavDataR(sub == cvInd)), 1), 1); ...
        zeros(size(cell2mat(IndividualSubjectBehavDataF(sub == cvInd)), 1), 1)];
    
    [B_behav, lambda] = logReg(trainDataBehav, trainLabels, [1 10], false);
    fprintf('Lambda behav = %d\n', lambda)
    
    P = [testDataBehav ones(size(testLabels, 1), 1)]*B_behav;
    P = exp(P)./(1 + exp(P));
    Yhat = double(P > 0.5);
    foldAccs{sub, 1} = Yhat == testLabels;
    
    fprintf('Behav accuracy = %0.2f\n', mean(foldAccs{sub, 1}));
    
    for p = 1:numPerms
        trainLabels = trainLabels(randperm(length(trainLabels)));
        testLabels = trainLabels(randperm(length(testLabels)));
        
        [B_behav, lambda] = logReg(trainDataBehav, trainLabels, [1 10], false);
%         fprintf('Lambda behav = %d\n', lambda)

        P = [testDataBehav ones(size(testLabels, 1), 1)]*B_behav;
        P = exp(P)./(1 + exp(P));
        Yhat = double(P > 0.5);
        foldPermAccs{sub, p} = Yhat == testLabels;
        
    end
end

save(sprintf('/Users/nrafidi/Documents/MATLAB/compEEG-data-rep/results/BehavPred_subsamp_%dfold.mat', numFolds), ...
            'foldAccs', 'foldPermAccs');
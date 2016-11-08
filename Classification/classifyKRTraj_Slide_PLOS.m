addpath ./GNB/
addpath ./logisticRegression/

subjects = {'AA', 'BB', 'DD', 'EE', 'F', 'GG', 'HH', 'JJ', ...
    'K', 'M', 'O', 'R', 'S', 'T', 'U', 'V', 'X', 'Y', 'Z', 'N'};
numSub = length(subjects);
compWinToUse = 24;%17, 24:28, 38

modelsToSave = [10, 14];

timePoints = -100:20:960;
numTimePoints = length(timePoints);
numRounds = 4;
algToUse = 'LR';

trajSubCorr = [];
trajSubInc = [];
for s = 1:numSub
    
    load(['../../compEEG-data/results/' subjects{s} ...
        '/KRanalysis_SlidingFeat_CWin' num2str(compWinToUse) '.mat']);
    labels = logical(krLabels);
    
    trajSubCorr = cat(1, trajSubCorr, krTraj(labels,:,:));
    trajSubInc = cat(1, trajSubInc, krTraj(~labels,:,:));
end

if ~isempty(modelsToSave)
    numModels = length(modelsToSave);
    models = cell(numModels, 1);
    labelsToUse = [ones(size(trajSubCorr, 1), 1); zeros(size(trajSubInc, 1), 1)];
    for m = 1:numModels
        models{m}.alg = algToUse;
        dataToUse = [squeeze(trajSubCorr(:,:,m));squeeze(trajSubInc(:,:,m))];
        if strcmp(algToUse, 'GNB')
        elseif strcmp(algToUse, 'LR')
            [B, lambda] = logReg(dataToUse, labelsToUse, [1 2], false);
            
            models{m}.time = timePoints(modelsToSave(m));
            models{m}.weights = B;
            models{m}.reg = lambda;
        else
        end
    end
    save(sprintf('/Users/nrafidi/Documents/MATLAB/compEEG-data/results/classifyKRTraj_models_%s_cWin%d.mat', algToUse, compWinToUse), 'models');
end

%%
numCVReps = 5;
numCVFolds = 5;
maxWinSize = 1;
rng('shuffle');

acc_RvF = nan(numCVReps, numCVFolds, maxWinSize, numTimePoints);
for r = 1:numCVReps
    disp(r);
    folds_Corr = crossvalind('Kfold', size(trajSubCorr, 1), numCVFolds);
    folds_Inc = crossvalind('Kfold', size(trajSubInc, 1), numCVFolds);
    for f = 1:numCVFolds
        trainData = [trajSubCorr(folds_Corr ~= f, :, :); ...
            trajSubInc(folds_Inc ~= f, :, :)];
        trainLabels = [ones(sum(folds_Corr ~= f), 1); ...
            zeros(sum(folds_Inc ~= f), 1)];
        testData = [trajSubCorr(folds_Corr == f, :, :); ...
            trajSubInc(folds_Inc == f, :, :)];
        testLabels = [ones(sum(folds_Corr == f), 1); ...
            zeros(sum(folds_Inc == f), 1)];
        for w = maxWinSize
            h = waitbar(0, sprintf('Fold %d, Rep %d', f, r));
            for t = 1:(numTimePoints - w + 1)
                waitbar(t/(numTimePoints-w+1), h);
                trainDataWin = squeeze(mean(trainData(:, :, t:(t+w-1)), 3));
                testDataWin = squeeze(mean(testData(:, :, t:(t+w-1)), 3));
                
                if strcmp(algToUse, 'GNB')
                    GNBmodel = nbayes_train(trainDataWin, trainLabels + 1, 1);
                    P = nbayes_apply(testDataWin, GNBmodel);
                    [~, predictedLabelIndices] = max(P, [], 2);
                    testLabelHat = GNBmodel.labelVocab(predictedLabelIndices) - 1;
                elseif strcmp(algToUse, 'LR')
                    [B, lambda] = logReg(trainDataWin, trainLabels, [1 2], false);
                    disp(lambda)
                    disp(mean(B))
                    P = [testDataWin ones(length(testDataWin), 1)]*B;
                    P = exp(P)./(1 + exp(P));
                    testLabelHat = double(P > 0.5);
                else
                    [SVMmodel, BoxConstraint] = fitAndTuneSVM(trainDataWin, trainLabels, 'rbf', 2);
                    disp(BoxConstraint)
                    testLabelHat = predict(SVMmodel, testDataWin);
                end
                
                acc_RvF(r, f, w, t) = mean(testLabelHat == testLabels);
            end
            close(h);
        end
    end
end
%%
save(sprintf('/Users/nrafidi/Documents/MATLAB/compEEG-data/results/classifyKRTraj_%s_cWin%d.mat', algToUse, compWinToUse), 'acc_RvF');

figure
acc_RvF = squeeze(mean(mean(acc_RvF, 1), 2));
plot(timePoints, acc_RvF);
title(algToUse);
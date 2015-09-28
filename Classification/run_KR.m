function [ ROC_X, ROC_Y, ROC_T, AUCs ] = run_KR( krTraj, krLabels, ...
    krResp, foldsToUse)
% runTupled_KR: Attempts to predict krLabels from the given krTraj for all
%   possible subsets of the trajectory covariates
%
% Inputs:
%   krTraj: n x p matrix of real-valued KR classifier trajectories (and
%       potentially extra information)
%   krLabels: n x 1 binary vector indicating which
%
% Outputs:
%   accuracies: n x t matrix of 2-fold classification accuracies for the
%       prediction of krLabels from krTraj for each of the t subset tuples of
%       the p covariates.

addpath ./logisticRegression/


AUCs = nan(2, 1);
ROC_X = cell(2, 1);
ROC_Y = cell(2, 1);
ROC_T = cell(2, 1);



% Determine which samples are such that all members of the tuple
% have correct responses
samplesToUse = all(krResp, 2);
% samplesToUse = true(size(krResp, 1), 1);
%         fprintf('Proportion of correct answers in training set = %d\n', ...
%             sum(krLabels(foldsToUse & samplesToUse))/sum(foldsToUse & samplesToUse));
for indTrain = 1:2
    if indTrain == 2
        foldsToUse = ~foldsToUse;
    end
    
    trainIndPos = find(foldsToUse & samplesToUse & krLabels);
    trainIndNeg = find(foldsToUse & samplesToUse & ~krLabels);
    
    numTrainSamp = min([length(trainIndPos), length(trainIndNeg)]);
    trainInd = [trainIndPos(1:numTrainSamp); trainIndNeg(1:numTrainSamp)];
    
    trainData = krTraj(trainInd, :);
    testData = krTraj(~foldsToUse & samplesToUse, :);
    
    numTest = size(testData, 1);
    B = logReg(trainData, krLabels(trainInd), 0, false);
    
    P = [testData ones(numTest, 1)]*B;
    P = exp(P)./(1 + exp(P));
    %             Yhat = double(P > 0.5);
    
    [ROC_X{indTrain}, ROC_Y{indTrain}, ...
        ROC_T{indTrain}, AUCs(indTrain)] = ...
        perfcurve(krLabels(~foldsToUse & samplesToUse), P, 1);
    
    %double(Yhat == krLabels(~foldsToUse));
    
end

end


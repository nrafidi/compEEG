function [ ROC_X, ROC_Y, ROC_T, AUCs, class_ests, true_labs ] = ...
    run_KR( krTraj, krLabels, foldsToUse)
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
class_ests = cell(2,1);
true_labs = cell(2,1);

% numSamps = size(krLabels, 1);
% foldsToUse = false(numSamps, 1);
% % foldsToUse(1:2:(numSamps/2-2)) = true;
% foldsToUse((numSamps/2-1):end) = true;


% Determine which samples are such that all members of the tuple
% have correct responses
% samplesToUse = all(krResp, 2);
% samplesToUse = true(size(krResp, 1), 1);
samplesToUse = true(size(foldsToUse, 1),1);
%         fprintf('Proportion of correct answers in training set = %d\n', ...
%             sum(krLabels(foldsToUse & samplesToUse))/sum(foldsToUse & samplesToUse));
for indTrain = 1:2
    if indTrain == 2
        foldsToUse = ~foldsToUse;
    end
    
    trainIndPos = find(foldsToUse & samplesToUse & krLabels);
    trainIndNeg = find(foldsToUse & samplesToUse & ~krLabels);
%     
%     numTrainSamp = min([length(trainIndPos), length(trainIndNeg)]);
%     trainInd = [trainIndPos(1:numTrainSamp); trainIndNeg(1:numTrainSamp)];
%     
%     if numTrainSamp < 8
%         keyboard;
%     end
    trainInd = [trainIndPos; trainIndNeg];

    trainData = krTraj(trainInd, :);
    testData = krTraj(~foldsToUse & samplesToUse, :);
    
    numTest = size(testData, 1);
    B = logReg(trainData, krLabels(trainInd), 0, false);
    
    P = [testData ones(numTest, 1)]*B;
    P = exp(P)./(1 + exp(P));
    class_ests{indTrain} = P;
    true_labs{indTrain} = krLabels(~foldsToUse & samplesToUse);
    %             Yhat = double(P > 0.5);
    
    [ROC_X{indTrain}, ROC_Y{indTrain}, ...
        ROC_T{indTrain}, AUCs(indTrain)] = ...
        perfcurve(true_labs{indTrain}, P, 1);
    
%     if AUCs(indTrain) < 0.4
%         keyboard;
%     end
    %double(Yhat == krLabels(~foldsToUse));
    
end

fprintf('Number training samples of %d gets AUC of %d\n', length(trainInd), AUCs(2));

end


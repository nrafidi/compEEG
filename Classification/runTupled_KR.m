function [ ROC_X, ROC_Y, ROC_T, AUCs ] = runTupled_KR( krTraj, krLabels, foldsToUse)
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

numCov = size(krTraj, 2);%/2
tuples = {};
totalTuples = 0;
for subsetSize = 1:numCov
    tuples = cat(1, tuples, {nchoosek(1:numCov, subsetSize)});
    %{cat(2, nchoosek(1:numCov, subsetSize), nchoosek(1:numCov, subsetSize) + numCov)});
    totalTuples = totalTuples + nchoosek(numCov, subsetSize);
end
numTuples = length(tuples);

AUCs = nan(2, totalTuples);
ROC_X = cell(2, totalTuples);
ROC_Y = cell(2, totalTuples);
ROC_T = cell(2, totalTuples);
indTuple = 1;

fprintf('Proportion of correct answers in training set = %d\n', ...
    sum(krLabels(foldsToUse))/sum(foldsToUse));
for iTuple = 1:numTuples
    tuple = tuples{iTuple};
    for iWithinTuple = 1:size(tuple, 1)
        for indTrain = 1:2
            if indTrain == 2
                foldsToUse = ~foldsToUse;
            end
            trainData = krTraj(foldsToUse, tuple(iWithinTuple, :));
            testData = krTraj(~foldsToUse, tuple(iWithinTuple, :));
            numTest = size(testData, 1);
            B = logReg(trainData, krLabels(foldsToUse), 0, false);
            
            P = [testData ones(numTest, 1)]*B;
            P = exp(P)./(1 + exp(P));
%             Yhat = double(P > 0.5);
            
            [ROC_X{indTrain, indTuple}, ROC_Y{indTrain, indTuple}, ...
                ROC_T{indTrain, indTuple}, AUCs(indTrain,indTuple)] = ...
                perfcurve(krLabels(~foldsToUse), P, 1);
            
            %double(Yhat == krLabels(~foldsToUse));
            
        end
        indTuple = indTuple + 1;
    end
end

end


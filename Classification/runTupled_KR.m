function [ accuracies ] = runTupled_KR( krTraj, krLabels )
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

[numSamp, numCov] = size(krTraj);
tuples = {};
totalTuples = 0;
for subsetSize = 1:numCov
    tuples = cat(1, tuples, {nchoosek(1:numCov, subsetSize)});
    totalTuples = totalTuples + nchoosek(numCov, subsetSize);
end
numTuples = length(tuples);

accuracies = nan(numSamp, totalTuples);
indTuple = 1;
for iTuple = 1:numTuples
    tuple = tuples{iTuple};
    for iWithinTuple = 1:size(tuple, 1)
        for indTrain = 1:2
            switch indTrain
                case 1
                    indTest = 2;
                case 2
                    indTest = 1;
            end
            trainData = krTraj(indTrain:2:end, tuple(iWithinTuple, :));
            testData = krTraj(indTest:2:end, tuple(iWithinTuple, :));
            numTest = size(testData, 1);
            B = logReg(trainData, krLabels(indTrain:2:end), 0, false);
            
            P = [testData ones(numTest, 1)]*B;
            P = exp(P)./(1 + exp(P));
            Yhat = double(P > 0.5);
            
            accuracies(indTest:2:end,indTuple) = double(Yhat == krLabels(indTest:2:end));
            
        end
        indTuple = indTuple + 1;
    end
end

end


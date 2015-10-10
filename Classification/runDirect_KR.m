function [ ROC_X, ROC_Y, ROC_T, AUCs] = ...
    runDirect_KR( krTraj, krLabels)
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


[ROC_X, ROC_Y, ...
        ROC_T, AUCs] = ...
        perfcurve(krLabels, krTraj, 0);


end


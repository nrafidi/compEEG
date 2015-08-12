function [ trueAcc, populationAcc ] = runSubjBootstrap_KR(...
    krTrajList, krLabelList, numDraws)
% runSubjBootstrap runs tupled KR for the subjects in the given krTrajList
%   and then creates a population histogram of results by sampling those
%   subjects with replacement numDraws times.

trueAcc = runTupled_KR(cell2mat(krTrajList), cell2mat(krLabelList));

populationAcc = cell(numDraws, 1);

numSub = length(krTrajList);
for iDraw = 1:numDraws
    subjectsToUse = datasample(1:numSub, numSub);
    populationAcc{iDraw} = runTupled_KR(cell2mat(krTrajList(subjectsToUse)), ...
        cell2mat(krLabelList(subjectsToUse)));
end

end


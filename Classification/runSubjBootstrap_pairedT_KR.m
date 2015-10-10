function [ p_true, p_pop, draws_pop, mean_diff_true, mean_diff_pop ] = runSubjBootstrap_pairedT_KR(...
    krTrajList, krLabelList, numDraws)
% runSubjBootstrap runs tupled KR for the subjects in the given krTrajList
%   and then creates a population histogram of results by sampling those
%   subjects with replacement numDraws times.

% pre-calculate folds so that you don't accidentally double dip
numSub = length(krTrajList);

trajList = cell2mat(krTrajList);
labelList = cell2mat(krLabelList);
numPos = sum(labelList);
numNeg = length(labelList) - numPos;
numToUse = min([numPos, numNeg]);
posInd = find(labelList);
negInd = find(~labelList);
[~, p_true] = ttest(trajList(posInd(1:numToUse), :), trajList(negInd(1:numToUse),:));

mean_diff_true = mean(trajList(posInd, :), 1) - mean(trajList(negInd,:), 1);

p_pop = nan(numDraws, 1);
mean_diff_pop = nan(numDraws, 1);
draws_pop = cell(numDraws, 1);

for iDraw = 1:numDraws
    subjectsToUse = datasample(1:numSub, numSub);
    draws_pop{iDraw} = subjectsToUse;
    trajList = cell2mat(krTrajList(subjectsToUse));
    labelList = cell2mat(krLabelList(subjectsToUse));
    numPos = sum(labelList);
    numNeg = length(labelList) - numPos;
    numToUse = min([numPos, numNeg]);
    posInd = find(labelList);
    negInd = find(~labelList);
    [~, p_pop(iDraw)] = ttest(trajList(posInd(1:numToUse), :), trajList(negInd(1:numToUse),:));
    mean_diff_pop(iDraw) = mean(trajList(posInd, :), 1) - mean(trajList(negInd,:), 1);
end

end


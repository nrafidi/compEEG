function [ ROC_X_true, ROC_Y_true, ROC_T_true, AUCs_true, ...
    ROC_X_pop, ROC_Y_pop, ROC_T_pop, AUCs_pop ] = runSubjBootstrap_KR(...
    krTrajList, krLabelList, krRespList, numDraws)
% runSubjBootstrap runs tupled KR for the subjects in the given krTrajList
%   and then creates a population histogram of results by sampling those
%   subjects with replacement numDraws times.

% pre-calculate folds so that you don't accidentally double dip
numSub = length(krTrajList);
foldsPerSub = cell(numSub, 1);
for iSub = 1:numSub
    numInst = size(krTrajList{iSub}, 1);
    foldsPerSub{iSub} = false(numInst, 1);
    foldsPerSub{iSub}(1:2:end) = true;
end

[ROC_X_true, ROC_Y_true, ROC_T_true, AUCs_true] = runTupled_KR(cell2mat(krTrajList), cell2mat(krLabelList), ...
    cell2mat(krRespList), cell2mat(foldsPerSub));

ROC_X_pop = cell(numDraws, 1);
ROC_Y_pop = cell(numDraws, 1);
ROC_T_pop = cell(numDraws, 1);
AUCs_pop = cell(numDraws, 1);


for iDraw = 1:numDraws
    subjectsToUse = datasample(1:numSub, numSub);
    [ROC_X_pop{iDraw}, ROC_Y_pop{iDraw}, ROC_T_pop{iDraw}, AUCs_pop{iDraw}] = runTupled_KR(cell2mat(krTrajList(subjectsToUse)), ...
        cell2mat(krLabelList(subjectsToUse)), ...
        cell2mat(krRespList(subjectsToUse)), ...
        cell2mat(foldsPerSub(subjectsToUse)));
end

end


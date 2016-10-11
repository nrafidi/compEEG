function [ ROC_X_true, ROC_Y_true, ROC_T_true, AUCs_true, ...
    ROC_X_pop, ROC_Y_pop, ROC_T_pop, AUCs_pop, draws_pop] ...
    = runSubjBootstrap_KR(...
    krTrajList, krLabelList, numDraws, doDirect, doTuple, trueInd)
% runSubjBootstrap runs tupled KR for the subjects in the given krTrajList
%   and then creates a population histogram of results by sampling those
%   subjects with replacement numDraws times.

% pre-calculate folds so that you don't accidentally double dip
numSub = length(krTrajList);
if ~doDirect
    foldsPerSub = cell(numSub, 1);
    for iSub = 1:numSub
        foldsPerSub{iSub} = rand(size(krTrajList{iSub}, 1), 1) > 0.5;
    end
end

if doDirect
    [ROC_X_true, ROC_Y_true, ROC_T_true, AUCs_true] = ...
        runDirect_KR(cell2mat(krTrajList), cell2mat(krLabelList));
else
    if doTuple
        [ROC_X_true, ROC_Y_true, ROC_T_true, AUCs_true] = ...
            runTupled_KR(cell2mat(krTrajList), cell2mat(krLabelList), ...
            cell2mat(foldsPerSub));
    else
        [ROC_X_true, ROC_Y_true, ROC_T_true, AUCs_true] = ...
            run_KR(cell2mat(krTrajList), cell2mat(krLabelList), ...
            cell2mat(foldsPerSub), trueInd);%, ...
        %     cell2mat(krRespList), cell2mat(foldsPerSub));
    end
end

ROC_X_pop = cell(numDraws, 1);
ROC_Y_pop = cell(numDraws, 1);
ROC_T_pop = cell(numDraws, 1);
AUCs_pop = cell(numDraws, 1);
draws_pop = cell(numDraws, 1);

for iDraw = 1:numDraws
    subjectsToUse = datasample(1:numSub, numSub);
    draws_pop{iDraw} = subjectsToUse;
    
    if doDirect
        [ROC_X_pop{iDraw}, ROC_Y_pop{iDraw}, ROC_T_pop{iDraw}, AUCs_pop{iDraw}] = ...
            runDirect_KR(cell2mat(krTrajList(subjectsToUse)), ...
            cell2mat(krLabelList(subjectsToUse)));
    else
        if doTuple
            [ROC_X_pop{iDraw}, ROC_Y_pop{iDraw}, ROC_T_pop{iDraw}, AUCs_pop{iDraw}] = ...
                runTupled_KR(cell2mat(krTrajList(subjectsToUse)), ...
                cell2mat(krLabelList(subjectsToUse)), ...
                cell2mat(foldsPerSub(subjectsToUse)));
        else
            [ROC_X_pop{iDraw}, ROC_Y_pop{iDraw}, ROC_T_pop{iDraw}, AUCs_pop{iDraw}, ...
                ~, ~] = run_KR(cell2mat(krTrajList(subjectsToUse)), ...
                cell2mat(krLabelList(subjectsToUse)), ...
                cell2mat(foldsPerSub(subjectsToUse)), trueInd);%, ...
            %         cell2mat(krRespList(subjectsToUse)), ...
            %         cell2mat(foldsPerSub(subjectsToUse)));
        end
    end
end

end


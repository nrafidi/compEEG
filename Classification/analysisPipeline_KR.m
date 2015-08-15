% Full KR Analysis pipeline

% Starts with Competition and KR EEG data and produces KR population
% accuracy and histogram
%
subjects = {'AA', 'BB', 'DD', 'EE', 'F', 'GG', 'HH', 'JJ', ...
    'K', 'M', 'N', 'O', 'R', 'S', 'T', 'U', 'V', 'X', 'Y', 'Z'};
numSub = length(subjects);

behaveDataRoot = '/Users/nrafidi/Documents/MATLAB/compEEG-data/behavioral/';

krTrajList = cell(numSub, 1);
krLabelList = cell(numSub, 1);
accuracyList = [];

numDraws = 100;

for s = 1:numSub
    sub = subjects{s};
    disp(sub);
    
    fname = ['/Users/nrafidi/Documents/MATLAB/compEEG-data/results/' ...
        sub '/KRanalysis.mat'];
    
    if ~exist(fname, 'file')
        
        try
            %% Predict KR
            itemTraj = runKRPrediction(sub);
            
            %% Combine with final quiz answers
            [itemTrajCorr, itemTrajInc, corrAnswers] = sortKRTraj(sub);
            title(sub);
            
            %% Sort for bootstrap analysis
            
            krTraj = [];
            krLabels = [];
            skippedItems = [];
            indCorr = 1;
            indInc = 1;
            for i = 1:length(corrAnswers)
                if corrAnswers(i)
                    if ~any(isnan(itemTrajCorr(indCorr,:)))
                        krTraj = cat(1, krTraj, itemTrajCorr(indCorr,:));
                        krLabels = cat(1, krLabels, 1);
                    else
                        skippedItems = cat(1, skippedItems, i);
                    end
                    indCorr = indCorr + 1;
                else
                    if ~any(isnan(itemTrajInc(indInc,:)))
                        krTraj = cat(1, krTraj, itemTrajInc(indInc,:));
                        krLabels = cat(1, krLabels, 0);
                    else
                        skippedItems = cat(1, skippedItems, i);
                    end
                    indInc = indInc + 1;
                end
            end
            
            save(fname, 'krTraj', 'krLabels', 'skippedItems');
        catch
            fprintf('Subject %s Failed\n', sub);
        end
    else
        load(fname);
    end
    
    load([behaveDataRoot '/' sub '/' sub '_answerTraj.mat']);
    
    numPotTraj = size(responseTraj, 1);
    newResponseTraj = [];
    responseTraj(isnan(responseTraj)) = 0; %#ok<*SAGROW>
    %~any(isnan(responseTraj(iPTraj,:))) && 
    for iPTraj = 1:numPotTraj
        if ~any(iPTraj == skippedItems);
            newResponseTraj = cat(1, newResponseTraj, responseTraj(iPTraj,:));
        end
    end
    
    
    krTraj = newResponseTraj;%cat(2, krTraj, newResponseTraj);
    
    fprintf('Proportion Correct for subject %s = %d\n', sub, ...
        sum(krLabels)/length(krLabels));
    
    krTrajList{s} = krTraj;
    krLabelList{s} = krLabels;
end

[ ROC_X_true, ROC_Y_true, ROC_T_true, AUCs_true, ...
    ROC_X_pop, ROC_Y_pop, ROC_T_pop, AUCs_pop ] = runSubjBootstrap_KR(...
    krTrajList, krLabelList, numDraws);

save ../../compEEG-data/results/KR_analysis_output_onlyAnsInfo.mat ...
    ROC_X_true ROC_Y_true ROC_T_true AUCs_true ...
    ROC_X_pop ROC_Y_pop ROC_T_pop AUCs_pop
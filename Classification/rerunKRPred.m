% Rerun KR Prediction

subjects = {'AA', 'BB', 'DD', 'EE', 'F', 'GG', 'HH', 'JJ', ...
    'K', 'M', 'N', 'O', 'R', 'S', 'T', 'U', 'V', 'X', 'Y', 'Z'};
numSub = length(subjects);

behaveDataRoot = '/Users/nrafidi/Documents/MATLAB/compEEG-data/behavioral/';

krTrajList = cell(numSub, 1);
krRespList = cell(numSub, 1);
krLabelList = cell(numSub, 1);
accuracyList = [];

numDraws = 100;

for s = 1:numSub
    sub = subjects{s};
    disp(sub);
    
    fname = ['/Users/nrafidi/Documents/MATLAB/compEEG-data/results/' ...
        sub '/KRanalysis_logit.mat'];
    
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
        
    end
    
end
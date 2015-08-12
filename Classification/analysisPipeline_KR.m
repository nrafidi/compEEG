% Full KR Analysis pipeline

% Starts with Competition and KR EEG data and produces KR population
% accuracy and histogram

subjects = {'AA', 'BB', 'DD', 'EE', 'F', 'G', 'GG', 'HH', 'JJ', ...
    'K', 'M', 'N', 'O', 'R', 'S', 'T', 'U', 'V', 'X', 'Y', 'Z'};
numSub = length(subjects);

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
            indCorr = 1;
            indInc = 1;
            for i = 1:length(corrAnswers)
                if corrAnswers(i)
                    if ~any(isnan(itemTrajCorr(indCorr,:)))
                        krTraj = cat(1, krTraj, itemTrajCorr(indCorr,:));
                        krLabels = cat(1, krLabels, 1);
                    end
                    indCorr = indCorr + 1;
                else
                    if ~any(isnan(itemTrajInc(indInc,:)))
                        krTraj = cat(1, krTraj, itemTrajInc(indInc,:));
                        krLabels = cat(1, krLabels, 0);
                    end
                    indInc = indInc + 1;
                end
            end
            
            save(fname, 'krTraj', 'krLabels');
        catch
            fprintf('Subject %s Failed\n', sub);
        end
    else
        load(fname);
    end
    
    
    krTrajList{s} = krTraj;
    krLabelList{s} = krLabels;
end

[ trueAcc, populationAcc ] = runSubjBootstrap_KR(...
    krTrajList, krLabelList, numDraws);
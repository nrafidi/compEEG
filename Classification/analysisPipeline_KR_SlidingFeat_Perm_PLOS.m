% Paths
addpath ../Preprocessing/
behaveDataRoot = '/Users/nrafidi/Documents/MATLAB/compEEG-data/behavioral/';
fPrefix = '/Users/nrafidi/Documents/MATLAB/compEEG-data/preproc-final/';

% Parameters
subjects = {'AA', 'BB', 'DD', 'F', 'EE', 'GG', 'HH', 'JJ', ...
    'K', 'M', 'O', 'R', 'S', 'T', 'U', 'V', 'X', 'Y', 'Z', 'N'};
numSub = length(subjects);
fSuffix = '_Vis_BP2-200_N60_Ref_Epochs_Base_ICA1-2.mat';
KRanswer = importdata('/Users/nrafidi/Documents/MATLAB/compEEG-stim/KR_answer.txt');

compWinToUse = 17;
krWinToUse = [650, 670];
numPerms = 100;
for s = 1:numSub
    sub = subjects{s};
    fname = ['/Users/nrafidi/Documents/MATLAB/compEEG-data/results/' ...
        sub '/KRanalysis_SlidingFeat_Perm_CWin' num2str(compWinToUse) '.mat'];
    
    if ~exist(fname, 'file')
        tic
        loadFname = [fPrefix sub '/CompEEG__KR_' sub fSuffix];
        
        krTraj = [];
        krLabels = [];
        [featData, labels, winTime, ~] = extractFeatures(loadFname);
        featData = double(featData);
        %% Predict KR
        itemTraj = runKRPrediction_SlidingFeat_Perm(sub, featData, labels, compWinToUse, krWinToUse);
        for p = 2:numPerms
            itemTraj = cat(4, itemTraj, runKRPrediction_SlidingFeat_Perm(sub, featData, labels, compWinToUse, krWinToUse));
        end
        
        %% Combine with final quiz answers
        load(['../../compEEG-data/results/answers/' sub '_answers.mat']);
        numQ = length(answerList);
        corrAnswers = false(numQ, 1);
        for a = 1:numQ
            corrAnswers(a) = strcmpi(answerList{a}, KRanswer{a});
        end
        
        itemTrajCorr = itemTraj(corrAnswers, :, :, :);
        itemTrajInc = itemTraj(~corrAnswers, :, :, :);
        
        skippedItems = [];
        indCorr = 1;
        indInc = 1;
        for i = 1:length(corrAnswers)
            if corrAnswers(i)
                if ~any(isnan(itemTrajCorr(indCorr,:,:,:)))
                    krTraj = cat(1, krTraj, itemTrajCorr(indCorr,:,:,:));
                    krLabels = cat(1, krLabels, 1);
                else
                    skippedItems = cat(1, skippedItems, i);
                end
                indCorr = indCorr + 1;
            else
                if ~any(isnan(itemTrajInc(indInc,:,:,:)))
                    krTraj = cat(1, krTraj, itemTrajInc(indInc,:,:,:));
                    krLabels = cat(1, krLabels, 0);
                else
                    skippedItems = cat(1, skippedItems, i);
                end
                indInc = indInc + 1;
            end
            
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
        responseTraj = newResponseTraj;
        save(fname, 'krTraj', 'krLabels', 'skippedItems', 'responseTraj', 'winTime');
        fprintf('Subject %s succeeded.\n', sub);
        toc
    end
end
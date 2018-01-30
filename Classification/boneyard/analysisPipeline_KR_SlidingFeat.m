subjects = {'AA', 'BB', 'DD', 'F', 'EE', 'GG', 'HH', 'JJ', ...
        'K', 'M', 'O', 'R', 'S', 'T', 'U', 'V', 'X', 'Y', 'Z', 'N'};

% subjects = {'M'};

numSub = length(subjects);
addpath ../Preprocessing/
behaveDataRoot = '/Users/nrafidi/Documents/MATLAB/compEEG-data/behavioral/';
fPrefix = '/Users/nrafidi/Documents/MATLAB/compEEG-data/preproc-final/';
% fSuffix = '_Vis_Hilbert-theta_Epochs.mat';
fSuffix = '_Vis_Epochs_ICA1-2.mat';
KRanswer = importdata('/Users/nrafidi/Documents/MATLAB/compEEG-stim/KR_answer.txt');
winSizeOptions = 50;
compWinToUse = 21 + 14;
numSizes = length(winSizeOptions);
for s = 1:numSub
    sub = subjects{s};
    disp(sub);
    loadFname = [fPrefix sub '/CompEEG__KR_' sub fSuffix];
    fname = ['/Users/nrafidi/Documents/MATLAB/compEEG-data/results/' ...
        sub '/KRanalysis_SlidingFeat_lateComp_preStim.mat'];
%         sub '/KRanalysis_SlidingFeat_Hilbert-theta.mat'];
        
    krTraj = cell(numSizes, 1);
    krLabels = cell(numSizes, 1);
    %     try
    options = struct;
    options.overLap = true;
    options.minTime = -300;
    for w = 1:numSizes
        options.erpWinSize = winSizeOptions(w);
        [featData, labels, ~] = extractFeatures(loadFname, options);
        featData = double(featData);
        %% Predict KR
        itemTraj = runKRPrediction_SlidingFeat(sub, featData, labels, compWinToUse);
        
        %% Combine with final quiz answers
        load(['../../compEEG-data/results/answers/' sub '_answers.mat']);
        numQ = length(answerList);
        corrAnswers = false(numQ, 1);
        for a = 1:numQ
            corrAnswers(a) = strcmpi(answerList{a}, KRanswer{a});
        end
        
        itemTrajCorr = itemTraj(corrAnswers, :, :);
        itemTrajInc = itemTraj(~corrAnswers, :, :);
        
        skippedItems = [];
        indCorr = 1;
        indInc = 1;
        for i = 1:length(corrAnswers)
            if corrAnswers(i)
                if ~any(isnan(itemTrajCorr(indCorr,:)))
                    krTraj{w} = cat(1, krTraj{w}, itemTrajCorr(indCorr,:, :));
                    krLabels{w} = cat(1, krLabels{w}, 1);
                else
                    skippedItems = cat(1, skippedItems, i);
                end
                indCorr = indCorr + 1;
            else
                if ~any(isnan(itemTrajInc(indInc,:)))
                    krTraj{w} = cat(1, krTraj{w}, itemTrajInc(indInc,:, :));
                    krLabels{w} = cat(1, krLabels{w}, 0);
                else
                    skippedItems = cat(1, skippedItems, i);
                end
                indInc = indInc + 1;
            end
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
%     keyboard;
    save(fname, 'krTraj', 'krLabels', 'skippedItems', 'responseTraj');
    
end
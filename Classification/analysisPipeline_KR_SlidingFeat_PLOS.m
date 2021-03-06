% Paths
addpath ../Preprocessing/
behaveDataRoot = '/Users/nrafidi/Documents/MATLAB/compEEG-data-rep/behavioral/';
fPrefix = '/Users/nrafidi/Documents/MATLAB/compEEG-data-rep/preproc-final/';

% Parameters
subjects = {'JJJ', 'III', 'KKK', 'BBB', 'GGG', 'HHH', 'AAA', 'CCC', 'DDD', 'EEE', 'FFF', 'MM', ...
    'OO', 'PP', 'QQ', 'RR', 'SS', 'TT', 'WW',...
    'YY'};

%VV and UU have files missing
% {'AA', 'BB', 'DD', 'F', 'EE', 'GG', 'HH', 'JJ', ...
%     'K', 'M', 'O', 'R', 'S', 'T', 'U', 'V', 'X', 'Y', 'Z', 'N'};
numSub = length(subjects);
fSuffixString = '%s_BP2-%d_N60_Ref_Epochs_Base_ICA1-2.mat';
KRanswer = importdata('/Users/nrafidi/Documents/MATLAB/compEEG-stim/KR_answer.txt');

compWinToUse = 18; %17, 24:29, 38
for s = 1:numSub
    sub = subjects{s};
    disp(sub);
    
    if strcmp(sub, 'NN')
        fSuffix = sprintf(fSuffixString, '_Vis', 127);
    elseif strcmp(sub, 'HHH')
        fSuffix = sprintf(fSuffixString, '', 200);
    else
        fSuffix = sprintf(fSuffixString, '_Vis', 200);
    end
    
    loadFname = [fPrefix sub '/CompEEG__KR_' sub fSuffix];
    fname = ['/Users/nrafidi/Documents/MATLAB/compEEG-data-rep/results/' ...
        sub '/KRanalysis_SlidingFeat_CWin' num2str(compWinToUse) '_Vis.mat'];
    
    krTraj = [];
    krLabels = [];
    [featData, labels, winTime, ~] = extractFeatures(loadFname);
    featData = double(featData);
    %% Predict KR
    itemTraj = runKRPrediction_SlidingFeat(sub, featData, labels, compWinToUse);
    
    %% Combine with final quiz answers
    load(['../../compEEG-data-rep/results/answers/' sub '_answers.mat']);
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
                krTraj = cat(1, krTraj, itemTrajCorr(indCorr,:, :));
                krLabels = cat(1, krLabels, 1);
            else
                skippedItems = cat(1, skippedItems, i);
            end
            indCorr = indCorr + 1;
        else
            if ~any(isnan(itemTrajInc(indInc,:)))
                krTraj = cat(1, krTraj, itemTrajInc(indInc,:, :));
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
    %     keyboard;
    save(fname, 'krTraj', 'krLabels', 'skippedItems', 'responseTraj', 'winTime');
    
end
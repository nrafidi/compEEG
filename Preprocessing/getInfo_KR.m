% subjects = {'AA', 'BB', 'DD', 'F', 'EE', 'GG', 'HH', 'JJ', ...
%     'K', 'M', 'O', 'R', 'S', 'T', 'U', 'V', 'X', 'Y', 'Z', 'N'};

subjects = {'N'};

numSub = length(subjects);
behaveDataRoot = '/Users/nrafidi/Documents/MATLAB/compEEG-data/behavioral/';
fPrefix = '/Users/nrafidi/Documents/MATLAB/compEEG-data/preproc-final/';
fSuffix = '_Vis_BP2-200_N60_Ref_Epochs_Base_ICA1-2.mat';
KRanswer = importdata('/Users/nrafidi/Documents/MATLAB/compEEG-stim/KR_answer.txt');
winSizeOptions = [50, 100];
numSizes = length(winSizeOptions);
options = struct;
options.overLap = true;
options.minTime = -100;
for s = 1:numSub
    sub = subjects{s};
    disp(sub);
    loadFname = sprintf('%s%s/CompEEG__KR_%s%s', fPrefix, sub, sub, fSuffix);
    fname = sprintf('%s/%s/KR_dataLabels_winSize%d-%d.mat', ...
        fPrefix, sub, min(winSizeOptions), max(winSizeOptions));
    
    krData = cell(numSizes, 1);
    krLabels = cell(numSizes, 1);
    
    for w = 1:numSizes
        %% Extract Features from KR data
        options.erpWinSize = winSizeOptions(w);
        [featData, labels, ~] = extractFeatures(loadFname, options);
        featData = double(featData);
        featData = featData(labels(:,1) == 1,:);
        labels = labels(labels(:,1)==1, 2);
        
        %% Combine with final quiz answers
        load(sprintf('../../compEEG-data/results/answers/%s_answers.mat', sub));
        numQ = length(answerList);
        corrAnswers = false(numQ, 1);
        for a = 1:numQ
            corrAnswers(a) = strcmpi(answerList{a}, KRanswer{a});
        end
        
        uniqueItems = 1:60;
        skippedItems = [];
        correctItems = uniqueItems(corrAnswers);
        
        for i = uniqueItems
            
            indItem = labels == i;
            if sum(indItem) == 4
                krData{w}{i} = featData(indItem, :);
                krLabels{w} = cat(1, krLabels{w}, corrAnswers(i));
            else
                skippedItems = cat(1, skippedItems, i);
            end
        end
    end
    load(sprintf('%s/%s/%s_answerTraj.mat', behaveDataRoot, sub, sub));
    
    numPotTraj = size(responseTraj, 1);
    newResponseTraj = [];
    responseTraj(isnan(responseTraj)) = 0; %#ok<*SAGROW>
    for iPTraj = 1:numPotTraj
        if ~any(iPTraj == skippedItems);
            newResponseTraj = cat(1, newResponseTraj, responseTraj(iPTraj,:));
        end
    end
    responseTraj = newResponseTraj;
    
    save(fname, 'krData', 'krLabels', 'skippedItems', 'responseTraj');
    fprintf('%s Complete\n', sub);
end
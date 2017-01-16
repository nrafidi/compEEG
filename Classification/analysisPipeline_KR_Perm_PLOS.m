% Full KR Analysis pipeline
function analysisPipeline_KR_Perm_PLOS
addpath ./logisticRegression/
subjects = {'AA', 'BB', 'DD', 'EE', 'F', 'GG', 'HH', 'JJ', ...
    'K', 'M', 'N', 'O', 'R', 'S', 'T', 'U', 'V', 'X', 'Y', 'Z'};
numSub = length(subjects);
eegDataRoot = '/Users/nrafidi/Documents/MATLAB/compEEG-data/preproc-final/';
resRoot = '/Users/nrafidi/Documents/MATLAB/compEEG-data/results/';
KRanswer = importdata('/Users/nrafidi/Documents/MATLAB/compEEG-stim/KR_answer.txt');

numDraws = 1000;
numPerms = 100;
compWinToUse = 17;
krWinToUse = 39;
numChans = 64;
numRounds = 4;
trueNumItems = 60;

dataFname = sprintf('%s/allSub_CompKRData_cWin%d_krWin%d.mat', ...
    eegDataRoot, compWinToUse, krWinToUse);
trueAccFname = sprintf('%s/trueAccBetas_cWin%d_krWin%d.mat', ...
    resRoot, compWinToUse, krWinToUse);
drawAccFname = '%s/subDrawsPerms/draw%d.mat';

if exist(dataFname, 'file')
    load(dataFname);
else
    compDataList = cell(numSub, 1);
    compLabelList = cell(numSub, 1);
    krDataList = cell(numSub, 1);
    krLabelList = cell(numSub, 1);
    for s = 1:numSub
        sub = subjects{s};
        
        % Load and store Competition Data
        load(sprintf('%s/%s/CompEEG_%s_Vis_BP2-200_N60_Ref_Epochs_Base_ICA1-2_Features_Overlap_Time.mat', ...
            eegDataRoot, sub, sub));
        startInd = (compWinToUse-1)*numChans + 1;
        endInd = compWinToUse*numChans;
        compDataList{s} = double(featData(:,startInd:endInd));
        compLabelList{s} = labels(:,1);
        
        % Load KR Data
        load(sprintf('%s/%s/CompEEG__KR_%s_Vis_BP2-200_N60_Ref_Epochs_Base_ICA1-2_Features_Overlap_Time.mat', ...
            eegDataRoot, sub, sub));
        startInd = (krWinToUse-1)*numChans + 1;
        endInd = krWinToUse*numChans;
        featData = double(featData(labels(:,1) == 1, startInd:endInd));
        itemLabels = labels(labels(:,1) == 1, 2);
        uniqueItems = unique(itemLabels);
        uniqueItems = uniqueItems(1:trueNumItems);
        numItems = length(uniqueItems);
        
        % Load KR Labels
        load(sprintf('%s/answers/%s_answers.mat', resRoot, sub));
        corrAnswers = false(numItems, 1);
        for a = 1:numItems
            corrAnswers(a) = strcmpi(answerList{a}, KRanswer{a});
        end
        
        % For which items do we have all 4 rounds of data? Store Data and
        % Labels for those
        itemInds = [];
        itemsToUse = false(numItems, 1);
        for i = 1:numItems
            findItem = find(itemLabels == uniqueItems(i))';
            if length(findItem) == numRounds
                itemInds = cat(1, itemInds, findItem);
                itemsToUse(i) = true;
            end
        end
        numItems = size(itemInds, 1);
        krDataList{s} = nan(numItems, numRounds, numChans);
        for i = 1:numItems
            krDataList{s}(i,:,:) = featData(itemInds(i,:), :);
        end
        krLabelList{s} = double(corrAnswers(itemsToUse));
        %         keyboard;
    end
    save(dataFname, ...
        'compDataList', 'krDataList', 'compLabelList', 'krLabelList');
end

% Compute the true AUC and save the true models
if exist(trueAccFname, 'file')
    load(trueAccFname);
else
    [trueAcc, betas] = computeAcc(compDataList, compLabelList, krDataList, ...
        krLabelList, false);
    save(trueAccFname, 'trueAcc', 'betas');
end

% For each draw, sample subjects with replacement
trueBootAccs = nan(numDraws, 1);
permBootAccs = nan(numDraws, numPerms);
for d = 1:numDraws
    if exist(sprintf(drawAccFname, resRoot, d), 'file')
        load(sprintf(drawAccFname, resRoot, d));
        trueBootAccs(d) = trueBootAcc;
        permBootAccs(d,:) = permBootAcc;
    else
        rng(d);
        fprintf('Draw %d/1000\n', d);
        subjectsToUse = datasample(1:numSub, numSub);
        trueBootAccs(d) = computeAcc(compDataList, compLabelList, krDataList, ...
            krLabelList, false, betas(subjectsToUse));
        for p = 1:numPerms
            permBootAccs(d, p) = computeAcc(compDataList, compLabelList, krDataList, ...
                krLabelList, true);
        end
        trueBootAcc = trueBootAccs(d);
        permBootAcc = permBootAccs(d,:);
        save(sprintf(drawAccFname, resRoot, d), 'trueBootAcc', 'permBootAcc');
    end
end
end

function [AUC, betas] = computeAcc(compDataList, compLabelList, krDataList, ...
    krLabelList, doPerm, varargin)

numSub = length(compDataList);
if nargout > 1
    betas = cell(numSub, 1);
end
if nargin > 5
    betas = varargin{1};
end
compProbs = [];
pooledLabels = [];
% For each subject, learn and apply the classifier
for s = 1:numSub
    if nargin > 5
        B = betas{s};
    else
        if doPerm
            numItems = length(compLabelList{s});
            B = logReg(compDataList{s}, compLabelList{s}(randPerm(numItems)), [1 2], false);
        else
            B = logReg(compDataList{s}, compLabelList{s}, [1 2], false);
        end
    end
    krData = krDataList{s};
    [numItems, numRounds, ~] = size(krData);
    roundProbs = nan(numItems, numRounds);
    for r = 1:numRounds
        roundProbs(:,r) = [squeeze(krData(:,r,:)) ones(numItems, 1)]*B;
    end
    
    compProbs = cat(1, compProbs, roundProbs);
    pooledLabels = cat(1, pooledLabels, krLabelList{s});
    
    if nargout > 1
        betas{s} = B;
    end
end

meanDiff = mean(compProbs(:,1:2), 2) - mean(compProbs(:,3:4), 2);

[~, ~, ~, AUC] = perfcurve(pooledLabels, meanDiff, 1);

end
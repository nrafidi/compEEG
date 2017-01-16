% trainPermCompClass_forKR(sub, krData, krLabels, compWin, krWin)

compWin = 17;
numPerms = 100;
subjects = {'AA', 'BB', 'DD', 'F', 'EE', 'GG', 'HH', 'JJ', ...
    'K', 'M', 'O', 'R', 'S', 'T', 'U', 'V', 'X', 'Y', 'Z', 'N'};
numSub = length(subjects);

addpath ./logisticRegression/
numChan = 64;
dataRoot = ['/Users/nrafidi/Documents/MATLAB/compEEG-data/preproc-final/' sub '/'];
fCPrefix = [dataRoot 'CompEEG_'];



for s = 1:numSub
    sub = subjects{s};
    fResPrefix = ['/Users/nrafidi/Documents/MATLAB/compEEG-data/results/' sub '/'];
    if ~exist(fResPrefix, 'dir')
        mkdir(fResPrefix);
    end
    fSuffix = '_Vis_BP2-200_N60_Ref_Epochs_Base_ICA1-2_Features_Overlap_Time.mat';
    compFname = [fCPrefix sub fSuffix];
    load(compFname);
    numSamp = size(featData, 1);
    featData = double(featData);
    [compData, mu, sigma] = zscore(featData);
    permB = cell(numPerm, 1);
    for p = 1:numPerm
        rng(p);
        compLabels = labels(randperm(numSamp),1); %#ok<*NODEF>
        
        compWinToUse = (compWin*numChan + 1):((compWin+1)*numChan);
        
        permB{p} = logReg(compData(:,compWinToUse), compLabels, [1 2], false);
        
        
    end
    
    save([fResPrefix, 'permB.mat'], 'permB');
end



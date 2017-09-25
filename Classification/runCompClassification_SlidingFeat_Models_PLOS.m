addpath ./logisticRegression/
addpath ../Preprocessing/

isRepExp = true;
if isRepExp
    subjects = {'III', 'JJJ', 'KKK', 'BBB', 'GGG', 'HHH', 'AAA', 'CCC', 'DDD', 'EEE', 'FFF', 'MM', ...
        'OO', 'PP', 'QQ', 'RR', 'SS', 'TT', 'WW',...
        'YY'};
    fPrefix = '/Users/nrafidi/Documents/MATLAB/compEEG-data-rep/preproc-final/';
    resDir = '/Users/nrafidi/Documents/MATLAB/compEEG-data-rep/results/';
else
    subjects =  {'L', 'P', 'W', 'CC', 'FF', 'H', 'I', 'J', ...
        'AA', 'BB', 'DD', 'EE', 'GG', 'HH', 'JJ', 'M', 'N',...
        'F', 'K', 'O', 'R', 'S', 'T', 'U', 'V', 'X', 'Y', 'Z'};
    fPrefix = '/Users/nrafidi/Documents/MATLAB/compEEG-data/preproc-final/';
    resDir = '/Users/nrafidi/Documents/MATLAB/compEEG-data/results/';
end
numSub = length(subjects);
numFeatWins = 54;
winSizeOptions = 50;%:50:200;
doZ = 1;

fSuffix = '_Vis_BP2-200_N60_Ref_Epochs_Base_ICA1-2_Features_Overlap_Time.mat';

subModels = cell(numSub, numFeatWins);

rng('shuffle');

for s = 1:length(subjects)
    
    sub = subjects{s};
    loadFname = [fPrefix sub '/CompEEG_' sub fSuffix];
    if ~exist(loadFname, 'file')
        fprintf('Subject %s failed.\n', sub);
    else
        load(loadFname);
        
        Y = labels(:,1);
        [N, W] = size(featData);
        numFeatWins = W/64;
        
        for p = 1:numFeatWins
            startInd = (p-1)*64 + 1;
            endInd = p*64;
            X = featData(:,startInd:endInd);
            
            if doZ
                [trainData, mu, sigma] = zscore(X);
            else
                trainData = X; %#ok<*UNRCH>
            end
            
            [B, ~] = logReg(trainData, Y, [1 2], false);
            
            
            subModels{s, p} = B;
            
        end
        subModel = subModels(s,:);
        if ~exist([resDir sub '/'], 'dir')
            mkdir([resDir sub '/']);
        end
        save([resDir sub '/CompEEG_CV_Slide_Models_PLOS.mat'],...
            'subModel');
        fprintf('Subject %s succeeded.\n', sub);
    end
end

save([resDir 'CompEEG_CV_Slide_Models_PLOS.mat'],...
    'subModels');
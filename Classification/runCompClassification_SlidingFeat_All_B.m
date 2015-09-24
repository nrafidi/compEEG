% Runs cross-validated classification within subject for the competition
% data

addpath ./logisticRegression/
addpath ../Preprocessing/

subjects = {'L', 'P', 'W', 'CC', 'FF', 'H', 'I', 'J', ...
    'AA', 'BB', 'DD', 'EE', 'GG', 'HH', 'JJ', 'M', 'N',...
    'F', 'K', 'O', 'R', 'S', 'T', 'U', 'V', 'X', 'Y', 'Z'};
numSub = length(subjects);
numFeatWins = 40;%#magic
winSizeOptions = 50;%:50:200;
doZ = 1;
fPrefix = '/Users/nrafidi/Documents/MATLAB/compEEG-data/preproc-final/';
fSuffix = '_Vis_BP2-200_N60_Ref_Epochs_Base_ICA1-2.mat';

subModels = cell(numSub, numFeatWins);

rng('shuffle');

for s = 1:length(subjects)
    
    sub = subjects{s};
    loadFname = [fPrefix sub '/CompEEG_' sub fSuffix];
    if ~exist(loadFname, 'file')
        fprintf('Subject %s failed.\n', sub);
    else
        options = struct;
        options.overLap = true;
        options.erpWinSize = winSizeOptions;
        [featData, labels, ~] = extractFeatures(loadFname, options);
        
        
        Y = labels(:,1);
        N = size(featData, 1);
        
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
        if ~exist(['/Users/nrafidi/Documents/MATLAB/compEEG-data/results/' sub '/'], 'dir')
            mkdir(['/Users/nrafidi/Documents/MATLAB/compEEG-data/results/' sub '/']);
        end
        save(['/Users/nrafidi/Documents/MATLAB/compEEG-data/results/' sub '/CompEEG_CV_Slide_Models.mat'],...
            'subModel');
        fprintf('Subject %s succeeded.\n', sub);
    end
end

save('/Users/nrafidi/Documents/MATLAB/compEEG-data/results/CompEEG_CV_Slide_Models.mat',...
    'subModels');
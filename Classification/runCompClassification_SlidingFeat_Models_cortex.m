function runCompClassification_SlidingFeat_Models_cortex(sub)

addpath ~/compEEG/Classification/logisticRegression/

fPrefix = '~/CompEEG/Data/';
resDir = '~/CompEEG/Results/';
numFeatWins = 54;
doZ = 1;

fSuffix = '_Vis_BP2-200_N60_Ref_Hilbert-theta_Epochs_Features_Overlap_Time.mat';

models = cell(numFeatWins);
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
        
        
        models{p} = B;
        
    end
    if ~exist([resDir sub '/'], 'dir')
        mkdir([resDir sub '/']);
    end
end


save([resDir 'CompEEG_' sub '_CV_Slide_Models_PLOS_theta.mat'],...
    'models');
end
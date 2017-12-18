% Runs cross-validated classification within subject for the competition
% data

isRepExp = true;

addpath ./logisticRegression/
addpath ../Preprocessing/

if isRepExp
    subjects = {'III', 'JJJ', 'KKK', 'BBB', 'GGG', 'HHH', 'AAA', 'CCC', 'DDD', 'EEE', 'FFF', 'MM', ...
        'OO', 'PP', 'QQ', 'RR', 'SS', 'TT', 'WW',...
        'YY'};
else
    subjects =  {'L', 'P', 'W', 'CC', 'FF', 'H', 'I', 'J', ...
        'AA', 'BB', 'DD', 'EE', 'GG', 'HH', 'JJ', 'M', 'N',...
        'F', 'K', 'O', 'R', 'S', 'T', 'U', 'V', 'X', 'Y', 'Z'};
end
numSub = length(subjects);

numPerm = 100;


doZ = 1;
if isRepExp
    fPrefix = '/Users/nrafidi/Documents/MATLAB/compEEG-data-rep/preproc-final/';
else
    fPrefix = '/Users/nrafidi/Documents/MATLAB/compEEG-data/preproc-final/';
end
fSuffix = '_Vis_BP2-%d_N60_Ref_Epochs_Base_ICA1-2_Features_Overlap_Time.mat';
numFolds = 5;
trueSubAccs = nan(numSub, numFolds);
permSubAccs = nan(numSub, numPerm, numFolds);

load comp5Fseed.mat
%%
numSamp = nan(length(subjects), 1);
    for s = 1:length(subjects)
        sub = subjects{s};
        
        
            if strcmp(sub, 'NN')
                fSuffix = sprintf(fSuffix, 127);
            else
                fSuffix = sprintf(fSuffix, 200);
            end
            loadFname = [fPrefix sub '/CompEEG_' sub fSuffix];
            if ~exist(loadFname, 'file')
%                 keyboard;
                fprintf('Subject %s failed.\n', sub);
                continue
            else
                load(loadFname);
                
                
                
                Y = labels(:,1);
                [N, W] = size(featData);
                numSamp(s) = N;
            end
    end
disp(mean(numSamp))
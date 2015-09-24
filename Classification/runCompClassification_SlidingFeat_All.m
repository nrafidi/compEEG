% Runs cross-validated classification within subject for the competition
% data

addpath ./logisticRegression/
addpath ../Preprocessing/

subjects = {'L', 'P', 'W', 'CC', 'FF', 'H', 'I', 'J', ...
    'AA', 'BB', 'DD', 'EE', 'GG', 'HH', 'JJ', 'M', 'N',...
    'F', 'K', 'O', 'R', 'S', 'T', 'U', 'V', 'X', 'Y', 'Z'};
numSub = length(subjects);
numFeatWins = 40;%magic
winSizeOptions = 50:50:200;
numSizes = length(winSizeOptions);

freqName = 'beta';

doZ = 1;
fPrefix = '/Users/nrafidi/Documents/MATLAB/compEEG-data/preproc-final/';
fSuffix = ['_Vis_Hilbert-' freqName '_Epochs.mat'];
subAccs = nan(numSub, numFeatWins, numSizes);
subModels = cell(numSub, numFeatWins, numSizes);

rng('shuffle');

for s = 1:length(subjects)
    
    sub = subjects{s};
    loadFname = [fPrefix sub '/CompEEG_' sub fSuffix];
    if ~exist(loadFname, 'file')
        fprintf('Subject %s failed.\n', sub);
    else
        options = struct;
        options.overLap = true;
        for w = 1:numSizes
            options.erpWinSize = winSizeOptions(w);
            options.specWinSize = winSizeOptions(w);
            [featData, labels, ~] = extractFeatures(loadFname, options);
            
            
            Y = labels(:,1);
            N = size(featData, 1);
            
            trainingMark = floor(0.85*N);
            trainInds = 1:trainingMark;
            testInds = (trainingMark+1):N;
            for p = 1:numFeatWins
                startInd = (p-1)*64 + 1;
                endInd = p*64;
                X = featData(:,startInd:endInd);
                
                if doZ
                    [trainData, mu, sigma] = zscore(X(trainInds, :));
                    testData = (X(testInds,:) - repmat(mu, length(testInds), 1))./repmat(sigma,length(testInds), 1);
                else
                    trainData = X(trainInds, :); %#ok<*UNRCH>
                    testData = X(testInds,:);
                end
                
                [B, lambda] = logReg(trainData, Y(trainInds), [1 2], false);
                
                P = [testData ones(length(testInds), 1)]*B;
                P = exp(P)./(1 + exp(P));
                Yhat = double(P > 0.5);
                
                subAccs(s, p, w) = sum(Yhat == Y(testInds))/length(Yhat);
                subModels{s, p, w} = B;
            end
            
        end
        subAcc = squeeze(subAccs(s,:,:));
        save(['/Users/nrafidi/Documents/MATLAB/compEEG-data/results/CompEEG_CV_' sub '_Hilbert-' freqName '_Slide_All.mat'],...
            'subAcc');
        fprintf('Subject %s succeeded.\n', sub);
        
        subModel = subModels(s,:);
        if ~exist(['/Users/nrafidi/Documents/MATLAB/compEEG-data/results/' sub '/'], 'dir')
            mkdir(['/Users/nrafidi/Documents/MATLAB/compEEG-data/results/' sub '/']);
        end
        save(['/Users/nrafidi/Documents/MATLAB/compEEG-data/results/' sub '/CompEEG_CV_Slide_Hilbert-' freqName '_Models.mat'],...
            'subModel');
    end
end
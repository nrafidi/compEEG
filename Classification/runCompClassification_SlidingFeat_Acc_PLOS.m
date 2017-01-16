% Runs cross-validated classification within subject for the competition
% data

addpath ./logisticRegression/
addpath ../Preprocessing/

subjects = {'L', 'P', 'W', 'CC', 'FF', 'H', 'I', 'J', ...
    'AA', 'BB', 'DD', 'EE', 'GG', 'HH', 'JJ', 'M', 'N',...
    'F', 'K', 'O', 'R', 'S', 'T', 'U', 'V', 'X', 'Y', 'Z'};
numSub = length(subjects);
numFeatWins = 55;
numFolds = 5;
doZ = 1;
fPrefix = '/Users/nrafidi/Documents/MATLAB/compEEG-data/preproc-final/';
fSuffix = '_Vis_BP2-200_N60_Ref_Epochs_Base_ICA1-2_ItemAvg_Features_Overlap_Time.mat';
subAccs = nan(numSub, numFeatWins, numFolds);

load comp5Fseed.mat
% rng('shuffle');
%%
for s = 1:length(subjects)
    tic
    sub = subjects{s};
    loadFname = [fPrefix sub '/CompEEG_' sub fSuffix];
    saveFname = ['/Users/nrafidi/Documents/MATLAB/compEEG-data/results/' sub '/CompEEG_5FCV-rsh_Slide_Accs.mat'];
    
    if exist(saveFname, 'file')
        load(saveFname);
        subAccs(s,:,:) = subAcc;
        fprintf('Subject %s Loaded.\n', sub);
    else
        if ~exist(loadFname, 'file')
            fprintf('Subject %s failed.\n', sub);
            continue
        else
            load(loadFname);
            
            Y = labels(:,1);
            [N, W] = size(featData);
            numFeatWins = W/64;
            for p = 1:numFeatWins
                startInd = (p-1)*64 + 1;
                endInd = p*64;
                X = featData(:,startInd:endInd);
                rng(comp5Fseed);
                folds = crossvalind('Kfold', N, numFolds);
                for f = 1:numFolds
                    testInds = (folds == f);
                    if doZ
                        [trainData, mu, sigma] = zscore(X(~testInds, :));
                        testData = (X(testInds,:) - repmat(mu, sum(testInds), 1))./repmat(sigma,sum(testInds), 1);
                    else
                        trainData = X(~testInds, :); %#ok<*UNRCH>
                        testData = X(testInds,:);
                    end
                    if f > 1
                        [B, lambda] = logReg(trainData, Y(~testInds), [1 2], false, B);
                    else
                        [B, lambda] = logReg(trainData, Y(~testInds), [1 2], false);
                    end
                    P = [testData ones(sum(testInds), 1)]*B;
                    P = exp(P)./(1 + exp(P));
                    Yhat = double(P > 0.5);
                    
                    subAccs(s, p, f) = sum(Yhat == Y(testInds))/length(Yhat);
                end
            end
            subAcc = squeeze(subAccs(s,:,:));
            save(saveFname, 'subAcc');
            fprintf('Subject %s succeeded.\n', sub);
        end
        toc
    end
end
save('/Users/nrafidi/Documents/MATLAB/compEEG-data/results/CompEEG_5FCV-rsh_Slide_Accs.mat',...
    'subAccs', 'winTime');
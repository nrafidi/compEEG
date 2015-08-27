% Runs cross-validated classification within subject for the competition
% data

addpath ./logisticRegression/
addpath ../Preprocessing/

subjects = {'L', 'P', 'W', 'CC', 'FF'};%H I J
%{'AA', 'BB', 'DD', 'EE', 'F', 'G', 'GG', 'HH', 'JJ', ...
%     'K', 'M', 'N', 'O', 'R', 'S', 'T', 'U', 'V', 'X', 'Y', 'Z'};
numSub = length(subjects);
numFolds = 5;
numFeatWins = 28;
doZ = 1;
fPrefix = '/Users/nrafidi/Documents/MATLAB/compEEG-data/preproc-final/';
fSuffix = '_Vis_BP2-200_N60_Ref_Epochs_Base_ICA1-2_Features.mat';

subAccs = nan(numSub, numFeatWins, numFolds);

rng('shuffle');

for s = 1:length(subjects)
    
    sub = subjects{s};
    if ~exist([fPrefix sub '/CompEEG_' sub fSuffix], 'file')
        [featData, labels] = preprocPipeline(sub, 'CompEEG');
    else
        load([fPrefix sub '/CompEEG_' sub fSuffix]);
    end
    
    
    Y = labels(:,1);
    N = size(featData, 1);
    for p = 1:numFeatWins
        startInd = (p-1)*64 + 1;
        endInd = p*64;
        X = featData(:,startInd:endInd);
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
%             tic
            if f > 1
                [B, lambda] = logReg(trainData, Y(~testInds), [1 2], false, B);
            else
                [B, lambda] = logReg(trainData, Y(~testInds), [1 2], false);
            end
%             toc
%             disp(lambda);
            P = [testData ones(sum(testInds), 1)]*B;
            P = exp(P)./(1 + exp(P));
            Yhat = double(P > 0.5);
            
            subAccs(s, p, f) = sum(Yhat == Y(testInds))/length(Yhat);
            
        end
        disp(squeeze(mean(subAccs(s, p, :), 3)));
    end
    subAcc = squeeze(subAccs(s,:,:));
    
    
    save(['/Users/nrafidi/Documents/MATLAB/compEEG-data/results/CompEEG_CV_' sub '_Def_Slide.mat'],...
        'subAcc');
    
end

% disp(squeeze(mean(subAccs, 3)));
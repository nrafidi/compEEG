% Runs cross-validated classification within subject for the competition
% data

addpath ./logisticRegression/
addpath ../Preprocessing/

subjects = {'L', 'P', 'W', 'CC', 'FF', 'H', 'I', 'J', ...
    'AA', 'BB', 'DD', 'EE', 'GG', 'HH', 'JJ', 'M', 'N',...
    'F', 'K', 'O', 'R', 'S', 'T', 'U', 'V', 'X', 'Y', 'Z'};
numSub = length(subjects);
numFeatWins = 54;
winSizeOptions = 50;

doZ = 1;
fPrefix = '/Users/nrafidi/Documents/MATLAB/compEEG-data/preproc-final/';
fSuffix = '_Vis_BP2-200_N60_Ref_Epochs_Base_ICA1-2_Features_Overlap_Time.mat';
subAccs = nan(numSub, numFeatWins, numFeatWins);

rng('shuffle');
%%
for s = 1:length(subjects)
    tic
    sub = subjects{s};
    loadFname = [fPrefix sub '/CompEEG_' sub fSuffix];
    resFname = ['/Users/nrafidi/Documents/MATLAB/compEEG-data/results/' sub '/CompEEG_CV_Slide_TGM_Accs.mat'];
    if ~exist(loadFname, 'file')
        fprintf('Subject %s failed.\n', sub);
        continue
    elseif exist(resFname, 'file')
        load(resFname);
        subAccs(s,:,:) = subAcc;
        fprintf('Subject %s loaded.\n', sub);
    else
        load(loadFname);
        
        Y = labels(:,1);
        [N, W] = size(featData);
        numFeatWins = W/64;
        rng('shuffle');
        labelShuff = randperm(N);
        featData = featData(labelShuff,:);
        Y = Y(labelShuff);
        trainingMark = floor(0.85*N);
        trainInds = 1:trainingMark;
        testInds = (trainingMark+1):N;
        h = waitbar(0, sub);
        for p = 1:numFeatWins
            startIndTrain = (p-1)*64 + 1;
            endIndTrain = p*64;
            for pp = 1:numFeatWins
                waitbar(((p-1)*numFeatWins + pp)/(numFeatWins^2), h);
                startIndTest = (pp-1)*64 + 1;
                endIndTest = pp*64;
                if doZ
                    [trainData, mu, sigma] = zscore(featData(trainInds, startIndTrain:endIndTrain));
                    testData = (featData(testInds,startIndTest:endIndTest) - repmat(mu, length(testInds), 1))./repmat(sigma,length(testInds), 1);
                else
                    trainData = featData(trainInds, startIndTrain:endIndTrain); %#ok<*UNRCH>
                    testData = featDatac(testInds,startIndTest:endIndTest);
                end
                
                [B, lambda] = logReg(trainData, Y(trainInds), [1 2], false);
                
                P = [testData ones(length(testInds), 1)]*B;
                P = exp(P)./(1 + exp(P));
                Yhat = double(P > 0.5);
                
                subAccs(s, p, pp) = sum(Yhat == Y(testInds))/length(Yhat);
            end
        end
        close(h);
        subAcc = squeeze(subAccs(s,:,:));
        save(resFname,'subAcc');
        fprintf('Subject %s succeeded.\n', sub);
        toc
    end
    
end

save('/Users/nrafidi/Documents/MATLAB/compEEG-data/results/CompEEG_CV_Slide_TGM_Accs.mat',...
    'subAccs', 'winTime');
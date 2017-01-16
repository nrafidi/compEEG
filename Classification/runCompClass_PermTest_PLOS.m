% Runs cross-validated classification within subject for the competition
% data

addpath ./logisticRegression/
addpath ../Preprocessing/

subjects = {'L', 'P', 'W', 'CC', 'FF', 'H', 'I', 'J', ...
    'AA', 'BB', 'DD', 'EE', 'GG', 'HH', 'JJ', 'M', 'N',...
    'F', 'K', 'O', 'R', 'S', 'T', 'U', 'V', 'X', 'Y', 'Z'};
numSub = length(subjects);

numPerm = 100;
winToUse = 17;

doZ = 1;
fPrefix = '/Users/nrafidi/Documents/MATLAB/compEEG-data/preproc-final/';
fSuffix = '_Vis_BP2-200_N60_Ref_Epochs_Base_ICA1-2_Features_Overlap_Time.mat';
numFolds = 5;
trueSubAccs = nan(numSub, numFolds);
permSubAccs = nan(numSub, numPerm, numFolds);

load comp5Fseed.mat
%%

for s = 1:length(subjects)
    sub = subjects{s};
    resultFname = ['/Users/nrafidi/Documents/MATLAB/compEEG-data/results/' ...
        sub '/CompEEG_5FCV_win' num2str(winToUse) '_permAccs.mat'];
    if exist(resultFname, 'file');
        trueSubAccs(s,:) = trueAcc;
        permSubAccs(s,:,:) = permAccs;
        fprintf('Subject %s Loaded\n', sub);
    else
        
        loadFname = [fPrefix sub '/CompEEG_' sub fSuffix];
        if ~exist(loadFname, 'file')
            fprintf('Subject %s failed.\n', sub);
            continue
        else
            load(loadFname);
            
            Y = labels(:,1);
            [N, W] = size(featData);
            numFeatWins = W/64;
            % Data Splitting
            % True acc for this splitting of the data
            startInd = (winToUse-1)*64 + 1;
            endInd = winToUse*64;
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
                
                trueSubAccs(s, f) = sum(Yhat == Y(testInds))/length(Yhat);
            end
        end
        %         h = waitbar(0, sub);
        tic
        for p = 1:numPerm
            %             waitbar(p/numPerm, h);
            %Now for the permutation
            rng('shuffle');
            Y = Y(randperm(N));
            
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
                
                permSubAccs(s, p, f) = sum(Yhat == Y(testInds))/length(Yhat);
            end
        end
        trueAcc = trueSubAccs(s,:);
        permAccs = permSubAccs(s,:,:);
        save(resultFname,...
            'trueAcc', 'permAccs');
        fprintf('Subject %s succeeded.\n', sub);
        
        toc
    end
    %     close(h)
end
save(['/Users/nrafidi/Documents/MATLAB/compEEG-data/results/CompEEG_5FCV_win' ...
    num2str(winToUse) '_permAccs.mat'],...
    'trueSubAccs', 'permSubAccs', 'winTime');
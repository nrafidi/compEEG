% Runs cross-validated classification within subject for the competition
% data

function runCompClass_PermTest_cortex(sub, winToUse)

addpath ~/compEEG/Classification/logisticRegression/

% loadFname = sprintf('~/CompEEG/Data/CompEEG_%s_Vis_BP2-200_N60_Ref_Hilbert-theta_Epochs_Features_Overlap_Time.mat', sub);
loadFname = sprintf('~/CompEEG/Data/CompEEG_%s_Vis_BP2-200_N60_Ref_Epochs_Base_ICA1-2_Features_Overlap_Time.mat', sub);
resultFname = sprintf('~/CompEEG/Results/CompEEG_%s_5FCV_win%d_permAccs100.mat', sub, winToUse);
numPerm = 100;
doZ = 1;
numFolds = 5;

load comp5Fseed.mat
%%
if ~exist(resultFname, 'file');
    if ~exist(loadFname, 'file')
        fprintf('Subject %s failed.\n', sub);
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
        disp(size(X, 1))
        rng(comp5Fseed);
        folds = crossvalind('Kfold', N, numFolds);
        trueAcc = nan(numFolds, 1);
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
            
            trueAcc(f) = sum(Yhat == Y(testInds))/length(Yhat);
        end
    end
    %             h = waitbar(0, sub);
    permAccs = nan(numPerm, numFolds);
    tic
    for p = 1:numPerm
        %                 waitbar(p/numPerm, h);
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
            
            permAccs(p, f) = sum(Yhat == Y(testInds))/length(Yhat);
        end
    end
    save(resultFname,...
        'trueAcc', 'permAccs', 'winTime');
    fprintf('Subject %s succeeded.\n', sub);
    
    toc
    %             close(h)
end

end
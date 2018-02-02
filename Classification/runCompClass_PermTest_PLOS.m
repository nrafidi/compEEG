% Runs cross-validated classification within subject for the competition
% data

isRepExp = false;

addpath ./logisticRegression/
addpath ../Preprocessing/

if isRepExp
    subjects = {'OO', 'III', 'JJJ', 'KKK', 'BBB', 'GGG', 'HHH', 'AAA', 'CCC', 'DDD', 'EEE', 'FFF', 'MM', ...
        'PP', 'QQ', 'RR', 'SS', 'TT', 'WW',...
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
fSuffix = '_Vis_BP2-%d_N60_Ref_Hilbert-theta_Epochs_Features_Overlap_Time.mat';
numFolds = 5;
trueSubAccs = nan(numSub, numFolds);
permSubAccs = nan(numSub, numPerm, numFolds);

load comp5Fseed.mat
%%
for winToUse = 1:40 %[1:6, 41:54]
    for s = 1:length(subjects)
        sub = subjects{s};
        
        if isRepExp
            resultDir = ['/Users/nrafidi/Documents/MATLAB/compEEG-data-rep/results/' sub '/'];
        else
            resultDir = ['/Users/nrafidi/Documents/MATLAB/compEEG-data/results/' sub '/'];
        end
        if ~exist(resultDir, 'dir')
            mkdir(resultDir);
        end
        
        resultFname = [resultDir '/CompEEG_5FCV_win' num2str(winToUse) '_permAccs.mat'];
        if exist(resultFname, 'file');
            load(resultFname);
            trueSubAccs(s,:) = trueAcc;
            permSubAccs(s,:,:) = permAccs;
            
            if s == length(subjects)
                if strcmp(sub, 'NN')
                    fSuffix = sprintf(fSuffix, 127);
                else
                    fSuffix = sprintf(fSuffix, 200);
                end
                loadFname = [fPrefix sub '/CompEEG_' sub fSuffix];
                load(loadFname);
            end
            fprintf('Subject %s Loaded\n', sub);
        else
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
                numFeatWins = W/64;
                % Data Splitting
                % True acc for this splitting of the data
                startInd = (winToUse-1)*64 + 1;
                endInd = winToUse*64;
                X = featData(:,startInd:endInd);
                disp(size(X, 1))
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
%             h = waitbar(0, sub);
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
                    
                    permSubAccs(s, p, f) = sum(Yhat == Y(testInds))/length(Yhat);
                end
            end
            trueAcc = trueSubAccs(s,:);
            permAccs = permSubAccs(s,:,:);
            save(resultFname,...
                'trueAcc', 'permAccs', 'winTime');
            fprintf('Subject %s succeeded.\n', sub);
            
            toc
%             close(h)
        end
        
    end
    if isRepExp
        save(['/Users/nrafidi/Documents/MATLAB/compEEG-data-rep/results/CompEEG_5FCV_win' ...
            num2str(winToUse) '_permAccs.mat'],...
            'trueSubAccs', 'permSubAccs', 'winTime');
    else
        save(['/Users/nrafidi/Documents/MATLAB/compEEG-data/results/CompEEG_5FCV_win' ...
            num2str(winToUse) '_permAccs.mat'],...
            'trueSubAccs', 'permSubAccs', 'winTime');
    end
end
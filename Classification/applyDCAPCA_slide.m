addpath logisticRegression/
dataRoot = '/Users/nrafidi/Documents/MATLAB/compEEG-data/preproc-final/';
resultRoot = '/Users/nrafidi/Documents/MATLAB/compEEG-data/results/';

type = 'centered_alltimes';

numT = 44;
numDim = 128;
numCompToUse = 10;

load(sprintf('%sallSubData_KR_5fold_withLabels.mat', dataRoot));
for fold = 1:5
    addpath(sprintf('%sfold%d/', dataRoot, fold));
    load(sprintf('S_dca_pca_%s_fold%d.mat', type, fold));
    numSub = length(S);
    meanToSub = nan(numSub, numDim);
    stdToDiv = nan(numSub, numDim);
    for s = 1:numSub
        if isfield(S(s), 'mean')
            meanToSub(s,:) = [S(s).mean; S(s).mean]';
        else
            meanToSub(s,:) = zeros(1, numDim);
        end
        if isfield(S(s), 'std')
            stdToDiv(s,:) = [S(s).std; S(s).std]';
        else
            stdToDiv(s,:) = ones(1, numDim);
        end
        
    end
%     acc_Reg = nan(numT, 1);
    acc_DCA = nan(numT, 1);
    acc_PCA = nan(numT, 1);
    
    testLabels = cell2mat(allSubTestLabels(:,fold));
    trainLabels = cell2mat(allSubTrainLabels(:,fold));
    
    numTestSamples = length(testLabels);
    h = waitbar(0, sprintf('Fold %d', fold));
    for t = 1:numT
        %     for c = 1:numComp
        waitbar(t/numT, h);
        allSubTrainData_DCA = cell(numSub, 1);
        allSubTestData_DCA = cell(numSub, 1);
        allSubTrainData_PCA = cell(numSub, 1);
        allSubTestData_PCA = cell(numSub, 1);
        for s = 1:numSub
            subTrainData = bsxfun(@rdivide, ...
                bsxfun(@minus, allSubTrainData{s, t, fold}, meanToSub(s,:)), ...
                stdToDiv(s,:)); %#ok<*SAGROW>
            subTestData = bsxfun(@rdivide, ...
                bsxfun(@minus, allSubTestData{s, t, fold}, meanToSub(s,:)), ...
                stdToDiv(s,:));
            allSubTrainData_DCA{s} = subTrainData*[S(s).U_dca(:, 1:numCompToUse);S(s).U_dca(:, 1:numCompToUse)];
            allSubTestData_DCA{s} = subTestData*[S(s).U_dca(:, 1:numCompToUse);S(s).U_dca(:, 1:numCompToUse)];
            allSubTrainData_PCA{s} = subTrainData*[S(s).U_pca(:, 1:numCompToUse);S(s).U_pca(:, 1:numCompToUse)];
            allSubTestData_PCA{s} = subTestData*[S(s).U_pca(:, 1:numCompToUse);S(s).U_pca(:, 1:numCompToUse)];
        end
        
        
        trainData = cell2mat(allSubTrainData(:,t,fold));
        testData = cell2mat(allSubTestData(:,t,fold));
        
        dca_trainData = cell2mat(allSubTrainData_DCA);
        dca_testData = cell2mat(allSubTestData_DCA);
        %         [dca_trainData, meanDCA, sigDCA] = zscore(trainData*U_dca);
        %         dca_testData = bsxfun(@rdivide, bsxfun(@minus, testData*U_dca, meanDCA), sigDCA);
        w_DCA = logReg(dca_trainData, trainLabels, 0, 0);
        
        P_DCA = [dca_testData ones(numTestSamples, 1)]*w_DCA;
        acc_DCA(t) = sum(double(P_DCA > 0) == testLabels)/numTestSamples;
        fprintf('DCA accuracy at %d = %d\n', t, acc_DCA(t));
        
        pca_trainData = cell2mat(allSubTrainData_PCA);
        pca_testData = cell2mat(allSubTestData_PCA);
        %         [pca_trainData, meanPCA, sigPCA] = zscore(trainData*U_pca);
        %         pca_testData = bsxfun(@rdivide, bsxfun(@minus, testData*U_pca, meanPCA), sigPCA);
        [w_PCA, lambda] = logReg(pca_trainData, trainLabels, 0, 0);
        %         disp(lambda);
        P_PCA = [pca_testData ones(numTestSamples, 1)]*w_PCA;
        acc_PCA(t) = sum(double(P_PCA > 0) == testLabels)/numTestSamples;
        fprintf('PCA accuracy at %d = %d\n', t, acc_PCA(t));
%         
%         [trainData, mu, sig] = zscore(trainData);
%         testData = bsxfun(@rdivide, bsxfun(@minus, testData, mu), sig);
%         w_reg = logReg(trainData, trainLabels, [1 2], 0);
%         P_reg = [testData ones(size(testData, 1), 1)]*w_reg;
%         acc_Reg(t) = sum(double(P_reg > 0) == testLabels)/numTestSamples;
%         fprintf('Reg accuracy at %d = %d\n', t, acc_Reg(t));
        %     end
    end
    
    close(h);
    
    save(sprintf('%sKRClassification_PCA-DCA_fold%d_%s_regLambda_numComp%d.mat', ...
        resultRoot, fold, type, numCompToUse), 'acc_DCA', 'acc_PCA');
end
addpath logisticRegression/
dataRoot = '/Users/nrafidi/Documents/MATLAB/compEEG-data/preproc-final/';
resultRoot = '/Users/nrafidi/Documents/MATLAB/compEEG-data/results/';

numT = 44;
% numComp = 9;
numFolds = 5;

load(sprintf('%sallSubData_KR_5fold_withPassLabels.mat', dataRoot));
for fold = 1:5
    addpath(sprintf('%sfold%d/', dataRoot, fold));
    
    acc_Reg = nan(numT, 1);
    acc_DCA = nan(numT, 1);%numComp);
    acc_PCA = nan(numT, 1);%numComp);
    
    testLabels = cell2mat(allSubTestLabels(:,fold));
    trainLabels = cell2mat(allSubTrainLabels(:,fold));
    
    numTestSamples = length(testLabels);
    h = waitbar(0, sprintf('Fold %d', fold));
    for t = 1:numT
        %     for c = 1:numComp
        waitbar(t/numT, h);
        load(sprintf('U_dca_pca_timepoint%d.mat', t));
        trainData = cell2mat(allSubTrainData(:,t,fold));
        testData = cell2mat(allSubTestData(:,t,fold));
        
        dca_trainData = trainData*U_dca;
        dca_testData = testData*U_dca;
%         [dca_trainData, meanDCA, sigDCA] = zscore(trainData*U_dca);
%         dca_testData = bsxfun(@rdivide, bsxfun(@minus, testData*U_dca, meanDCA), sigDCA);
        [w_DCA, lambda] = logReg(dca_trainData, trainLabels, 0, 0);
%         disp(lambda);
        P_DCA = [dca_testData ones(numTestSamples, 1)]*w_DCA;
        acc_DCA(t) = sum(double(P_DCA > 0) == testLabels)/numTestSamples;
        fprintf('DCA accuracy at %d = %d\n', t, acc_DCA(t));
        
        pca_trainData = trainData*U_pca;
        pca_testData = testData*U_pca;
%         [pca_trainData, meanPCA, sigPCA] = zscore(trainData*U_pca);
%         pca_testData = bsxfun(@rdivide, bsxfun(@minus, testData*U_pca, meanPCA), sigPCA);
        [w_PCA, lambda] = logReg(pca_trainData, trainLabels, 0, 0);
%         disp(lambda);
        P_PCA = [pca_testData ones(numTestSamples, 1)]*w_PCA;
        acc_PCA(t) = sum(double(P_PCA > 0) == testLabels)/numTestSamples;
        fprintf('PCA accuracy at %d = %d\n', t, acc_PCA(t));
        
        [trainData, mu, sig] = zscore(trainData);
        testData = bsxfun(@rdivide, bsxfun(@minus, testData, mu), sig);
        w_reg = logReg(trainData, trainLabels, [1 2], 0);
        P_reg = [testData ones(size(testData, 1), 1)]*w_reg;
        acc_Reg(t) = sum(double(P_reg > 0) == testLabels)/numTestSamples;
        fprintf('Reg accuracy at %d = %d\n', t, acc_Reg(t));
        %     end
    end
    
    close(h);
    
    save(sprintf('%sKRClassification_Reg-PCA-DCA_fold%d_regZ_regLambda.mat', resultRoot, fold), 'acc_Reg', 'acc_DCA', 'acc_PCA');
end
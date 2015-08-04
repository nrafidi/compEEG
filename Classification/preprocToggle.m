% Script for trying different combinations of preprocessing steps to see
% which ones matter for realtime analysis

clear
clc
eeglab;
addpath ./logisticRegression/
addpath ../Preprocessing/
subjects = {'H', 'I', 'J', 'L', 'P', 'W', 'CC', 'FF'};
experi = 'CompEEG';

% Toggle visual inspection, filters (2),  blink removal
results = nan(length(subjects), 16);
options = struct;
% Iterate through subjects
for s = 1:length(subjects)
    sub = subjects{s};
    
    for doVis = 0:1
        options.isVis = doVis;
        for doBP = 0:1
            if doBP
                options.HP = 2;
                options.LP = 200;
            else
                options.HP = nan;
            end
            for doN = 0:1
                if doN
                    options.N = 60;
                else
                    options.N = nan;
                end
                for runICA = 0:1
                    options.runICA = runICA;
                    
                    [featdata, labels] = ...
                        preprocPipeline(sub, experi, options);
                    
                    trainData = featdata(1:2:end,:);
                    trainLabels = labels(1:2:end, 1);
                    testData = featdata(2:2:end,:);
                    testLabels = labels(2:2:end, 1);
                    
                    
                    [trainData, mu, sigma] = zscore(trainData);
                    testData = (testData - ...
                        repmat(mu, size(testData,1), 1))./repmat(sigma,size(testData,1), 1);
                    
                    
                    [B, lambda] = logReg(trainData, trainLabels, [1 2], false);
                    
                    P = [testData ones(size(testData,1), 1)]*B;
                    P = exp(P)./(1 + exp(P));
                    Yhat = double(P > 0.5);
                    
                    ind = 8*doVis + 4*doBP + 2*doN + runICA + 1;
                    
                    results(s, ind) = sum(Yhat == testLabels)/length(Yhat);
                    
                end
            end
        end
        
    end
end

save ../../compEEG-data/results/preprocToggle.mat results subjects
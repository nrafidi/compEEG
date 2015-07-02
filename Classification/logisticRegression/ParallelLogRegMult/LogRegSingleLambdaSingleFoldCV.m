function [Weights, cvErr] = ...
    LogRegSingleLambdaSingleFoldCV( ...
        Features, ...
        Labels, ...
        lambda, ...
        numFolds, ...
        indexFold, ...
        regBias)
    
    folds = crossvalind('Kfold', size(Features, 1), numFolds);
    
    TrainFeatures = Features(folds ~= indexFold, :);
    TestFeatures = Features(folds == indexFold, :);
    
    TrainLabels = Labels(folds ~= indexFold, :);
    TestLabels = Labels(folds == indexFold, :);
    
    Weights = GradDesc(TrainFeatures, TrainLabels, lambda, regBias);
    %Something wrong with SAG - doesn't currently converge
    %Weights = StochasticAverageGradient(TrainFeatures, TrainLabels, lambda, regBias);
    PredictedLabels = LogPredict(TestFeatures, Weights, regBias);
    Predictedhat = zeros(size(PredictedLabels));
    Predictedhat(PredictedLabels >= .5) = 1;
    Predictedhat(PredictedLabels < .5) = 0;
    cvErr = sum(sum(abs(TestLabels - Predictedhat)))./(sum(folds==indexFold));
    
end
    
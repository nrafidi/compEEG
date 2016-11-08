function [SVMmodel, BoxConstraint] = fitAndTuneSVM(data, labels, kernel, numCVFolds)

boxConstraints = [1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 1e1, 1e2, 1e3, 1e4];
nB = length(boxConstraints);
N = size(data, 1);
folds = crossvalind('Kfold', N, numCVFolds);
meanAccPerFold = nan(numCVFolds, nB);
classNames = unique(labels);

for b = 1:nB
    tic
    for f = 1:numCVFolds
        SVMmodel = fitcsvm(data(folds ~= f,:),labels(folds ~= f),'KernelFunction',kernel,...
            'BoxConstraint',boxConstraints(b),'ClassNames',classNames);
        testLabelHat = predict(SVMmodel, data(folds == f,:));
        meanAccPerFold(f, b) = mean(testLabelHat == labels(folds == f));
    end
    toc
end

meanAccPerFold = mean(meanAccPerFold, 1);
[~, bestBox] = max(meanAccPerFold);
BoxConstraint = boxConstraints(bestBox);
SVMmodel = fitcsvm(data, labels, 'KernelFunction', kernel, ...
    'BoxConstraint', BoxConstraint, 'ClassNames', classNames);

end
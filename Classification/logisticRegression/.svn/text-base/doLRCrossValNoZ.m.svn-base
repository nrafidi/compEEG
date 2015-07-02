function errs = doLRCrossValNoZ(X, Y, chooseLambda, regBias, folds, ...
    numFolds, doSave, fname)
% doLRCrossValNoZ runs cross validated logistic regression with L2 penalty
% and returns errors on each fold
% X = data feature set (NxF)
% Y = data labels (Nx?) If Y is a matrix and not a vector, multiple
% logistic regression is used.
% chooseLambda: a 2D vector, if chooseLambda(1) == 1, lambda
% will be found by chooseLambda(2)-fold cross-validation
% regBias: if true, adds a feature of all 1s to X, and bias term will be
% learned with L2 penalty. If not, learns bias separately.
% folds = fold assignments of data
% numFolds = does a numFolds-fold cross-validation
% if doSave, will save errs and the learned weights to fname

if regBias
    X(:,end+1) = 1;
    numFeat = size(X,2);
else
    numFeat = size(X,2)+1;
end
errs = -1*ones(numFolds,1);

if size(Y,2) > 1
    betas = zeros(numFolds, numFeat, size(Y,2));
else
    betas = zeros(numFolds, numFeat);
end
for f = 1:numFolds
    if size(Y,2) > 1
        weights = logRegMult(X(folds~=f, :), Y(folds~=f,:), chooseLambda, regBias);
    else
        weights = logReg(X(folds~=f, :), Y(folds~=f,:), chooseLambda, regBias);
    end
       
    if size(Y,2) > 1
        betas(f, :, :) = weights;
    else
        betas(f, :) = weights;
    end
    Yhat = logPred(X(folds == f,:), weights, regBias);
    errs(f) = sum(abs(Y(folds==f,:) - Yhat))/length(Yhat);
end

if doSave
    save(fname, 'errs', 'betas');
end

end


% Predicts Yhat from X and weights. If regBias is true, we assume that 1s
% have already been added to X. If not, we add the 1s now.
function Yhat = logPred(X, weights, regBias)

if ~regBias
    X(:,end+1) = 1;
end
Yhat= exp(X*weights)./(1 + exp(X*weights));
Yhat(Yhat >= .5) = 1;
Yhat(Yhat < .5) = 0;
end
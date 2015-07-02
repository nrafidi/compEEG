function errs = doLRCrossValWinZ(X, Y, chooseLambda, regBias, folds, ...
    numFolds, doSave, fname)
%doLRCrossValWinZ runs cross validated logistic regression with L2 penalty
% and returns errors on each fold, zscoring X within each cross-validation
% fold
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

errs = -1*ones(numFolds,1);

numFeat = size(X,2)+1;

if size(Y,2) > 1
    betas = zeros(numFolds, numFeat, size(Y,2));
else
    betas = zeros(numFolds, numFeat);
end
for f = 1:numFolds
    %Zscore the training fold
    [Xtr, mu, sig] = zscore(X(folds~=f, :));
    mu = repmat(mu, sum(folds==f), 1);
    sig = repmat(sig, sum(folds==f),1);
    
    if regBias
        Xtr(:,end+1) = 1; %#ok<*AGROW>
    end
    
    if size(Y,2) > 1
        weights = logRegMult(Xtr, Y(folds~=f,:), chooseLambda, regBias);
    else
        weights = logReg(Xtr, Y(folds~=f,:), chooseLambda, regBias);
    end
    
    if size(betas, 2) ~= size(weights,1)
        fprintf('ruh roh\n');
    end
    
    if size(Y,2) > 1
        betas(f, :, :) = weights;
    else
        betas(f, :) = weights;
    end
    % Adjust the test fold
    Xts = (X(folds==f, :) - mu)./sig;
    Xts(:, end+1) = 1;
    Yhat = logPred(Xts, weights);
    errs(f) = sum(abs(Y(folds==f,:) - Yhat))/length(Yhat);
end

if doSave
    save(fname, 'errs', 'betas');
end

end

% Predicts Yhat from X and weights.
function Yhat = logPred(X, weights)
Yhat= exp(X*weights)./(1 + exp(X*weights));
Yhat(Yhat >= .5) = 1;
Yhat(Yhat < .5) = 0;
end
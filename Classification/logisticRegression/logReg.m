%NOTE: this code is a work in progress. Feel free to edit it, but please
%comment what you change and sign with your initials -NSR

% TO DO: different lambdas for each feature? Uncertain
function [weights, lambda] = logReg(X, Y, chooseLambda, regBias, varargin)
% logReg learns a logistic regression between X and Y.
% chooseLambda: if
% if regBias is true, then it will be assumed that a feature of ones has
% already been added to X and the bias will be calculated just like the
% other weights. If not, then the bias will be calculated without the ridge
% penalty, assuming there is no feature of 1s in X.

if chooseLambda(1) == 1
    % Choose lambda with k-fold cross validation (specified by the second
    % entry of chooseLambda
    folds = crossvalind('Kfold', size(X,1), chooseLambda(2));
    lambdas = fliplr([1e-6 1e-5 1e-4 1e-3 1e-2 .1 1 10 100 1000]);
    cverr = zeros(size(lambdas));
    lambFolds = chooseLambda(2);
    if nargin > 4
        weights_fold = varargin{1};
    else
        weights_fold = zeros(size(X, 2) + 1, 1);
    end
    for l = 1:length(lambdas) %was par
        lamb = lambdas(l);
        for k = 1:lambFolds
%             tic
            weights_fold = gradDesc(X(folds~=k,:), Y(folds~=k), lamb, regBias, weights_fold(1:(end-1)));
%             toc
            Y_fold = logPred(X(folds ==k,:), weights_fold, regBias);
            Yhat = zeros(size(Y_fold));
            Yhat(Y_fold >= .5) = 1;
            Yhat(Y_fold < .5) = 0;
            cverr(l) = cverr(l) + sum(abs(Y(folds==k) - Yhat))./(k*sum(folds==k));
        end
    end
    [~, ind] = min(cverr);
    lambda = lambdas(ind);
else
    lambda = chooseLambda(1);
end

weights = gradDesc(X, Y, lambda, regBias);

end

function weights = gradDesc(X, Y, lambda, regBias, varargin)
[N, F] = size(X);
% step = 1e-2;
% step = 1;
eps = (1/sqrt(N))*0.01;
step = eps/10;
tr_prev = 10e6*N;
tr = tr_prev/N;
if nargin > 4
    weights = varargin{1};
else
    weights = zeros(F, 1);
end
times = 0;
if regBias
    while abs(tr - tr_prev) > eps
        times = times + 1;
        pred_err = (Y - logPred(X, weights, regBias));
        pred_err = repmat(pred_err, 1, F);
        weights = weights + step*sum(X.*pred_err, 1)' - step*lambda*weights;
        tr_prev = tr;
        tr = norm(Y - logPred(X, weights, regBias));
        if tr < tr_prev
            step = step*1.1;
        else
            step = step*0.1;
        end
    end
else
    weights(end+1) = 0;
    weights_prev = weights;
    while abs(tr - tr_prev) > eps
        times = times + 1;
        
        % Accelerated code
        v = weights + ((times - 1)/(times +2))*(weights - weights_prev);
        weights_prev = weights;
        step = 1;
        weights_plus = v - step*grad(X, Y, v, regBias, lambda);
        while g(X, Y, weights_plus, regBias, lambda) > g(X, Y, v, regBias, lambda) ...
                + grad(X, Y, v, regBias, lambda)'*(weights_plus - v) + (1/(2*step))*norm(weights_plus - v)^2
            step = 0.5*step; %some beta
            weights_plus = v - step*grad(X, Y, v, regBias, lambda);
        end
        weights = weights_plus;
        % Non-accelerated code
        %                 pred_err = (Y - logPred(X, weights, regBias));
        %                 weights(F+1) = weights(F+1) + step*sum(pred_err);
        %                 pred_err = repmat(pred_err, 1, F);
        %                 weights(1:F) = weights(1:F) + step*sum(X.*pred_err, 1)' - step*lambda*weights(1:F);
        %                 if tr < tr_prev
        %                     step = step*1.1;
        %                 else
        %                     step = step*0.1;
        %                 end
        
        tr_prev = tr;
        tr = norm(Y - logPred(X, weights, regBias));
        
        if times > 1000
            break
        end
    end
    %         disp(tr)
    %         disp(times)
end
end

% Weights is an Fx1 vector, where F is the number of features in X
function Yhat = logPred(X, weights, regBias)

Xhat = X;
if ~regBias
    Xhat(:,end+1) = 1;
end
Yhat= exp(Xhat*weights)./(1 + exp(Xhat*weights));

end

function val = grad(X, Y, weights, regBias, lambda)
% Assumes you're passing in one column of the weight matrix, and one column
% of Y
F = size(X, 2);
val = zeros(F+1, 1);

est = logPred(X, weights, regBias);

if length(est) ~= length(Y)
    keyboard;
end

pred_err = (Y - est);
if ~regBias
    val(F+1) = -sum(pred_err);
end
pred_err_j = repmat(pred_err, 1, F);
val(1:F) = -sum(X.*pred_err_j, 1)' - lambda*weights(1:F);
end

function val = g(X, Y, weights, regBias, lambda)
if ~regBias
    X(:, end+1) = 1;
end
% Assumes Y and weights are vectors
val = sum(-Y.*(X*weights) + log(1 + exp(X*weights))) + lambda*norm(weights)^2;
end
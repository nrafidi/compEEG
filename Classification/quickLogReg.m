X = squeeze(mean(data, 3));

Y = labels(:,1);


trainSet = 13:180;
testSet = 181:191;

B = mnrfit(X(trainSet, :), Y(trainSet) + 1);

Yhat = mnrval(B, X(testSet,:));

Yhat_thresh = double(Yhat(:,2) > 0.5);
err = sum(Yhat_thresh ~= Y(testSet))/length(testSet)
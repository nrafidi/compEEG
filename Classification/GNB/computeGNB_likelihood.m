function [ loglikelihood ] = computeGNB_likelihood( GNBmodel, data, labels )

numSamp = size(data, 1);

likelihood_allclass = nbayes_apply(data, GNBmodel);

likelihood_perclass = nan(numSamp, 1);

for i = 1:numSamp
    likelihood_perclass(i) = likelihood_allclass(i, labels(i));
end

loglikelihood = sum(likelihood_perclass);


end


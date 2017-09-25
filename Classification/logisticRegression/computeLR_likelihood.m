function [ loglikelihood ] = computeLR_likelihood( B, data, labels )

numSamp = size(data, 1);

P = [data ones(numSamp, 1)]*B;
P = exp(P)./(1 + exp(P));

loglikelihood = sum(labels.*log10(P) + (1-labels).*log10(1-P));


end


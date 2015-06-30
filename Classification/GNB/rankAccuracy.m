% rankAcc = rankAccuracy(logprobs, correctLabels, GNBmodel)
%
% Input: logprobs: a mxn array, where the i,j entry gives the log probability of label j for example
%  i correctLabels: a mx1 column vector that contains the correct labels for the m examples (e.g.,
%  [2 1 4 3]') GNBmodel: a GNB model, as trained by nbayes_train.m (used by this code to determine
%  mapping between labels and columns of logprobs)
%
% Output:
%  a mx1 column vector of rank accuracies, from 0 to 1.0, for each of the m predictions
%
% HISTORY:
% 3/10/2014 created by Tom

function rankAcc = rankAccuracy(logprobs, correctLabels, GNBmodel)
    [Y predictedLabelIndexes]=sort(logprobs,2,'descend');
    nLabels=size(logprobs,2);
    rankAcc=zeros(length(correctLabels),1);
    for m=1:length(correctLabels)
        rank=find(correctLabels(m)==GNBmodel.labelVocab(predictedLabelIndexes(m,:)));
        rankAcc(m)=1-((rank-1)/(nLabels-1));
    end
    

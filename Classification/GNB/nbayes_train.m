% nbayesModel =  nbayes_train(exampleMatrix, labels, poolVarEstimate, <classProbabilities)
%
% Train a Gaussian Naive Bayes classifier 
%
% INPUTS:
% examples : a mxn matrix with one training example per row
% labels : a 1xm column vector of integer labels from 1 to L
%
%  alternatively, you can input probabilistic labels by defining this variable to be a structure
%  with two fields:
%     labels.labelVocab: column vector of unique names of classes (e.g., [1 3]')
%     labels.classProbs: matrix where i,j entry gives probability that ith example belongs to jth class
%
% poolVarEstimate : if 1 pool estimate of sigma across all labels, else don't
%
% RETURNS:
% A nbayesModel structure, with fields 
%   labelVocab : a 1xk row vector of unique labels 
%   perLabelModels : a 1xlength(labelVocab) list of structs, each defining the P(Y,X=label)
%      distribution for one label. The ith model in this list is for the ith label in labelVocab.
%      Each of the members of perLabelModels is itself a struct, which looks like this:
%               label: 0
%          labelPrior: 0.5000
%                  mu: [3.5000 39.5000 7]
%               sigma: [0.7071 6.3640 1.4142]
%
% History: created March 2014 by Tom
%  4/14/2014 Tom added ability to handle labels as input structure assigning probabilistic labels

function rslt =  nbayes_train(examples, labels, poolVarEstimate)
    if isstruct(labels)
        rslt.labelVocab=labels.labelVocab;
    else
        rslt.labelVocab=unique(labels);
    end
    perLabelModels=[];
    for i=1:length(rslt.labelVocab) % train model for each possible label
        lab=rslt.labelVocab(i);
        rslt.perLabelModels(i).label=lab;
        if isstruct(labels) % using probabilistic labels
            totalClassProbMass=sum(labels.classProbs(:,i));
            rslt.perLabelModels(i).labelPrior = totalClassProbMass/size(examples,1);
            rslt.perLabelModels(i).mu = (labels.classProbs(:,i)'*examples)./totalClassProbMass;
            meanCtrEx=examples - repmat(rslt.perLabelModels(i).mu,size(examples,1),1);
            rslt.perLabelModels(i).sigma = sqrt((labels.classProbs(:,i)'*(meanCtrEx.^2))./totalClassProbMass);
        else % using simple hard labels
            labIdxs=find(lab == labels);
            rslt.perLabelModels(i).labelPrior = length(labIdxs)/size(examples,1);
            rslt.perLabelModels(i).mu = mean(examples(labIdxs,:), 1);
            rslt.perLabelModels(i).sigma = std(examples(labIdxs,:), [], 1);
        end
    end

    if poolVarEstimate % replace label-specific sigma estimates by their mean
        sigmas=zeros(1,size(examples,2));
        for i=1:length(rslt.labelVocab) % sum label-specific sigmas, div by nlabels
            sigmas=sigmas+rslt.perLabelModels(i).sigma;
        end
        sigmas=sigmas/length(rslt.labelVocab);
        
        for i=1:length(rslt.labelVocab)
            rslt.perLabelModels(i).sigma=sigmas;
        end
    end

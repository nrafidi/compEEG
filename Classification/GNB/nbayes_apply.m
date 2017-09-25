% rslt =  nbayes_apply(examples, model, <selectedFeatures>)
%
% Apply a given Gaussian Naive Bayes classifier to a set of examples
%
% INPUTS:
% examples : a mxn matrix with one training example per row
% model : a NBayes classifier trained using nbayes_train
% <selectedFeatures> : (optional), if provided, a set of indices into the columns of examples, selecting a subset of features that will be used
%
% RETURNS:
%  an m x number-of-labels matrix, where element i,j gives ln [P(y_j) \prod_k P(X_ik | y_j)]
%      Note sorting candidate y_j|X_i within a row gives same as sorting on P(y_j | X_i)
%
% Example:  see EXAMPLE.m for simple script to train a GNBmodel, then apply it using this funciton
%
% History:
% 3/4/2014 created by Tom
% 3/8/2014 added optional argument <selectedFeatures> to allow using just a subset of mus, sigmas


function rslt =  nbayes_apply(examples, model, varargin)
    
    if nargin > 2
        selectedFeats=varargin{1};
    else
        selectedFeats=0;
    end
    
    rslt=zeros(size(examples,1),length(model.labelVocab));
    
    
    for k=1:length(model.perLabelModels) % enter rslt row for ith label
        mu= model.perLabelModels(k).mu;
        sigma= model.perLabelModels(k).sigma;
        labelPrior= model.perLabelModels(k).labelPrior;
        for i=1:size(examples,1)
            rslt(i,k)=loggaussian(examples(i,:),mu,sigma,selectedFeats,1);
        end
    end
% % %     log(labelPrior)+

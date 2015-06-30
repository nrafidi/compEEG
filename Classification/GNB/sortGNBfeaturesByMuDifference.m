% sortedFeats = sortGNBfeaturesByMuDifference(GNBmodel, <useSigmas (default=0)>)
%
% Inputs:
%  GNBmodel = a GNBmodel trained using the function nbayes_train.m
%  useSigmas = an optional argument (default=0).  If 0, ranks features by distance between class-conditional 
%              means.  If 1, it also uses class conditional standard deviations.
%
% Assumes input GNBmodel is for a N-class classification.  Returns the indices of the
% GNBmodel.perLabelModels(1).mu, sorted by the average magnitude of the pairwise difference between their means:
%      GNBmodel.perLabelModels(j).mu - GNBmodel.perLabelModels(k).mu, averaged over all j ~= k.
%
% If optional argument useSigma is 1 (default is 0), then normalize the distance by the standard
% deviations of the models j and k.
%
% Intended to be used as a simple feature selection method.  
%
% History
% 3/8/2014 Created by Tom
% 3/9/2014 Tom generalized to case where GNBmodel has more than two class labels
% 3/20/2014 Tom added optional argument useSigmas.  but tests seem to work best with useSigma=0...


function sortedFeats = sortGNBfeaturesByMuDifference(varargin)
    
    GNBmodel=varargin{1};
    if nargin>1
        useSigmas=varargin{2};
    else
        useSigmas=0; % the default
    end
        
    totalMuDiff=zeros(size(GNBmodel.perLabelModels(1).mu));
    for k=1:length(GNBmodel.labelVocab)-1
        for m=k+1:length(GNBmodel.labelVocab)
            muDiff=abs(GNBmodel.perLabelModels(k).mu-GNBmodel.perLabelModels(m).mu);
            if useSigmas
                % here we get the mean sigma for models k and m, and nomralize distances by these.
                sigmas= 0.5*(GNBmodel.perLabelModels(k).sigma + GNBmodel.perLabelModels(m).sigma);
                muDiff = muDiff ./ sigmas;
                
                % here we normalize the distance between means by the std deviations, averaging the two.
                % muDiff= 0.5*((muDiff./GNBmodel.perLabelModels(k).sigma)+(muDiff./GNBmodel.perLabelModels(m).sigma));
            end
            totalMuDiff=totalMuDiff+muDiff;
        end
    end
    
 [vals sortedFeats]=sort(totalMuDiff, 'descend'); 
    

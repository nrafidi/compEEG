% loggaussian(x,mu,sigma,<seleectedFeatures 0>,<ignoreNormalizer 0>)
%
% Inputs:
%  x: row vector
%  mu: row vector of means
%  signma: row vector of std deviations
%  <seleectedFeatures> = (optional) a subset of the indices of mu, in which case only these will be
%            used to calculate result. If 0 or unsupplied, all of the dimensions of mu will be used.
%  <ignoreNormalizer> = (optional) if 1, then avoid computing the normalization term. This helps
%
% Output:
%  ln of probability of x under a naive Gaussian distribution (with a diagonal covariance matrix)
%
% Example: loggaussian([2 3], [3 3], [1 1])  = -2.3379
%
% History:
% Created 3/2014 by Tom 
% 3/8/2014 added optional argument <selectedFeatures> to allow using just a subset of mus, sigmas
% 11/9/2015 Modified calculation of normalizer to avoid underflow, vectorized and generally cleaned up - Nicole
function logp = loggaussian(x,mu,sigma, varargin)
    
    if isempty(varargin) || isequal(0,varargin{1})
        selectedFeats= 1:length(x);
    else
        selectedFeats=varargin{1};
    end    
    
    if length(varargin)>1
        ignoreNormalizer=varargin{2};
    else
        ignoreNormalizer=0; % default
    end
    
    if ignoreNormalizer
        logp=0;
    else
        logp= length(selectedFeats)*log((1/sqrt(2*pi)));
    end
    if any(sigma == 0)
%         fprintf('Zero variance???\n');
%         keyboard;
        sigma = sigma + eps;
    end
    logp = logp - sum(log(sigma(selectedFeats)) + 0.5*((x(selectedFeats) - mu(selectedFeats))./sigma(selectedFeats)).^2);
    

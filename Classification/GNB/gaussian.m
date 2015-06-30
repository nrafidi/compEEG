% gaussian(x,mu,sigma)
%
% Inputs:
%  x: row vector
%  mu: row vector of means
%  signma: row vector of std deviations
% Output:
%  probability of x under a naive Gaussian distribution (with a diagonal covariance matrix)
%
% Example: gaussian([2 3], [3 3], [1 1])
%
% History:
% Created 3/2014 by Tom

function pr = gaussian(x,mu,sigma)
    normalizer=0.3989;  % this is 1/sqrt(2*pi)
    pr=1;
    for i=1:length(x)
        pr= pr* (normalizer/sigma(i)) * exp(-0.5 * (x(i)-mu(i))^2/sigma(i));
    end
    
    

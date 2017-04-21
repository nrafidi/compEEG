function [risk, odd]=odds(varargin)
%ODDS
% This function calculates the Risk Ratio and the Odds Ratio (OR) on a 2x2
% input matrix. Both ratios are computed with confidence intervals. If
% confidence interval of OR doesn't encompass the value OR=1, then the
% function computes the Bayesian Credibility Assessment of the test. If the
% test is credible, the function calculates the Association Parameter Phi.
% The association parameter Phi=sqrt(chisquare/N).
% The routine coumputes the Power and, if necessary, the sample sizes needed
% to achieve a power=0.80 using a modified asymptotic normal method with
% continuity correction as described by Hardeo Sahai and Anwer Khurshid in
% Statistics in Medicine, 1996, Vol. 15, Issue 1: 1-21.
%
% Syntax: 	ODDS(X,FL)
%      
%     Inputs:
%           X - 2x2 data matrix composed like this: 
%.............................................Cases...Controls
%                                              ___________
%Treated (or exposed to risk factor)          |  A  |  B  |
%                                             |_____|_____|
%Placebo (or not exposed to risk factors )    |  C  |  D  |
%                                             |_____|_____|
%                                                
%           ALPHA - Significance level (default=0.05).
%
%     Outputs:
%           - Risk Ratio with Confidence interval.
%           - Absolute Risk Reduction.
%           - Relative Risk Reduction.
%           - Odds Ratio wirh Confidence interval
%           - Critical Odds Ratio (Bayesian Credibility Assessment)
%           - Phi (association parameter)
%           - Power and sample sizes calculation
%
%      Example: 
% ..............................Cancer..Controls
%                                ___________
% Passive smoke exposed         |  25 |  21 |
%                               |_____|_____|
% Passive smoke not exposed     |  7  |  27 |
%                               |_____|_____|
%
% Data matrix must be x=[25 21; 7 27];
%
% Calling on Matlab the function: odds(x)
% answer is:
%
% Significance level: 95%
%  
% Risk Ratio: 2.1021<2.6398<4.0597
% Absolute risk reduction: 33.8%
% Relative risk reduction: 62.1%
%  
% Odds Ratio: 1.6662<4.5918<12.6544
% Phi: 0.3149
% Moderate positive association (risk factor)
%  
% Bayesian Credibility Assessment
% Critical Odds Ratio: 2.4664
% OR>COR. Test is credible at the 95%
%  
% alpha = 0.0500  n1 = 46  n2 = 34
% Z1-b = 1.1924  Power (2-tails) = 0.8834
%
%           Created by Giuseppe Cardillo
%           giuseppe.cardillo-edta@poste.it
%
% To cite this file, this would be an appropriate format:
% Cardillo G. (2007) Odds: compute odds and risk ratio on a 2x2 matrix. 
% http://www.mathworks.com/matlabcentral/fileexchange/15347


%Input error handling
global x
args=cell(varargin);
nu=numel(args);
if isempty(nu)
    error('Warning: Matrix of data is missed...')
elseif nu>2
    error('Warning: Max two input data are required')
end
default.values = {[],0.05};
default.values(1:nu) = args;
[x alpha] = deal(default.values{:});
%check if x is a 2x2 matrix
if ~isequal(size(x),[2 2])
    error('Input matrix must be a 2x2 matrix')
end
if ~all(isfinite(x(:))) || ~all(isnumeric(x(:)))
    error('Warning: all X values must be numeric and finite')
end
if ~isequal(x(:),round(x(:)))
    error('Warning: X data matrix values must be whole numbers')
end
if ~isempty(x(x==0))
    error('It is impossible to proceed. Zeros are present in the matrix.')
end
if nu==2 %if necessary check alpha
    if ~isscalar(alpha) || ~isnumeric(alpha) || ~isfinite(alpha) || isempty(alpha)
        error('Warning: it is required a numeric, finite and scalar ALPHA value.');
    end
    if alpha <= 0 || alpha >= 1 %check if alpha is between 0 and 1
        error('Warning: ALPHA must be comprised between 0 and 1.')
    end
end
clear args default nu

% fprintf('Significance level: %d%%\n', (1-alpha)*100)
% disp(' ')
Za=-realsqrt(2)*erfcinv(2-alpha);
R=sum(x,2); %sum of the rows
p=x(:,1)./R; 

%risk ratio (RR)
rr=p(1)/p(2); 
% rrse=realsqrt(abs(sum((1-p)./x(:,1))))'; %standard error of log(RR)
rrse=realsqrt(sum(1./x(:,1)-1./R));
rrci=exp(reallog(rr)+([-1 1].*(Za*rrse))); %RR confidence interval
d=abs(diff(p)); %absolute risk reduction
rrr=d/p(1); %relative risk reduction
% fprintf('Risk Ratio: %0.4f<%0.4f<%0.4f\n',rrci(1),rr,rrci(2))
% if(rrci(1)<=1 && rrci(2)>=1)
%    disp('Confidence interval encompasses RR=1. None significative association.')
% end
% fprintf('Absolute risk reduction: %0.1f%%\n',d*100)
% fprintf('Relative risk reduction: %0.1f%%\n',rrr*100)
% disp(' ')

%odd ratio (OR)
or=prod(diag(x))/prod(diag(rot90(x))); 
orse=realsqrt(sum(x(:).^-1)); %standard error of log(OR)
orci=exp(reallog(or)+([-1 1].*(Za*orse))); %OR confidence interval
% fprintf('Odds Ratio: %0.4f<%0.4f<%0.4f\n',orci(1),or,orci(2))
% if(orci(1)<=1 && orci(2)>=1)
%    disp('Confidence interval encompasses OR=1. None significative association.')
% else
%     Phi;
%     disp(' ')
%     disp('Bayesian Credibility Assessment')
%     orci=reallog(orci); 
%     cor=exp(-diff(orci)^2/(4*realsqrt(prod(orci)))); %Critical odds ratio (COR)
%     if or<1
%         fprintf('Critical Odds Ratio: %0.4f\n',cor)
%         if or<cor
%             fprintf('OR<COR. Test is credible at the %d%%\n',(1-alpha)*100)
%         else
%             fprintf('OR>=COR. Test isn''t credible at the %d%%\n',(1-alpha)*100)
%         end
%     else
%         cor=1/cor; %correct cor
%         fprintf('Critical Odds Ratio: %0.4f\n',cor)
%         if or>cor
%             fprintf('OR>COR. Test is credible at the %d%%\n',(1-alpha)*100)
%         else
%             fprintf('OR<=COR. Test isn''t credible at the %d%%\n',(1-alpha)*100)
%         end
%     end
% end
% disp(' ')
% 
% %power (Asymptotic normal method)
% k=R(2)/R(1);
% q=1-p;
% pm=(p(1)+k*p(2))/(k+1);
% qm=1-pm;
% Z1_b=(realsqrt(R(1)*d^2)-Za*realsqrt((1+1/k)*pm*qm))/realsqrt(p(1)*q(1)+p(2)*q(2)/k);
% pwr=0.5*erfc(-Z1_b/realsqrt(2));
% fprintf('alpha = %0.4f  n1 = %d  n2 = %d\n',alpha,R)
% fprintf('Z1-b = %0.4f  Power (2-tails) = %0.4f\n',Z1_b,pwr)
% if pwr<0.8
%     %sample size (Modified Asymptotic normal method with continuity correction)
%     nstar=(Za*realsqrt(pm*qm*(1+1/k))-realsqrt(2)*erfcinv(1.6)*realsqrt(p(1)*q(1)+p(2)*q(2)/k))^2/d^2;
%     n1=round(nstar/4*(1+realsqrt(1+2*(k+1)/(k*d*nstar)))^2);
%     n2=round(k*n1);
%     disp(' ')
%     disp('To achieve a recommended Power=0.80')
%     fprintf('n1 = %d (add %d subjects to exposed row)\n',n1,n1-R(1))
%     fprintf('n2 = %d (add %d subjects to not exposed row)\n',n2,n2-R(2))
% end

risk.mean=rr;
risk.ci=rrci;

odd.mean=or;
odd.ci=orci;
return

function Phi
global x
Phi=(det(x)-sum(x(:))/2)/realsqrt(prod(sum(x))*prod(sum(x,2)));
fprintf('Phi: %0.4f\n',Phi)
switch sign(Phi)
    case -1
        s='negative association (protective factor)';
    case 1
        s='positive association (risk factor)';
end
Phi=abs(Phi);
if Phi<=0.3
    p='Weak ';
elseif (Phi>0.3 && Phi<=0.7)
    p='Moderate ';
else
    p='Strong ';
end
disp([p s])
return
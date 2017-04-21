function forestplot(response, predictor, subgroup, varargin)
% FORESTPLOT generates a forest plot to demonstrate the effects of a 
% predictor in multiple subgroups or across multiple studies. A forest plot
% is a graphical display designed to illustrate the relative strength of 
% treatment effects in multiple quantitative scientific studies addressing 
% the same question.
%
% forestplot(response, predictor, subgroup, subgroup_text, predictor_text, limx, 'ArgumentName',ArgumentValue, ...)
% response  -   a Nx1 binary vector, where N is the number patients, "1"  
%               means an event/disease happened on this patient, "0"
%               means the opposite and "NaN" means unknown. When response 
%               is a NxM matrix, each column is associated with one study.
%               Pad with "NaN" in the end if the patient sizes are 
%               different across studies. 
% preditor  -   a Nx1 binary vector, where "1" means this patient is
%               exposed to this treament or is within desire predictor
%               range (e.g. low dose), "0" mean the opposite and "NaN" means
%               unknown. When response is a NxM matrix, each column is 
%               associated with one study.Pad with "NaN" if the patient 
%               sizes are different across studies.
% subgroup  -   a NxM binary vector, where "1" means this patient belongs
%               to this subgroup "0" mean the opposite and "NaN" means
%               unknown.
%
% Optional 'ArgumentName' includes:
% subgroup_text     -   a 1xM cell vector of strings, where each strings is 
%                       a description of one subgroup/study. For example, 
%                       {'Group 1', 'Group 2', 'Group 3'} or {'Study 1',
%                       'Study 2', 'Study 3'}
% predictor_text    -   a 1x2 cell vector of strings, where the first
%                       string is the description of the non-treatment
%                       group and the second is for the treatment group.
%                       For example, {'Favor high dose', 'Favor low dose'}
% limx  -   a 1x2 vector to specify the limits of x coordinate. The plot is
%           in logarithm scale.
% 'stat'    -   the statistics (odds ratio or relative risk) used to 
%               generate the forest plot. (options = 'or', 'rr', default = 
%               'or')
% For example, 
%   response=round(rand(100, 1)); 
%   preditor=round(rand(100, 1)); 
%   subgroup=round(rand(100, 3));
%   forestplot(response, preditor, subgroup, 'stat', 'rr');
% plots a forest plot.
% 
%   Developed under Matlab version 7.10.0.499 (R2010a)
%   Created by Qi An
%   anqi2000@gmail.com


%   QA 6/1/2012 initial skeleton
%   QA 7/23/2012 updated description; added more error handling
%   QA 2/13/2013 Fixed a couple of bugs; simplified the input structure


%% default values and error handlings
if nargin<3
    error('At least 3 inputs should be specified! See help for more information.');
end
if size(response, 2)>1 % meta-analysis
    if size(response, 2)==size(predictor, 2) && size(response, 2)==size(subgroup, 2)
    else
        error('When conduting meta-analysis, the number of columns should be the same for response, predictor and subgroup!');
    end
end
for i=1:size(subgroup, 2)
    subgroup_text{i}=['Group ' num2str(i)];
end
predictor_text={'favor no treatment', 'favor treatment'};
limx=[];
stat='or';

narg=nargin-3; % exclude first three inputs
if narg==0 % three inputs
elseif rem(narg,2)==1 % even number of inputs
    error('Optional inputs should be specified in pairs! see help for more information');
else
    for i=1:2:(narg) % assign input values
        switch lower(varargin{i})
            case 'subgroup_text'
                subgroup_text=varargin{i+1};
            case 'predictor_text'
                predictor_text=varargin{i+1};
            case 'limx'
                limx=limx;
            case 'stat'
                if strcmp(varargin{i+1}, 'or')==1
                    stat='or';
                elseif strcmp(varargin{i+1}, 'rr')==1
                    stat='rr';
                else
                    error('Unknow statistics!')
                end
            otherwise 
                error('Unknown input!')
        end
    end
end

    
N=size(response, 1); % number of patients
Ns=size(subgroup, 2); % number of subgroups


%% calculate statistics
rr=[]; % relative risk
odd=[]; % odds ratio
for i=1:Ns
    idxSub=find(subgroup(:,i)==1); % indices of subgroup patients
    
    if size(response, 2)==size(predictor, 2) && size(response, 2)==size(subgroup, 2) % meta-analysis
        idx1=find(response(idxSub, i)==1 & predictor(idxSub, i)==1); % within this study, preditor is high and response is yes
        idx2=find(response(idxSub, i)==1 & predictor(idxSub, i)==0); % within this study, preditor is low and response is yes
        idx3=find(response(idxSub, i)==0 & predictor(idxSub, i)==1); % within this study, preditor is high and response is no
        idx4=find(response(idxSub, i)==0 & predictor(idxSub, i)==0); % within this study, preditor is low and response is no
    elseif size(response, 2)==1 && size(predictor, 2)==1 % subgroup analysis
        idx1=find(response(idxSub, 1)==1 & predictor(idxSub, 1)==1); % within this subgroup, preditor is high and response is yes
        idx2=find(response(idxSub, 1)==1 & predictor(idxSub, 1)==0); % within this subgroup, preditor is low and response is yes
        idx3=find(response(idxSub, 1)==0 & predictor(idxSub, 1)==1); % within this subgroup, preditor is high and response is no
        idx4=find(response(idxSub, 1)==0 & predictor(idxSub, 1)==0); % within this subgroup, preditor is low and response is no
    else
        error('size(response, 2) and size(response, 2) are not equal!')
    end

    
    mtx=[length(idx1), length(idx3); length(idx2), length(idx4)]; % confusion matrix
    if sum(mtx(:)==0)>0 % 0 patient in any one of the category
        rr{i}.mean=nan;
        rr{i}.ci=nan(1,2);
        odd{i}.mean=nan;
        odd{i}.ci=nan(1,2);
    else
        [rr{i}, odd{i}]=odds(mtx, 0.05); % calculate relative risk and odds ratio
    end
end

if strcmp(stat, 'or')==1
    val=odd;
elseif strcmp(stat, 'rr')==1
    val=rr;
end


%% plot figures
figure;
hold on;
for i=1:Ns
    semilogx([val{i}.ci(1),val{i}.ci(2)], [i, i] ,'k-', 'LineWidth', 3); % plot the confidence interval
    if val{i}.mean<1 & ~isnan(val{i}.mean)
        plot(val{i}.mean, i, 'ks','markerSize',ceil(1/val{i}.mean)+3,'MarkerEdgeColor','r','MarkerFaceColor','r'); % plot the boxes
    elseif val{i}.mean>=1 & ~isnan(val{i}.mean)
        plot(val{i}.mean, i, 'ks','markerSize',ceil(val{i}.mean)+3,'MarkerEdgeColor','r','MarkerFaceColor','r'); % plot the boxes
    end
end
ylim([0 Ns+1]);
set(gca, 'ytick', 1:Ns+1, 'yticklabel', [subgroup_text {''} ], 'fontsize', 14, 'xscale', 'log');
if ~isempty(limx)
    xlim(limx);
else
    limx=xlim;
end
plot([1;1], [0, Ns], 'k--');
xrange=diff(log10(limx)); % the range of x coordinate
text(10^(-0.1*xrange), Ns+0.75, predictor_text{1}, 'HorizontalAlignment', 'right', 'fontsize', 14);
text(10^(0.1*xrange), Ns+0.75, predictor_text{2}, 'HorizontalAlignment', 'left', 'fontsize', 14);
if strcmp(stat, 'or')==1
    xlabel('Odds Ratio', 'fontsize', 14);
    text(limx(2)*5, Ns+0.75, 'OR [95% CI]', 'fontsize', 14, 'HorizontalAlignment', 'right');
elseif strcmp(stat, 'rr')==1
    xlabel('Relative Risk', 'fontsize', 14);
    text(limx(2)*5, Ns+0.75, 'RR [95% CI]', 'fontsize', 14, 'HorizontalAlignment', 'right');
end
for i=1:Ns
    text(limx(2)*5, i, sprintf('%4.2f [%4.2f %4.2f]', val{i}.mean, val{i}.ci(1), val{i}.ci(2)), 'HorizontalAlignment', 'right', 'fontsize', 14); % display odds ratio
end
xlim([limx(1), limx(2)*5])
hold off;
   

    
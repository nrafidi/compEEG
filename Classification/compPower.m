%Calculate number of subjects needed to replicate result
resDir = '/Users/nrafidi/Documents/MATLAB/compEEG-data/results/';
power = 0.9;
%% Competition Result

load([resDir 'CompEEG_5FCV_win17_permAccs.mat']);

nullMean = mean(mean(permSubAccs, 3), 2);
nullStd = std(reshape(permSubAccs, 28, []), [], 2);
altMean = mean(trueSubAccs, 2);

subZ = (altMean - nullMean)./nullStd;

% Alpha calculation
meanSub = mean(trueSubAccs, 2);
meanPerm = squeeze(mean(permSubAccs, 3));

pvals = sum(meanPerm > repmat(meanSub, 1, 100), 2)/100;
chi_vals = -2.*log(pvals);
group_pval_comp = 1 - chi2cdf(sum(chi_vals),2*length(pvals));

if group_pval_comp == 0
    alpha = 0.001;
end

sampsizepwr('t', [mean(subZ), std(subZ)], 0, power, [], 'Tail', 'left', 'Alpha', 0.001)

%% KR Result

load([resDir 'KR_analysis_output_600-700ms-MeanDiff_singleItem_Direct_cWin17.mat']);

altMean = AUCs_true;

load([resDir 'KR_analysis_output_Perm_singleItem_Direct_cWin17.mat']);

nullDist = cell2mat(AUCs_true);

nullMean = mean(nullDist);
nullStd = std(nullDist);

pval_kr = sum(altMean < nullDist)/100;

if pval_kr == 0
    alpha = 0.001;
end

sampSizeKR = sampsizepwr('t', [nullMean, nullStd], altMean, power, [], 'Tail', 'right', 'Alpha', alpha);
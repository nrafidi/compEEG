load ../../compEEG-data/results/KR_analysis_output_subMeans_650-800ms-R4_meanTime_subCV.mat

plotKRAnalysis;

bad_draws = find(populationAccHist < 0.4);
numBad = length(bad_draws);
class_preds_bad = class_ests_pop(bad_draws);

mean_pred_bad = nan(numBad, 2);
std_pred_bad = nan(numBad, 2);

for i = 1:numBad
    mean_pred_bad(i,1) = mean(class_preds_bad{i}{1});
    mean_pred_bad(i,2) = mean(class_preds_bad{i}{2});
    std_pred_bad(i,1) = std(class_preds_bad{i}{1});
    std_pred_bad(i,2) = std(class_preds_bad{i}{2});
end

good_draws = find(populationAccHist > 0.4);
numGood = length(good_draws);
class_preds_good = class_ests_pop(good_draws);

mean_pred_good = nan(numGood, 2);
std_pred_good = nan(numGood, 2);

for i = 1:numGood
    mean_pred_good(i,1) = mean(class_preds_good{i}{1});
    mean_pred_good(i,2) = mean(class_preds_good{i}{2});
    std_pred_good(i,1) = std(class_preds_good{i}{1});
    std_pred_good(i,2) = std(class_preds_good{i}{2});
end
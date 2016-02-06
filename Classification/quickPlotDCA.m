resultRoot = '/Users/nrafidi/Documents/MATLAB/compEEG-data/results/';
times = -100:20:760;
numT = length(times);

type = 'zscored_windowed';
numComp = 10;

accs_DCA = nan(numT, 5);
accs_PCA = nan(numT, 5);
% accs_Reg = nan(numT, 5);
for fold = 1:5
    load(sprintf('%sKRClassification_PCA-DCA_fold%d_%s_regLambda_numComp%d.mat', resultRoot, fold, type, numComp));
%     if fold > 1
%         accs_Reg(:, fold) = acc_Reg;
%     else
%         accs_Reg(:, fold) = acc_Reg/205;
%     end
%     load(sprintf('%sKRClassification_PCA-DCA_fold%d_lowerEps.mat', resultRoot, fold));
    accs_DCA(:, fold) = acc_DCA;
    accs_PCA(:, fold) = acc_PCA;
    
end
f = figure;
plot(times, mean(accs_DCA, 2));
hold on
plot(times, mean(accs_PCA, 2));
plot(times, mean(accs_Reg, 2));
line([-50, -50], [0.48 0.64], 'Color', 'k');
line([min(times), max(times)], [0.5 0.5], 'Color', 'k', 'LineStyle', '--');
xlim([min(times), max(times)]);
ylim([0.48 0.64]);
xlabel('Time Relative to Stimulus Onset');
ylabel('Classification Accuracy');
title(sprintf('Subsequent memory prediction accuracy\n50ms window averages, Mean(R1R2) - Mean(R3R4)\nZscored, Windowed, %d comp', numComp));
legend({'DCA', 'PCA', 'Reg'});

saveas(f, sprintf('%sfigures/Reg_DCA_PCA_w50_meanConcat_%s_numComp%d.fig', resultRoot, type, numComp));
saveas(f, sprintf('%sfigures/Reg_DCA_PCA_w50_meanConcat_%s_numComp%d.pdf', resultRoot, type, numComp));
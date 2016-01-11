% Process and plot KR analysis results

meanTrueAcc = mean(AUCs_true, 1);
numDraws = length(AUCs_pop);

numTuples = size(AUCs_true, 2);
populationAccHist = nan(numDraws, numTuples);

for iDraw = 1:numDraws    
    populationAccHist(iDraw,:) = mean(AUCs_pop{iDraw}, 1);
end

% 
% titles = {'R1', 'R2', 'R3', 'R4', 'R1, R2', 'R1, R3', 'R1, R4', ...
%     'R2, R3', 'R2, R4', 'R3, R4', 'R1, R2, R3', 'R1, R2, R4', ...
%     'R1, R3, R4', 'R2, R3, R4', 'R1, R2, R3, R4'};

titles = {'Mean R1R2 - Mean R3R4'};%{};%{'350-600ms-R3R4'};%{'Slope From R1-R4 at 350-600ms'};%;%;{'650-800ms - R4'}
for iTuple = 1:numTuples
    if numTuples > 1
    subplot(4, 4, iTuple);
    end
    hist(populationAccHist(:,iTuple));
    percentAbove = sum(populationAccHist(:,iTuple) > 0.5)/length(populationAccHist(:,iTuple));
    title(sprintf('%s\nPercent above chance = %0.0f', titles{iTuple}, 100*percentAbove));
    xlim([0.35, 0.65]);
    hold on;
    line([meanTrueAcc(iTuple) meanTrueAcc(iTuple)], ...
        [0 600]);
    if iTuple == numTuples
        legend({'Population', 'True'}, 'Location', 'eastoutside');
    end
    hold off;
end
% 
% [minAUC, minInd] = min(populationAccHist);
% figure;
% plot(ROC_X_pop{minInd}{1}, ROC_Y_pop{minInd}{1});
% title(num2str(minAUC));
% 
% probabilities = class_ests_pop{minInd}{1};
% labels = true_labs_pop{minInd}{1};
% thresholds = min(probabilities):0.01:max(probabilities);
% N = sum(labels);
% false_pos = nan(length(thresholds), 1);
% true_pos = nan(length(thresholds), 1);
% for t = 1:length(thresholds)
%     est_pos = double(probabilities > thresholds(t));
%     
%     false_pos(t) = sum(est_pos(~logical(labels)));
%     true_pos(t) = sum(est_pos(logical(labels)));
%     
% end
% figure;
% plot(false_pos, true_pos);
% 
% 
% class_accs = nan(numDraws, 1);
% for i = 1:numDraws
%     ests = [class_ests_pop{i}{1}; class_ests_pop{i}{2}];
%     labs = [true_labs_pop{i}{1}; true_labs_pop{i}{2}];
%     
%     class_accs(i) = sum(double(ests > 0.5) == labs)/length(labs);
% end
% figure;
% histogram(class_accs, 20);
% hold on;
% ests = [class_ests_true{1}; class_ests_true{2}];
% labs = [true_labs_true{1}; true_labs_true{2}];
% true_accs = sum(double(ests > 0.5) == labs)/length(labs);
% line([true_accs true_accs], [0 200], 'Color', 'r');


% suptitle(sprintf('Bootstrap results for different tuples of KR data\nOnly Using Correctly Answered Instances'));
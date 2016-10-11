subjects = {'AA', 'BB', 'DD', 'EE', 'F', 'GG', 'HH', 'JJ', ...
    'K', 'M', 'N', 'O', 'R', 'S', 'T', 'U', 'V', 'X', 'Y', 'Z'};
numSub = length(subjects);

meanRespTraj = nan(numSub, 4);
meanRespTrajCorr = nan(numSub, 4);
meanRespTrajInc = nan(numSub, 4);
meanCorr = nan(numSub, 1);
krTrajList = cell(numSub, 1);
krLabelList = cell(numSub, 1);
for s = 1:numSub
    load(sprintf('/Users/nrafidi/Documents/MATLAB/compEEG-data/behavioral/%s/%s_answerTraj.mat', subjects{s}, subjects{s}));
    [~, ~, corrAnswers] = sortKRTraj(subjects{s});
    
    numTraj = size(responseTraj, 1);
    noNaN = false(numTraj, 1);
    for t = 1:numTraj
        if ~any(isnan(responseTraj(t,:)))
            noNaN(t) = true;
        end
    end
    
    krTrajList{s} = responseTraj(noNaN, :);
    krLabelList{s} = corrAnswers(noNaN);
    meanRespTraj(s,:) = nanmean(responseTraj, 1);
%     figure;
%     bar(meanRespTraj(s,:));
%     title(num2str(mean(corrAnswers)));
    meanRespTrajCorr(s,:) = nanmean(responseTraj(corrAnswers,:), 1);
    meanRespTrajInc(s,:) = nanmean(responseTraj(~corrAnswers,:), 1);
    meanCorr(s) = nanmean(corrAnswers, 1);
end
%%
% scatter(1:4, mean(meanRespTraj, 1));
f1 = figure;
totalMean = [mean(meanRespTraj, 1) mean(meanCorr)];
totalStd = [std(meanRespTraj) std(meanCorr)];
bar(1:5, totalMean);
hold on
set(gca, 'XTickLabel',{'R1', 'R2', 'R3', 'R4', 'T'});
errorbar(1:5, totalMean, totalStd, 'r.');
xlabel('Round of Testing');
ylabel('Fraction Correct');
title(sprintf('Average Performance Across Subjects\nin Each Session 2 Round and in Session 3'));
set(gca, 'fontsize', 16);
set(gcf, 'color', 'w');
export_fig(f1, '/Users/nrafidi/Documents/MATLAB/compEEG-data/results/figures/Behavioral/avgPerf.pdf');

f2 = figure;
corrIncMean = [mean(meanRespTrajCorr, 1);mean(meanRespTrajInc, 1)];
bar(1:4, corrIncMean');
hold on
set(gca, 'XTickLabel',{'R1', 'R2', 'R3', 'R4'});
errorbar(0.855:3.855, corrIncMean(1,:), std(meanRespTrajCorr), 'r.');
errorbar(1.145:4.145, corrIncMean(2,:), std(meanRespTrajInc), 'k.');
legend('Subsequently Remembered', 'Subsequently Forgotten');
xlabel('Round of Testing');
ylabel('Fraction Correct');
title(sprintf('Average Performance Across Subjects\nDivided by Session 3 Performance'));
set(gca, 'fontsize', 16);
set(gcf, 'color', 'w');
export_fig(f2, '/Users/nrafidi/Documents/MATLAB/compEEG-data/results/figures/Behavioral/avgPerf_div.pdf');
%%
doDirect = false;
doTuple = true;
numDraws = 1000;
tuplesToPlot = [1:5, 8, 10, 11, 14, 15];
namesOfTuples = {'R1', 'R2', 'R3', 'R4', 'R1, R2', 'R2, R3', 'R3, R4', ...
    'R1, R2, R3', 'R2, R3, R4', 'R1, R2, R3, R4'};
numTuples = length(tuplesToPlot);
[ ROC_X_true, ROC_Y_true, ROC_T_true, AUCs_true, ...
ROC_X_pop, ROC_Y_pop, ROC_T_pop, AUCs_pop, draws_pop] ...
= runSubjBootstrap_KR(...
krTrajList, krLabelList, numDraws, doDirect, doTuple);

AUCs_true = mean(AUCs_true(:,tuplesToPlot), 1);
AUC_pop_hist = cell2mat(cellfun(@(x) mean(x(:,tuplesToPlot), 1), AUCs_pop, 'UniformOutput', false));

%%
f3 = figure;
for t = 1:numTuples
    subplot(2, 5, t);
    histogram(AUC_pop_hist(:,t), 'FaceColor', 'k', 'EdgeColor', 'k');
    hold on
    line([AUCs_true(t), AUCs_true(t)], [0, floor(numDraws/5)], 'Color', 'r');
%     if t == 1
%         legend('Bootstrap AUC', 'True AUC');
%     end
    xlim([0.5, 0.9]);
    ylim([0, floor(numDraws/5)]);
    title(namesOfTuples{t});
    set(gca, 'fontsize', 13);
end
suptitle('AUC for Subsets of Session 2 Data Predicting Session 3 Results');
set(gcf, 'color', 'w');
screenSize = get(0,'Screensize');
set(gcf, 'Position', screenSize); % Maximize figure
export_fig(f3, '/Users/nrafidi/Documents/MATLAB/compEEG-data/results/figures/Behavioral/predPerf_tuple_bootstrap.pdf');
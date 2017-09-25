subjectsR = {'JJJ', 'III', 'KKK', 'BBB', 'GGG', 'HHH', 'AAA', 'CCC', 'DDD', 'EEE', 'FFF', 'MM', ...
    'OO', 'PP', 'QQ', 'RR', 'SS', 'TT', 'WW',...
    'YY'};

subjectsO = {'AA', 'BB', 'DD', 'EE', 'F', 'GG', 'HH', 'JJ', ...
    'K', 'M', 'N', 'O', 'R', 'S', 'T', 'U', 'V', 'X', 'Y', 'Z'};

numSubR = length(subjectsR);
numSubO = length(subjectsO);
numSub = numSubO + numSubR;

behaveDataRootR = '/Users/nrafidi/Documents/MATLAB/compEEG-data-rep/behavioral/';
resDirR = '/Users/nrafidi/Documents/MATLAB/compEEG-data-rep/results/';

behaveDataRootO = '/Users/nrafidi/Documents/MATLAB/compEEG-data/behavioral/';
resDirO = '/Users/nrafidi/Documents/MATLAB/compEEG-data/results/';

meanRespTraj = nan(numSub, 4);
meanRespTrajCorr = nan(numSub, 4);
meanRespTrajInc = nan(numSub, 4);
meanCorr = nan(numSub, 1);
krTrajList = cell(numSub, 1);
krLabelList = cell(numSub, 1);
configurationPerf = zeros(16,2);
for s = 1:numSubR
    load(sprintf('%s%s/%s_answerTraj.mat', behaveDataRootR, subjectsR{s}, subjectsR{s}));
    corrAnswers = sortKRTrajBehave(subjectsR{s}, resDirR);
    
    numTraj = size(responseTraj, 1);
    noNaN = false(numTraj, 1);
    for t = 1:numTraj
        if ~any(isnan(responseTraj(t,:)))
            noNaN(t) = true;
        end
    end
    
    krTrajList{s} = responseTraj(noNaN, :);
    krLabelList{s} = corrAnswers(noNaN);
    
    for it = 1:size(krTrajList{s},1)
        if krTrajList{s}(it,1) == 0
            if krTrajList{s}(it,2) == 0
                if krTrajList{s}(it,3) == 0
                    if krTrajList{s}(it,4) == 0
                        configurationPerf(1,krLabelList{s}(it,1)+1) = configurationPerf(1,krLabelList{s}(it,1)+1)+1;
                    else
                        configurationPerf(2,krLabelList{s}(it,1)+1) = configurationPerf(2,krLabelList{s}(it,1)+1)+1;
                    end
                else
                    if krTrajList{s}(it,4) == 0
                        configurationPerf(3,krLabelList{s}(it,1)+1) = configurationPerf(3,krLabelList{s}(it,1)+1)+1;
                    else
                        configurationPerf(4,krLabelList{s}(it,1)+1) = configurationPerf(4,krLabelList{s}(it,1)+1)+1;
                    end
                end
            else
                if krTrajList{s}(it,3) == 0
                    if krTrajList{s}(it,4) == 0
                        configurationPerf(5,krLabelList{s}(it,1)+1) = configurationPerf(5,krLabelList{s}(it,1)+1)+1;
                    else
                        configurationPerf(6,krLabelList{s}(it,1)+1) = configurationPerf(6,krLabelList{s}(it,1)+1)+1;
                    end
                else
                    if krTrajList{s}(it,4) == 0
                        configurationPerf(7,krLabelList{s}(it,1)+1) = configurationPerf(7,krLabelList{s}(it,1)+1)+1;
                    else
                        configurationPerf(8,krLabelList{s}(it,1)+1) = configurationPerf(8,krLabelList{s}(it,1)+1)+1;
                    end
                end
            end
        else
            if krTrajList{s}(it,2) == 0
                if krTrajList{s}(it,3) == 0
                    if krTrajList{s}(it,4) == 0
                        configurationPerf(9,krLabelList{s}(it,1)+1) = configurationPerf(9,krLabelList{s}(it,1)+1)+1;
                    else
                        configurationPerf(10,krLabelList{s}(it,1)+1) = configurationPerf(10,krLabelList{s}(it,1)+1)+1;
                    end
                else
                    if krTrajList{s}(it,4) == 0
                        configurationPerf(11,krLabelList{s}(it,1)+1) = configurationPerf(11,krLabelList{s}(it,1)+1)+1;
                    else
                        configurationPerf(12,krLabelList{s}(it,1)+1) = configurationPerf(12,krLabelList{s}(it,1)+1)+1;
                    end
                end
            else
                if krTrajList{s}(it,3) == 0
                    if krTrajList{s}(it,4) == 0
                        configurationPerf(13,krLabelList{s}(it,1)+1) = configurationPerf(13,krLabelList{s}(it,1)+1)+1;
                    else
                        configurationPerf(14,krLabelList{s}(it,1)+1) = configurationPerf(14,krLabelList{s}(it,1)+1)+1;
                    end
                else
                    if krTrajList{s}(it,4) == 0
                        configurationPerf(15,krLabelList{s}(it,1)+1) = configurationPerf(15,krLabelList{s}(it,1)+1)+1;
                    else
                        configurationPerf(16,krLabelList{s}(it,1)+1) = configurationPerf(16,krLabelList{s}(it,1)+1)+1;
                    end
                end
            end
        end
    end
    
    meanRespTraj(s,:) = nanmean(responseTraj, 1);
    %     figure;
    %     bar(meanRespTraj(s,:));
    %     title(num2str(mean(corrAnswers)));
    meanRespTrajCorr(s,:) = nanmean(responseTraj(corrAnswers,:), 1);
    meanRespTrajInc(s,:) = nanmean(responseTraj(~corrAnswers,:), 1);
    meanCorr(s) = nanmean(corrAnswers, 1);
end

for s = 1:numSubO
    load(sprintf('%s%s/%s_answerTraj.mat', behaveDataRootO, subjectsO{s}, subjectsO{s}));
    corrAnswers = sortKRTrajBehave(subjectsO{s}, resDirO);
    
    numTraj = size(responseTraj, 1);
    noNaN = false(numTraj, 1);
    for t = 1:numTraj
        if ~any(isnan(responseTraj(t,:)))
            noNaN(t) = true;
        end
    end
    
    krTrajList{s+numSubR} = responseTraj(noNaN, :);
    krLabelList{s+numSubR} = corrAnswers(noNaN);
    
    for it = 1:size(krTrajList{s},1)
        if krTrajList{s}(it,1) == 0
            if krTrajList{s}(it,2) == 0
                if krTrajList{s}(it,3) == 0
                    if krTrajList{s}(it,4) == 0
                        configurationPerf(1,krLabelList{s}(it,1)+1) = configurationPerf(1,krLabelList{s}(it,1)+1)+1;
                    else
                        configurationPerf(2,krLabelList{s}(it,1)+1) = configurationPerf(2,krLabelList{s}(it,1)+1)+1;
                    end
                else
                    if krTrajList{s}(it,4) == 0
                        configurationPerf(3,krLabelList{s}(it,1)+1) = configurationPerf(3,krLabelList{s}(it,1)+1)+1;
                    else
                        configurationPerf(4,krLabelList{s}(it,1)+1) = configurationPerf(4,krLabelList{s}(it,1)+1)+1;
                    end
                end
            else
                if krTrajList{s}(it,3) == 0
                    if krTrajList{s}(it,4) == 0
                        configurationPerf(5,krLabelList{s}(it,1)+1) = configurationPerf(5,krLabelList{s}(it,1)+1)+1;
                    else
                        configurationPerf(6,krLabelList{s}(it,1)+1) = configurationPerf(6,krLabelList{s}(it,1)+1)+1;
                    end
                else
                    if krTrajList{s}(it,4) == 0
                        configurationPerf(7,krLabelList{s}(it,1)+1) = configurationPerf(7,krLabelList{s}(it,1)+1)+1;
                    else
                        configurationPerf(8,krLabelList{s}(it,1)+1) = configurationPerf(8,krLabelList{s}(it,1)+1)+1;
                    end
                end
            end
        else
            if krTrajList{s}(it,2) == 0
                if krTrajList{s}(it,3) == 0
                    if krTrajList{s}(it,4) == 0
                        configurationPerf(9,krLabelList{s}(it,1)+1) = configurationPerf(9,krLabelList{s}(it,1)+1)+1;
                    else
                        configurationPerf(10,krLabelList{s}(it,1)+1) = configurationPerf(10,krLabelList{s}(it,1)+1)+1;
                    end
                else
                    if krTrajList{s}(it,4) == 0
                        configurationPerf(11,krLabelList{s}(it,1)+1) = configurationPerf(11,krLabelList{s}(it,1)+1)+1;
                    else
                        configurationPerf(12,krLabelList{s}(it,1)+1) = configurationPerf(12,krLabelList{s}(it,1)+1)+1;
                    end
                end
            else
                if krTrajList{s}(it,3) == 0
                    if krTrajList{s}(it,4) == 0
                        configurationPerf(13,krLabelList{s}(it,1)+1) = configurationPerf(13,krLabelList{s}(it,1)+1)+1;
                    else
                        configurationPerf(14,krLabelList{s}(it,1)+1) = configurationPerf(14,krLabelList{s}(it,1)+1)+1;
                    end
                else
                    if krTrajList{s}(it,4) == 0
                        configurationPerf(15,krLabelList{s}(it,1)+1) = configurationPerf(15,krLabelList{s}(it,1)+1)+1;
                    else
                        configurationPerf(16,krLabelList{s}(it,1)+1) = configurationPerf(16,krLabelList{s}(it,1)+1)+1;
                    end
                end
            end
        end
    end
    
    
    meanRespTraj(s+numSubR,:) = nanmean(responseTraj, 1);
    %     figure;
    %     bar(meanRespTraj(s,:));
    %     title(num2str(mean(corrAnswers)));
    meanRespTrajCorr(s+numSubR,:) = nanmean(responseTraj(corrAnswers,:), 1);
    meanRespTrajInc(s+numSubR,:) = nanmean(responseTraj(~corrAnswers,:), 1);
    meanCorr(s+numSubR) = nanmean(corrAnswers, 1);
end

%% Configuration plots

% Total prevalence of each configuration
f = figure;
bar(sum(configurationPerf, 2))
xlim([0.5, 16.5])
set(gca, 'xtick', 1:16);
set(gca, 'xticklabels', {'0000', '0001', '0010', '0011', '0100', '0101', ...
    '0110', '0111', '1000', '1001', '1010', '1011', '1100', '1101', '1110', '1111'});
title('Prevalence of Response Configurations During Session 2');
xlabel('Response Configuration from Session 2')
ylabel('Number of Items')
set(gca, 'FontSize', 16);
set(f, 'Color', 'w');
set (f, 'Units', 'normalized', 'Position', [0,0,1,1]);
export_fig(f, [resDirR 'figures/Behavioral/configPrev.pdf']);
f = figure;
bar(configurationPerf)
xlim([0, 17])
set(gca, 'xtick', 1:16);
set(gca, 'xticklabels', {'0000', '0001', '0010', '0011', '0100', '0101', ...
    '0110', '0111', '1000', '1001', '1010', '1011', '1100', '1101', '1110', '1111'});
legend({'Forgotten', 'Remembered'});
title(sprintf('Prevalence of Response Configurations During Session 2\nSplit by Session 3 Performance'));
xlabel('Response Configuration from Session 2')
ylabel('Number of Items')
set(gca, 'FontSize', 16);
set(f, 'Color', 'w');
set (f, 'Units', 'normalized', 'Position', [0,0,1,1]);
export_fig(f, [resDirR 'figures/Behavioral/configPrevSplit.pdf']);
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
set(gca, 'fontsize', 18);
set(gcf, 'color', 'w');
export_fig(f1, [resDirR 'figures/Behavioral/avgPerf_pooled.pdf']);

f2 = figure;
corrIncMean = [mean(meanRespTrajCorr, 1);mean(meanRespTrajInc, 1)];
bar(1:4, corrIncMean');
hold on
set(gca, 'XTickLabel',{'R1', 'R2', 'R3', 'R4'});
errorbar(0.855:3.855, corrIncMean(1,:), std(meanRespTrajCorr), 'r.');
errorbar(1.145:4.145, corrIncMean(2,:), std(meanRespTrajInc), 'k.');
legend('Subsequently Remembered', 'Subsequently Forgotten', 'Location', 'NorthWest');
xlabel('Round of Testing');
ylabel('Fraction Correct');
title(sprintf('Average Performance Across Subjects\nGrouped by Session 3 Performance'));
set(gca, 'fontsize', 18);
set(gcf, 'color', 'w');
export_fig(f2, [resDirR 'figures/Behavioral/avgPerf_div_pooled.pdf']);
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
    set(gca, 'fontsize', 16);
end
suptitle('AUC for Subsets of Session 2 Data Predicting Session 3 Results');
set(gcf, 'color', 'w');
screenSize = get(0,'Screensize');
set(gcf, 'Position', screenSize); % Maximize figure
export_fig(f3, [resDirR 'figures/Behavioral/predPerf_tuple_bootstrap_pooled.pdf']);
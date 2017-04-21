% Full KR Analysis pipeline

% Starts with Competition and KR EEG data and produces KR population
% accuracy and histogram
subjects = {'III', 'JJJ', 'KKK', 'BBB', 'GGG', 'HHH', 'AAA', 'CCC', 'DDD', 'EEE', 'FFF', 'MM', ...
    'OO', 'PP', 'QQ', 'RR', 'SS', 'TT', 'WW',...
    'YY'};

%VV and UU have files missing

% {'AA', 'BB', 'DD', 'EE', 'F', 'GG', 'HH', 'JJ', ...
%     'K', 'M', 'N', 'O', 'R', 'S', 'T', 'U', 'V', 'X', 'Y', 'Z'};
numSub = length(subjects);
behaveDataRoot = '/Users/nrafidi/Documents/MATLAB/compEEG-data-rep/behavioral/';

krTrajList = cell(numSub, 1);
krRespList = cell(numSub, 1);
krLabelList = cell(numSub, 1);
withinSubAUC = nan(numSub, 1);
accuracyList = [];
numDraws = 1000;
doAUC = true;
doSubMeans = false;
savePlots = true;

featureSetNames = {'660ms-MeanDiff', '640ms-MeanDiff', '680ms-MeanDiff', '600-700ms-MeanDiff','400-500ms-MeanDiff', '350-400ms-Slope', '350-400ms-MeanDiff', '650-800ms-R4', '650-800ms-R3R4', ...
    '350-400ms-R3R4'};
featureSetTimeInds = {39, 38, 40, 36:41, 26:31, 23:26, 23:26, 38:46, 38:46, 23:26};
featureSetRInds = {1:4, 1:4, 1:4, 1:4, 1:4, 1:4, 1:4, 4, 3:4, 3:4};
trueInd = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1];
for compWinToUse = 18%[17, 24, 26, 38]
    fprintf('Comp Win %d\n', compWinToUse);
    for featureSet = 1
        fprintf('%s\n', featureSetNames{featureSet});
        
        for s = 1:numSub
            sub = subjects{s};
            
            fname = ['/Users/nrafidi/Documents/MATLAB/compEEG-data-rep/results/' ...
                sub '/KRanalysis_SlidingFeat_CWin' num2str(compWinToUse) '_Vis.mat'];
            
            load(fname);
            
            if doSubMeans
                krRespList{s} = responseTraj(:, featureSetRInds{featureSet});
                
                krTrajToAddCorr = squeeze(mean(mean(krTraj(logical(krLabels), featureSetRInds{featureSet}, featureSetTimeInds{featureSet}), 3), 1));
                krTrajToAddInc = squeeze(mean(mean(krTraj(~logical(krLabels), featureSetRInds{featureSet}, featureSetTimeInds{featureSet}), 3), 1));
                if strcmp(featureSetNames{featureSet}, '350-400ms-MeanDiff')
                    krTrajToAddCorr = squeeze(mean(krTrajToAddCorr(:,1:2), 2) - mean(krTrajToAddCorr(:,3:4), 2));
                    krTrajToAddInc = squeeze(mean(krTrajToAddInc(:,1:2), 2) - mean(krTrajToAddInc(:,3:4), 2));
                elseif strcmp(featureSetNames{featureSet}, '350-400ms-Slope')
                    krTrajToAddCorr = polyfit(1:4, krTrajToAddCorr, 1);
                    krTrajToAddCorr = krTrajToAddCorr(1);
                    krTrajToAddInc = polyfit(1:4, krTrajToAddInc, 1);
                    krTrajToAddInc = krTrajToAddInc(1);
                elseif strcmp(featureSetNames{featureSet}, '650-800ms-R3R4') || ...
                        strcmp(featureSetNames{featureSet}, '350-600ms-R3R4')
                    krTrajToAddCorr = krTrajToAddCorr(1) - krTrajToAddCorr(2);
                    krTrajToAddInc = krTrajToAddInc(1) - krTrajToAddInc(2);
                end
                
                krTrajList{s} = [krTrajToAddCorr; krTrajToAddInc];
                krLabelList{s} = [1;0];
            else
                krRespList{s} = responseTraj(:, featureSetRInds{featureSet});
                
                krTrajToAdd = squeeze(mean(krTraj(:, featureSetRInds{featureSet}, featureSetTimeInds{featureSet}), 3));
                if strcmp(featureSetNames{featureSet}, '350-400ms-MeanDiff') || ...
                        strcmp(featureSetNames{featureSet}, '400-500ms-MeanDiff') || ...
                        strcmp(featureSetNames{featureSet}, '600-700ms-MeanDiff') || ...
                        strcmp(featureSetNames{featureSet}, '660ms-MeanDiff') || ...
                        strcmp(featureSetNames{featureSet}, '640ms-MeanDiff') || ...
                        strcmp(featureSetNames{featureSet}, '680ms-MeanDiff')
                    krTrajToAdd = squeeze(mean(krTrajToAdd(:,1:2), 2) - mean(krTrajToAdd(:,3:4), 2));
                elseif strcmp(featureSetNames{featureSet}, '350-400ms-Slope')
                    old_krTrajToAdd = krTrajToAdd;
                    numTraj = size(krTrajToAdd, 1);
                    krTrajToAdd = nan(numTraj, 1);
                    for iTraj = 1:numTraj
                        slopeToAdd = polyfit(1:4, ...
                            old_krTrajToAdd(iTraj,:), 1);
                        krTrajToAdd(iTraj) = slopeToAdd(1);
                    end
                elseif strcmp(featureSetNames{featureSet}, '650-800ms-R3R4') || strcmp(featureSetNames{featureSet}, '350-400ms-R3R4')
                    krTrajToAdd = squeeze(mean(krTrajToAdd, 2));
                end
                percCorr = sum(krLabels)/length(krLabels);
                if (percCorr < 0.6) && (percCorr > 0.4)
                    [~, ~, ~, withinSubAUC(s)] = perfcurve(krLabels, krTrajToAdd, 0);
                end
                krTrajList{s} = krTrajToAdd;
                krLabelList{s} = krLabels;
            end
        end
        
        if doAUC
            [ ROC_X_true, ROC_Y_true, ROC_T_true, AUCs_true, ...
                ROC_X_pop, ROC_Y_pop, ROC_T_pop, AUCs_pop, draws_pop] = ...
                runSubjBootstrap_KR(...
                krTrajList, krLabelList, numDraws, true, false, trueInd(featureSet));
            AUCs_true = mean(AUCs_true);
            if AUCs_true < 0.5
                fprintf('Flipped\n');
                AUCs_true = 1 - AUCs_true;
                AUCs_pop = cellfun(@(x) 1 - x, AUCs_pop, 'UniformOutput', false);
            end
            AUCs_pop = cell2mat(AUCs_pop);
            f = figure;
            hist(AUCs_pop)
            hold on
            line([AUCs_true, AUCs_true], [0, 300]);
            xlim([0.4, 0.7]);
            hold off
            xlabel('AUC');
            ylabel('Number of draws');
            legend({'Population', 'True'});
            set(gcf, 'Color', 'w');
            set(gca, 'FontSize', 16);
%             title(sprintf('Single Item AUC, Train on Comp Window %.0f\n%s\nPercent above chance = %.2f', winTime(compWinToUse), featureSetNames{featureSet}, (sum(AUCs_pop > 0.5)/numDraws)*100));
            title(sprintf('AUC Histogram of Subject-level Bootstrap\n660ms Post Stimulus Onset\nPercent above chance = %.1f',  (sum(AUCs_pop > 0.5)/numDraws)*100));
            if savePlots
                export_fig(f, sprintf('/Users/nrafidi/Documents/MATLAB/compEEG-data-rep/results/figures/singleItem_AUC_%s_cWin%d_Vis.fig', featureSetNames{featureSet}, compWinToUse));
                export_fig(f, sprintf('/Users/nrafidi/Documents/MATLAB/compEEG-data-rep/results/figures/singleItem_AUC_%s_cWin%d_Vis.pdf', featureSetNames{featureSet}, compWinToUse));
                save(['../../compEEG-data-rep/results/KR_analysis_output_' featureSetNames{featureSet} '_singleItem_Direct_cWin' num2str(compWinToUse) '_Vis.mat'], ...
                    'ROC_X_true', 'ROC_Y_true', 'ROC_T_true', 'AUCs_true', ...
                    'ROC_X_pop', 'ROC_Y_pop', 'ROC_T_pop', 'AUCs_pop', 'draws_pop');
            end
        else
            [ p_true, p_pop, draws_pop, mean_diff_true, mean_diff_pop ] = runSubjBootstrap_pairedT_KR(...
                krTrajList, krLabelList, numDraws);
            
            %         save(['../../compEEG-data-rep/results/KR_analysis_output_pairedT_' featureSetNames{featureSet} '_meanTime.mat'], ...
            %             'p_true', 'p_pop', 'draws_pop', 'mean_diff_true', 'mean_diff_pop');
        end
    end
end
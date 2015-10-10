% Full KR Analysis pipeline

% Starts with Competition and KR EEG data and produces KR population
% accuracy and histogram
% Examine subject M
subjects = {'AA', 'BB', 'DD', 'EE', 'F', 'GG', 'HH', 'JJ', ...
    'K', 'M', 'N', 'O', 'R', 'S', 'T', 'U', 'V', 'X', 'Y', 'Z'};
numSub = length(subjects);
doSubMeans = 0;
behaveDataRoot = '/Users/nrafidi/Documents/MATLAB/compEEG-data/behavioral/';

krTrajList = cell(numSub, 1);
krRespList = cell(numSub, 1);
krLabelList = cell(numSub, 1);
withinSubAUC = nan(numSub, 1);
accuracyList = [];

numDraws = 1500;
doAUC = false;

featureSetNames = {'Slope', 'MeanDiff', '650-800ms-R4', '650-800ms-R3R4', ...
    '350-600ms-R3R4', 'SlopeCorrOnly'};
featureSetTimeInds = {18:26, 18:26, 33:40, 33:40, 18:26, 18:26};
featureSetRInds = {1:4, 1:4, 4, 3:4, 3:4, 1:4};

for featureSet = 6%:-1:2
    fprintf('%s\n', featureSetNames{featureSet});
    
    for s = 1:numSub
        sub = subjects{s};
        disp(sub);
        
        fname = ['/Users/nrafidi/Documents/MATLAB/compEEG-data/results/' ...
            sub '/KRanalysis_SlidingFeat_lateComp.mat'];
        
        if ~exist(fname, 'file')
            error('This file should exist');
            try %#ok<*UNRCH>
                %% Predict KR
                itemTraj = runKRPrediction(sub);
                
                %% Combine with final quiz answers
                [itemTrajCorr, itemTrajInc, corrAnswers] = sortKRTraj(sub);
                title(sub);
                
                %% Sort for bootstrap analysis
                
                krTraj = [];
                krLabels = [];
                skippedItems = [];
                indCorr = 1;
                indInc = 1;
                for i = 1:length(corrAnswers)
                    if corrAnswers(i)
                        if ~any(isnan(itemTrajCorr(indCorr,:)))
                            krTraj = cat(1, krTraj, itemTrajCorr(indCorr,:));
                            krLabels = cat(1, krLabels, 1);
                        else
                            skippedItems = cat(1, skippedItems, i);
                        end
                        indCorr = indCorr + 1;
                    else
                        if ~any(isnan(itemTrajInc(indInc,:)))
                            krTraj = cat(1, krTraj, itemTrajInc(indInc,:));
                            krLabels = cat(1, krLabels, 0);
                        else
                            skippedItems = cat(1, skippedItems, i);
                        end
                        indInc = indInc + 1;
                    end
                end
                
                save(fname, 'krTraj', 'krLabels', 'skippedItems');
            catch
                fprintf('Subject %s Failed\n', sub);
            end
        else
            load(fname);
        end
        %
        %     load([behaveDataRoot '/' sub '/' sub '_answerTraj.mat']);
        %
        %     numPotTraj = size(responseTraj, 1);
        %     newResponseTraj = [];
        %     responseTraj(isnan(responseTraj)) = 0; %#ok<*SAGROW>
        %     %~any(isnan(responseTraj(iPTraj,:))) &&
        %     for iPTraj = 1:numPotTraj
        %         if ~any(iPTraj == skippedItems);
        %             newResponseTraj = cat(1, newResponseTraj, responseTraj(iPTraj,:));
        %         end
        %     end
        %
        
        %     krTraj = newResponseTraj;%cat(2, krTraj, newResponseTraj);
        
        %     fprintf('Proportion Correct for subject %s = %d\n', sub, ...
        %         sum(krLabels)/length(krLabels));
        
        if doSubMeans
            krRespList{s} = responseTraj(:, featureSetRInds{featureSet});
            
            krTrajToAddCorr = squeeze(mean(mean(krTraj{1}(logical(krLabels{1}), featureSetRInds{featureSet}, featureSetTimeInds{featureSet}), 3), 1));
            krTrajToAddInc = squeeze(mean(mean(krTraj{1}(~logical(krLabels{1}), featureSetRInds{featureSet}, featureSetTimeInds{featureSet}), 3), 1));
            if strcmp(featureSetNames{featureSet}, 'MeanDiff')
                krTrajToAddCorr = squeeze(mean(krTrajToAddCorr(:,1:2), 2) - mean(krTrajToAddCorr(:,3:4), 2));
                krTrajToAddInc = squeeze(mean(krTrajToAddInc(:,1:2), 2) - mean(krTrajToAddInc(:,3:4), 2));
            elseif strcmp(featureSetNames{featureSet}, 'Slope')
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
            
            krTrajToAdd = squeeze(mean(krTraj{1}(:, featureSetRInds{featureSet}, featureSetTimeInds{featureSet}), 3));
            krLabels = krLabels{1};
            if strcmp(featureSetNames{featureSet}, 'MeanDiff')
                krTrajToAdd = squeeze(mean(krTrajToAdd(:,1:2), 2) - mean(krTrajToAdd(:,3:4), 2));
            elseif strcmp(featureSetNames{featureSet}, 'Slope')
                old_krTrajToAdd = krTrajToAdd;
                numTraj = size(krTrajToAdd, 1);
                krTrajToAdd = nan(numTraj, 1);
                for iTraj = 1:numTraj
                    slopeToAdd = polyfit(1:4, ...
                        old_krTrajToAdd(iTraj,:), 1);
                    krTrajToAdd(iTraj) = slopeToAdd(1);
                end
            elseif strcmp(featureSetNames{featureSet}, 'SlopeCorrOnly')
                old_krTrajToAdd = krTrajToAdd;
                numTraj = size(krTrajToAdd, 1);
                krTrajToAdd = [];
                old_krLabels = krLabels;
                krLabels = [];
                for iTraj = 1:numTraj
                    if all(responseTraj(iTraj, 3:end), 2)
                        corrAns = find(responseTraj);
                        slopeToAdd = polyfit(corrAns(1):4, ...
                            old_krTrajToAdd(iTraj,corrAns(1):4), 1);
                        krTrajToAdd = cat(1, krTrajToAdd, slopeToAdd(1));
                        krLabels = cat(1, krLabels, old_krLabels(iTraj));
                    end
                end
            end
            percCorr = sum(krLabels)/length(krLabels);
            if (percCorr < 0.6) && (percCorr > 0.4)
            [~, ~, ~, withinSubAUC(s)] = perfcurve(krLabels, krTrajToAdd, 0);
            end
            krTrajList{s} = krTrajToAdd;
            krLabelList{s} = krLabels;
        end
        %         if size(responseTraj, 1) ~= size(krTrajToAdd, 1)
        %             keyboard;
        %         end
    end
    %         keyboard;
    
    if doAUC
        [ ROC_X_true, ROC_Y_true, ROC_T_true, AUCs_true, ...
            ROC_X_pop, ROC_Y_pop, ROC_T_pop, AUCs_pop, draws_pop] = ...
            runSubjBootstrap_KR(...
            krTrajList, krLabelList, numDraws);
        
        save(['../../compEEG-data/results/KR_analysis_output_' featureSetNames{featureSet} '_meanTime_Direct.mat'], ...
            'ROC_X_true', 'ROC_Y_true', 'ROC_T_true', 'AUCs_true', ...
            'ROC_X_pop', 'ROC_Y_pop', 'ROC_T_pop', 'AUCs_pop', 'draws_pop');
    else
        [ p_true, p_pop, draws_pop, mean_diff_true, mean_diff_pop ] = runSubjBootstrap_pairedT_KR(...
            krTrajList, krLabelList, numDraws);
        
        save(['../../compEEG-data/results/KR_analysis_output_pairedT_' featureSetNames{featureSet} '_meanTime.mat'], ...
            'p_true', 'p_pop', 'draws_pop', 'mean_diff_true', 'mean_diff_pop');
    end
end
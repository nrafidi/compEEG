% Full KR Analysis pipeline

% Starts with Competition and KR EEG data and produces KR population
% accuracy and histogram
% Examine subject M
subjects = {'AA', 'BB', 'DD', 'EE', 'F', 'GG', 'HH', 'JJ', ...
    'K', 'M', 'N', 'O', 'R', 'S', 'T', 'U', 'V', 'X', 'Y', 'Z'};
numSub = length(subjects);

behaveDataRoot = '/Users/nrafidi/Documents/MATLAB/compEEG-data/behavioral/';

krTrajList = cell(numSub, 1);
krRespList = cell(numSub, 1);
krLabelList = cell(numSub, 1);
accuracyList = [];

numDraws = 100;

featureSetNames = {'350-500ms-All', 'MeanDiff', '650-800ms-R4'};
featureSetTimeInds = {18:26, 18:26, 33:40};
featureSetRInds = {1:4, 1:4, 4};

for featureSet = 1:3
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
        krRespList{s} = responseTraj(:, featureSetRInds{featureSet});
        
        krTrajToAdd = squeeze(mean(krTraj{1}(:, featureSetRInds{featureSet}, featureSetTimeInds{featureSet}), 3));
        
        if ~strcmp(featureSetNames{featureSet}, 'MeanDiff')
%             krTrajToAdd = reshape(krTrajToAdd, size(krTrajToAdd, 1), []);
        else
            krTrajToAdd = squeeze(mean(krTrajToAdd(:,1:2), 2) - mean(krTrajToAdd(:,3:4), 2));
        end
        
        krTrajList{s} = krTrajToAdd;
        krLabelList{s} = krLabels{1};
        
        if size(responseTraj, 1) ~= size(krTrajToAdd, 1)
            keyboard;
        end
    end
    
    [ ROC_X_true, ROC_Y_true, ROC_T_true, AUCs_true, ...
        ROC_X_pop, ROC_Y_pop, ROC_T_pop, AUCs_pop ] = runSubjBootstrap_KR(...
        krTrajList, krLabelList, krRespList, numDraws);
    
    save(['../../compEEG-data/results/KR_analysis_output_onlyCorrAns_' featureSetNames{featureSet} '_meanTime.mat'], ...
        'ROC_X_true', 'ROC_Y_true', 'ROC_T_true', 'AUCs_true', ...
        'ROC_X_pop', 'ROC_Y_pop', 'ROC_T_pop', 'AUCs_pop');
    
end
% Analyize and summarize behavioral data

%
subjects = { 'AA', 'BB', 'DD', 'EE', 'F', 'GG', 'HH', 'JJ', ...
    'K', 'M', 'N', 'O', 'R', 'S', 'T', 'U', 'V', 'X', 'Y', 'Z'};
numSub = length(subjects);
item = 1;

behaveDataRoot = '/Users/nrafidi/Documents/MATLAB/compEEG-data/behavioral/';
figureRoot = '/Users/nrafidi/Documents/MATLAB/compEEG-data/results/figures/singleTrial/';
stimRoot = '/Users/nrafidi/Documents/MATLAB/compEEG-stim/';
fid = fopen([stimRoot 'KR_test.txt']);
swahili = textscan(fid, '%s');
swahili = swahili{1};

potentialResponses = [0 0 0 1; 0 0 1 1;  0 1 1 1; 1 1 1 1];
respNames = {'R4', 'R3R4', 'R2R3R4', 'R1R2R3R4'};

% for r = 1:length(respNames)
%     if ~exist([figureRoot respNames{r} '_corr'], 'dir')
%         mkdir([figureRoot respNames{r} '_corr']);
%     end
%     if ~exist([figureRoot respNames{r} '_inc'], 'dir')
%         mkdir([figureRoot respNames{r} '_inc']);
%     end
% end

numPotResp = size(potentialResponses, 1);
krTrajByResponse = cell(numPotResp, 2);
for s = 1:numSub
    sub = subjects{s};
    if ~exist([figureRoot '/' sub '/'], 'dir')
        mkdir([figureRoot '/' sub '/']);
    end
    load([behaveDataRoot '/' sub '/' sub '_answerTraj.mat']);
    load(['/Users/nrafidi/Documents/MATLAB/compEEG-data/results/' ...
        sub '/KRanalysis_SlidingFeat_lateComp.mat']);
    
    krTraj = squeeze(mean(krTraj{1}(:, :, 18:26), 3));
    krLabels = krLabels{1};
    percCorr = floor(100*sum(krLabels)/length(krLabels));
    krLabels = logical(krLabels);
%     responseTraj(skippedItems,:) = []; %#ok<*SAGROW>
    subSwahili = swahili;
    subSwahili(skippedItems) = [];
    
    corrItems =find(krLabels);
    incItems = find(~krLabels);
    
    for cItem = 1:length(corrItems)
        trajID = responseTraj(corrItems(cItem,:), :);
        rowID = ismember(potentialResponses, trajID, 'rows');
        if any(rowID)
            f = figure;
            plot(1:4, krTraj(corrItems(cItem),:), 'b');
            hold on;
            for i = 1:4
                if trajID(i)
                    scatter(i, krTraj(corrItems(cItem),i), 'MarkerEdgeColor', 'blue', ...
                        'Marker', 'o', 'SizeData', 100);
                else
                    scatter(i, krTraj(corrItems(cItem),i), 'MarkerEdgeColor', 'blue', ...
                        'Marker', 'x', 'SizeData', 100);
                end
            end
            xlim([0.9, 4.1]);
            ylim([-0.1 1.1]);
            title(sprintf('%s %s\n%d', sub, subSwahili{corrItems(cItem)}, percCorr));
            hold off;
            saveas(f, [figureRoot '/' sub '/' sub '_' respNames{rowID} '_corr_' subSwahili{corrItems(cItem)} '.png']);
            close(f);
        end
    end
    
    for iItem = 1:length(incItems)
        trajID = responseTraj(incItems(iItem,:), :);
        rowID = ismember(potentialResponses, trajID, 'rows');
        if any(rowID)
            f = figure;
            plot(1:4, krTraj(incItems(iItem),:), 'r');
            hold on;
            for i = 1:4
                if trajID(i)
                    scatter(i, krTraj(incItems(iItem),i), 'MarkerEdgeColor', 'red', ...
                        'Marker', 'o', 'SizeData', 100);
                else
                    scatter(i, krTraj(incItems(iItem),i), 'MarkerEdgeColor', 'red', ...
                        'Marker', 'x', 'SizeData', 100);
                end
            end
            xlim([0.9, 4.1]);
            ylim([-0.1 1.1]);
            title(sprintf('%s %s\n%d', sub, subSwahili{incItems(iItem)}, percCorr));
            hold off;
            saveas(f, [figureRoot '/' sub '/' sub '_' respNames{rowID} '_cinc_' subSwahili{incItems(iItem)} '.png']);
            close(f);
        end
    end
    
end

% Analyize and summarize behavioral data

%
subjects = {'AA'};%, 'BB', 'DD', 'EE', 'F', 'GG', 'HH', 'JJ', ...
%     'K', 'M', 'N', 'O', 'R', 'S', 'T', 'U', 'V', 'X', 'Y', 'Z'};
numSub = length(subjects);
item = 1;

behaveDataRoot = '/Users/nrafidi/Documents/MATLAB/compEEG-data/behavioral/';
stimRoot = '/Users/nrafidi/Documents/MATLAB/compEEG-stim/';
fid = fopen([stimRoot 'KR_test.txt']);
swahili = textscan(fid, '%s');
swahili = swahili{1};

potentialResponses = [0 0 0 0; 0 0 0 1; 0 0 1 0; 0 0 1 1; 0 1 0 0; 0 1 0 1; ...
    0 1 1 0; 0 1 1 1; 1 0 0 0; 1 0 0 1; 1 0 1 0; 1 0 1 1; 1 1 0 0; 1 1 0 1; ...
    1 1 1 0; 1 1 1 1];
numPotResp = size(potentialResponses, 1);
krTrajByResponse = cell(numPotResp, 2);
for s = 1:numSub
    sub = subjects{s};
    
    load([behaveDataRoot '/' sub '/' sub '_answerTraj.mat']);
    load(['/Users/nrafidi/Documents/MATLAB/compEEG-data/results/' ...
        sub '/KRanalysis.mat']);
    
    krLabels = logical(krLabels);
    responseTraj(skippedItems,:) = []; %#ok<*SAGROW>
    subSwahili = swahili;
    subSwahili(skippedItems) = [];
    
    corrItems =find(krLabels);
    incItems = find(~krLabels);
    
    figure;
    plot(1:4, krTraj(corrItems(1),:), 'b');
    hold on;
    for i = 1:4
        if responseTraj(corrItems(1), i)
            scatter(i, krTraj(corrItems(1),i), 'MarkerEdgeColor', 'blue', ...
                'Marker', 'o', 'SizeData', 100);
        else
            scatter(i, krTraj(corrItems(1),i), 'MarkerEdgeColor', 'blue', ...
                'Marker', 'x', 'SizeData', 100);
        end
    end
    xlim([0.9, 4.1]);
    ylim([-0.1 1.1]);
    title(subSwahili{corrItems(1)});
    hold off;
    
    figure;
    plot(1:4, krTraj(incItems(2),:), 'r');
    hold on;
    for i = 1:4
        if responseTraj(incItems(2), i)
            scatter(i, krTraj(incItems(2),i), 'MarkerEdgeColor', 'red', ...
                'Marker', 'o', 'SizeData', 100);
        else
            scatter(i, krTraj(incItems(2),i), 'MarkerEdgeColor', 'red', ...
                'Marker', 'x', 'SizeData', 100);
        end
    end
    xlim([0.9, 4.1]);
    ylim([-0.1 1.1]);
    title(subSwahili{incItems(2)});
    hold off;
%     
%     for r = 1:numPotResp
%         indOfResp = ismember(responseTraj, potentialResponses(r,:), 'rows');
%         
%         krTrajByResponse{r, 1} = cat(1, krTrajByResponse{r, 1}, krTraj(indOfResp & krLabels, :));
%         krTrajByResponse{r, 2} = cat(1, krTrajByResponse{r, 2}, krTraj(indOfResp & ~krLabels, :));
%     end
    
end
% %%
% titles = {'None', 'R4', 'R3', 'R3, R4', 'R2', 'R2, R4', 'R2, R3', 'R2, R3, R4', ...
%     'R1', 'R1, R4', 'R1, R3', 'R1, R3, R4', 'R1, R2', ...
%     'R1, R2, R4', 'R1, R2, R3', 'R1, R2, R3, R4'};
% numPlots = sqrt(numPotResp);
% numTot = 0;
% for traj = 1:numPotResp
%     numCorrItems = size(krTrajByResponse{traj, 1}, 1);
%     numIncItems = size(krTrajByResponse{traj, 2}, 1);
%     numTot = numTot + numCorrItems + numIncItems;
%     titleForPlot = sprintf([titles{traj} ...
%         '\nC = %d I = %d'], numCorrItems, numIncItems);
%     
%     subplot(numPlots, numPlots, traj);
%     hold on;
%     if numCorrItems > 0
%         plot(1:4, mean(krTrajByResponse{traj, 1}, 1), 'b');
%         if numCorrItems > 1
%             errorbar(1:4, mean(krTrajByResponse{traj, 1}, 1), ...
%                 std(krTrajByResponse{traj, 1}, 1), 'b');
%         end
%     else
%         plot(1:4, zeros(1, 4), 'b');
%     end
%     
%     if numIncItems > 0
%         plot(1:4, mean(krTrajByResponse{traj, 2}, 1), 'r');
%         if numIncItems > 1
%             errorbar(1:4, mean(krTrajByResponse{traj, 2}, 1), ...
%                 std(krTrajByResponse{traj, 2}, 1), 'r');
%         end
%     else
%         plot(1:4, zeros(1, 4), 'r');
%     end
%     if traj == 13
%         legend({'Remembered', 'Forgotten'}, 'Location', 'north');
%     end
%     title(titleForPlot);
%     if traj == numPotResp
%         xlabel('Round');
%     end
%     if traj == 1
%         ylabel('Competition Classifier Output');
%     end
%     xlim([0.9 4.1]);
%     ylim([-0.1 1.1]);
% end

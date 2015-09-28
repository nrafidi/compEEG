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

titles = {'350-500ms All Rounds'};
for iTuple = 1:numTuples
    if numTuples > 1
    subplot(4, 4, iTuple);
    end
    hist(populationAccHist(:,iTuple));
    title(titles{iTuple});
    xlim([0.3, 0.8]);
    hold on;
    line([meanTrueAcc(iTuple) meanTrueAcc(iTuple)], ...
        [0 30]);
    if iTuple == numTuples
        legend({'Population', 'True'}, 'Location', 'eastoutside');
    end
    hold off;
end

% suptitle(sprintf('Bootstrap results for different tuples of KR data\nOnly Using Correctly Answered Instances'));
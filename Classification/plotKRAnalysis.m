% Process and plot KR analysis results

meanTrueAcc = mean(trueAcc, 1);
numDraws = length(populationAcc);

numTuples = size(trueAcc, 2);
populationAccHist = nan(numDraws, numTuples);

for iDraw = 1:numDraws    
    populationAccHist(iDraw,:) = mean(populationAcc{iDraw}, 1);
end


titles = {'R1', 'R2', 'R3', 'R4', 'R1, R2', 'R1, R3', 'R1, R4', ...
    'R2, R3', 'R2, R4', 'R3, R4', 'R1, R2, R3', 'R1, R2, R4', ...
    'R1, R3, R4', 'R2, R3, R4', 'R1, R2, R3, R4'};
for iTuple = 1:numTuples
    subplot(4, 4, iTuple);
    hist(populationAccHist(:,iTuple));
    title(titles{iTuple});
    xlim([0.4, 0.8]);
    hold on;
    line([meanTrueAcc(iTuple) meanTrueAcc(iTuple)], ...
        [0 30]);
    if iTuple == numTuples
        legend({'Population', 'True'}, 'Location', 'eastoutside');
    end
    hold off;
end

suptitle('Bootstrap results for different tuples of KR data');
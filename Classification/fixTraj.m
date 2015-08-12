krTraj = [];
krLabels = [];
indCorr = 1;
indInc = 1;
for i = 1:length(corrAnswers)
    if corrAnswers(i)
        if ~any(isnan(itemTrajCorr(indCorr,:)))
            krTraj = cat(1, krTraj, itemTrajCorr(indCorr,:));
            krLabels = cat(1, krLabels, 1);
        end
        indCorr = indCorr + 1;
    else
        if ~any(isnan(itemTrajInc(indInc,:)))
            krTraj = cat(1, krTraj, itemTrajInc(indInc,:));
            krLabels = cat(1, krLabels, 0);
        end
        indInc = indInc + 1;
    end
end
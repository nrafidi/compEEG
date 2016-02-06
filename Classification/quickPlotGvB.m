resultRoot = '/Users/nrafidi/Documents/MATLAB/compEEG-data/results/';

times = -100:20:760;
numT = length(times);
numF = 5;

load(sprintf('%sKR_Reg_slide.mat', resultRoot));
regAcc = accPerFold;
regWeights = weightVecs;

load(sprintf('%sKR_gVb_slide.mat', resultRoot));
gvbAcc = accPerFold;
gvbWeights = weightVecs;

minAcc = min(min([regAcc, gvbAcc]));
maxAcc = max(max([regAcc, gvbAcc]));

f1 = figure;
plot(times, mean(regAcc, 2));
hold on
plot(times, mean(gvbAcc, 2));
xlim([min(times), max(times)]);
ylim([minAcc-0.01, maxAcc+0.01]);
line([min(times), max(times)], [0.5, 0.5], 'LineStyle', '--', 'Color', 'k');
legend({'Remembered v Forgotten', 'Good Subject v Bad Subject', 'Chance'});
line([-40, -40], [minAcc - 0.01, maxAcc + 0.01], 'Color', 'k');
xlabel('Time relative to Stimulus Onset (ms)');
ylabel('Accuracy');
title(sprintf('Accuracy over Time\nMean Concat'));

saveas(f1, sprintf('%sfigures/gVb_slide.fig', resultRoot));
saveas(f1, sprintf('%sfigures/gVb_slide.pdf', resultRoot));

rdm_reg = nan(numT, numT, numF);
rdm_gvb = nan(numT, numT, numF);

for f = 1:numF
    rdm_reg(:,:,f) = squareform(pdist(cell2mat(regWeights(:,f)')'));
    rdm_gvb(:,:,f) = squareform(pdist(cell2mat(gvbWeights(:,f)')'));
end

f2 = figure;
imagesc(times, times, squeeze(mean(rdm_reg, 3)));
hold on;
line([min(times) max(times)], [-40, -40], 'Color', 'r');
line([-40, -40], [min(times) max(times)], 'Color', 'r');
colorbar;
xlabel('Time relative to stimulus onset (ms)');
ylabel('Time relative to stimulus onset (ms)');
title(sprintf('Average Weight RDM\nRemembered v Forgotten\nMean Concat'));
saveas(f2, sprintf('%sKR_reg_rdm.fig', resultRoot));
saveas(f2, sprintf('%sKR_reg_rdm.pdf', resultRoot));

f3 = figure;
imagesc(times, times, squeeze(mean(rdm_gvb, 3)));
hold on;
line([min(times) max(times)], [-40, -40], 'Color', 'r');
line([-40, -40], [min(times) max(times)], 'Color', 'r');
colorbar;
xlabel('Time relative to stimulus onset (ms)');
ylabel('Time relative to stimulus onset (ms)');
title(sprintf('Average Weight RDM\nGood Subject v Bad Subject\nMean Concat'));
saveas(f3, sprintf('%sKR_gvb_rdm.fig', resultRoot));
saveas(f3, sprintf('%sKR_gvb_rdm.pdf', resultRoot));
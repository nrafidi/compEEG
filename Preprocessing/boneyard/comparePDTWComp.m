% Compare PDTW to no PDTW (raw signal)

dataDir = '/Users/nrafidi/Documents/MATLAB/compEEG-data/results/';

fileName = [dataDir 'CompEEG_CV_Slide_Accs%s.mat'];

load(sprintf(fileName, '_PDTW'));

subAccs_p = subAccs;

load(sprintf(fileName, ''));

numSub = size(subAccs, 1);

for s = 1:numSub
    figure
    hold on
    plot(subAccs(s,:), 'b');
    plot(subAccs_p(s,:), 'r');
end
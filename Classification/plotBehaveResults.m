% Analyize and summarize behavioral data

% Full KR Analysis pipeline

% Starts with Competition and KR EEG data and produces KR population
% accuracy and histogram
%

isRepExp = true;
if isRepExp
    subjects = {'JJJ', 'III', 'KKK', 'BBB', 'GGG', 'HHH', 'AAA', 'CCC', 'DDD', 'EEE', 'FFF', 'MM', ...
        'OO', 'PP', 'QQ', 'RR', 'SS', 'TT', 'WW',...
        'YY'};
else
    subjects = {'AA', 'BB', 'DD', 'EE', 'F', 'GG', 'HH', 'JJ', ...
        'K', 'M', 'N', 'O', 'R', 'S', 'T', 'U', 'V', 'X', 'Y', 'Z'};
end
numSub = length(subjects);

if isRepExp
    behaveDataRoot = '/Users/nrafidi/Documents/MATLAB/compEEG-data-rep/behavioral/';
    resDir = '/Users/nrafidi/Documents/MATLAB/compEEG-data-rep/results/';
else
    behaveDataRoot = '/Users/nrafidi/Documents/MATLAB/compEEG-data/behavioral/';
    resDir = '/Users/nrafidi/Documents/MATLAB/compEEG-data/results/';
end

responseTrajList = [];
responseTrajList_C = [];
responseTrajList_I = [];
for s = 1:numSub
    sub = subjects{s};
    load([behaveDataRoot '/' sub '/' sub '_answerTraj.mat']);
    responseTrajList = cat(1, responseTrajList, responseTraj);
    
    load([resDir ...
        sub '/KRanalysis.mat']);
    krLabels = logical(krLabels);
    
    responseTrajList_C = cat(1, responseTrajList_C, responseTraj(krLabels,:));
    responseTrajList_I = cat(1, responseTrajList_I, responseTraj(~krLabels,:));
    
end

meanTraj = nanmean(responseTrajList, 1);
stdTraj = nanstd(responseTrajList, 1);

meanTraj_C = nanmean(responseTrajList_C, 1);
stdTraj_C = nanstd(responseTrajList_C, 1);

meanTraj_I = nanmean(responseTrajList_I, 1);
stdTraj_I = nanstd(responseTrajList_I, 1);

figure;
plot(1:4, meanTraj);
errorbar(1:4, meanTraj, stdTraj);
xlabel('Round');
ylabel(sprintf('Proportion correct answers\n(across subjects)'));

figure
plot(1:4, meanTraj_C);
errorbar(1:4, meanTraj_C, stdTraj_C);
hold on;
plot(1:4, meanTraj_I, 'r');
errorbar(1:4, meanTraj_I, stdTraj_I, 'r');
legend({'Remembered', 'Forgotten'})
xlabel('Round');
ylabel(sprintf('Proportion correct answers\n(across subjects)'));
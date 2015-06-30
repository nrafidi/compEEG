% Sorts KR predictions by final test performance (need to load manually
% into workspace bc fuck MATlAB)

% Need KR_sub_answer.csv

subjects = 'CE';
numSub = length(subjects);
fResPrefix = '/Users/nrafidi/Documents/MATLAB/compEEG-data/results/CompEEG_KR_';
KRanswer = importdata('/Users/nrafidi/Documents/MATLAB/compEEG-stim/KR_answer.txt');

for s = 1:numSub
    sub = subjects(s);
    load([fResPrefix sub '_itemTraj.mat']);
    
    subAnswers = KRsubanswer(s+1,:);
    numQ = length(subAnswers);
    corrAnswers = false(numQ, 1);
    for a = 1:numQ
        corrAnswers(a) = strcmpi(subAnswers{a}, KRanswer{a});
    end
    
    itemTrajCorr = itemTraj(corrAnswers, :);
    itemTrajInc = itemTraj(~corrAnswers, :);
    
    figure;
    plot(1:4, [nanmean(itemTrajCorr); nanmean(itemTrajInc)]);
    legend({'Correct', 'Incorrect'});
    
    save([fResPrefix sub '_itemTraj_Answers.mat'], 'itemTrajCorr', 'itemTrajInc', 'corrAnswers');
end
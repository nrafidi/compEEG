function [itemTrajCorr, itemTrajInc, corrAnswers] = sortKRTraj(sub)
% Sorts KR predictions by final test performance

fResPrefix = ['/Users/nrafidi/Documents/MATLAB/compEEG-data/results/' sub '/CompEEG__KR_'];
KRanswer = importdata('/Users/nrafidi/Documents/MATLAB/compEEG-stim/KR_answer.txt');

load([fResPrefix sub '_itemTraj.mat']);

if max(itemTraj) > 1
    itemTraj = exp(itemTraj)./(1 + exp(itemTraj));
end

load(['../../compEEG-data/results/answers/' sub '_answers.mat']);
numQ = length(answerList);
corrAnswers = false(numQ, 1);
for a = 1:numQ
    corrAnswers(a) = strcmpi(answerList{a}, KRanswer{a});
end

itemTrajCorr = itemTraj(corrAnswers, :);
itemTrajInc = itemTraj(~corrAnswers, :);

figure;
plot(1:4, [nanmean(itemTrajCorr); nanmean(itemTrajInc)]);
legend({'Correct', 'Incorrect'});

save([fResPrefix sub '_itemTraj_Answers_logit.mat'], 'itemTrajCorr', 'itemTrajInc', 'corrAnswers');

end
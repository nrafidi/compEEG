function corrAnswers = sortKRTrajBehave(sub, resDir)
% Sorts KR predictions by final test performance

fResPrefix = [resDir sub '/CompEEG__KR_'];
KRanswer = importdata('/Users/nrafidi/Documents/MATLAB/compEEG-stim/KR_answer.txt');

load([resDir 'answers/' sub '_answers.mat']);
numQ = length(answerList);
corrAnswers = false(numQ, 1);
for a = 1:numQ
    corrAnswers(a) = strcmpi(answerList{a}, KRanswer{a});
end


save([fResPrefix sub '_corrAnswers.mat'], 'corrAnswers');

end
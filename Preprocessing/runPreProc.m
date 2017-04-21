subjects = {'L', 'P', 'W', 'CC', 'FF', 'H', 'I', 'J', ...
'AA', 'BB', 'DD', 'EE', 'GG', 'HH', 'JJ', ...
        'X', 'Z', 'N'};
% Didn't do F, K, M, O, R, S, T, U, V, Y, N
numSub = length(subjects);
options.isVis = false;
for s = 17:numSub
    sub = subjects{s};
    disp(sub);
    preprocPipeline(sub, 'CompEEG', options);
    disp('CompEEG Complete');
%     preprocPipeline(sub, 'CompEEG__KR');
%     disp('CompEEG__KR Complete');
end
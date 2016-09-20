subjects = {'AA', 'BB', 'DD', 'F', 'EE', 'GG', 'HH', 'JJ', ...
        'K', 'M', 'O', 'R', 'S', 'T', 'U', 'V', 'X', 'Y', 'Z', 'N'};
numSub = length(subjects);

for s = 1:numSub
    sub = subjects{s};
    disp(sub);
    preprocPipeline(sub, 'CompEEG');
    disp('CompEEG Complete');
    preprocPipeline(sub, 'CompEEG__KR');
    disp('CompEEG__KR Complete');
end
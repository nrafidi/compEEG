subjects = {'AA', 'BB', 'DD', 'EE', 'F', 'GG', 'HH', 'JJ', ...
    'K', 'M', 'N', 'O', 'R', 'S', 'T', 'U', 'V', 'X', 'Y', 'Z'};

for s = 1:length(subjects)
    try
        reRunFeatures(subjects{s}, 'CompEEG');
%         reRunFeatures(subjects{s}, 'CompEEG__KR');
        fprintf(sprintf('Subject %s complete.\n', subjects{s}));
    catch
        fprintf(sprintf('Subject %s failed.\n', subjects{s}));
    end
end
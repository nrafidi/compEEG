% subjects = {'AA', 'BB', 'DD', 'EE', 'F', 'GG', 'HH', 'JJ', ...
%     'K', 'M', 'N', 'O', 'R', 'S', 'T', 'U', 'V', 'X', 'Y', 'Z'};

%Special fix needed for N, K, R
subjects = {'R'};
% subjects = {'Z', 'Y', 'X', 'V', 'U', 'T', 'S', 'R', 'O', 'F', 'M', 'JJ', ...
%     'HH', 'EE', 'DD'};


% subjects = {'F'};
% {'AA', 'FF', ...
%     'BB', 'DD', 'EE', 'GG', 'HH', 'JJ', 'M', 'N',...
%      'K', 'O', 'R', 'S', 'T', 'U', 'V', 'X', 'Y', 'Z'};
    
bpFinds = {[4, 8]};%, [8, 12], [12, 30]};
bpNames = {'theta'};%, 'alpha', 'beta'};
options = struct;
for s = 1:length(subjects)
    for b = 1:length(bpFinds)
        options.bpFind = bpFinds{b};
        options.bpName = bpNames{b};
        try
            preprocFreqDecomp(subjects{s}, 'CompEEG__KR', options);
            %         reRunFeatures(subjects{s}, 'CompEEG');
            %         reRunFeatures(subjects{s}, 'CompEEG__KR');
            fprintf(sprintf('Subject %s complete.\n', subjects{s}));
        catch
            fprintf(sprintf('Subject %s failed.\n', subjects{s}));
        end
    end
end
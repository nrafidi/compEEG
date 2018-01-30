subjects_rep = {'JJJ', 'III', 'KKK', 'BBB', 'GGG', 'HHH', 'AAA', ...
    'CCC', 'DDD', 'EEE', 'FFF', 'MM', ...
    'OO', 'PP', 'QQ', 'RR', 'SS', 'TT', 'WW',...
    'YY'};

subjects = {'AA', 'BB', 'DD', 'F', 'EE', 'GG', 'HH', 'JJ', ...
    'K', 'M', 'O', 'R', 'S', 'T', 'U', 'V', 'X', 'Y', 'Z', 'N'};

subjects_comp_only = {'L', 'P', 'W', 'CC', 'FF', 'H', 'I', 'J'};

data_dir = '/Users/nrafidi/Documents/MATLAB/compEEG-data/preproc-final/';
data_dir_rep = '/Users/nrafidi/Documents/MATLAB/compEEG-data-rep/preproc-final/';
final_fname_comp = '%s/%s/CompEEG_%s_Vis_BP2-200_N60_Ref_Hilbert-theta_Epochs_Features_Overlap_Time.mat';
final_fname_kr = '%s/%s/CompEEG__KR_%s_Vis_BP2-200_N60_Ref_Hilbert-theta_Epochs_Features_Overlap_Time.mat';

for i = 1:length(subjects_rep)
    sub = subjects_rep{i};
    if ~exist(sprintf(final_fname_comp, data_dir_rep, sub, sub), 'file')
        preprocFreqDecomp(sub, 'CompEEG', 1);
%         print('CompEEG Complete');
%         pause;
    end
    if ~exist(sprintf(final_fname_kr, data_dir_rep, sub, sub), 'file')
        preprocFreqDecomp(sub, 'CompEEG__KR', 1);
%         print('KR Complete');
%         pause;
    end
end
close all
for i = 1:length(subjects)
    sub = subjects{i};
    if ~exist(sprintf(final_fname_comp, data_dir, sub, sub), 'file')
        preprocFreqDecomp(sub, 'CompEEG', 0);
%         print('CompEEG Complete');
%         pause;
    end
    if ~exist(sprintf(final_fname_kr, data_dir, sub, sub), 'file')
        preprocFreqDecomp(sub, 'CompEEG__KR', 0);
%         print('KR Complete');
%         pause;
    end
end
close all
for i = 1:length(subjects_comp_only)
    sub = subjects_comp_only{i};
    if ~exist(sprintf(final_fname_comp, data_dir, sub, sub), 'file')
        preprocFreqDecomp(sub, 'CompEEG', 0);
%         print('CompEEG Complete');
%         pause;
    end
end
close all
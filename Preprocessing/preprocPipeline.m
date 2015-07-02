%% Input Parameters

% Subject ID
sub = 'C';
% Experiment Name
exp = 'CompEEG'; %alt: CompEEG_KR
% Path to data files
dataRoot = '/Users/nrafidi/Documents/MATLAB/compEEG-data/';
% High Pass Filter
HP = 2;
% Low Pass Filter
LP = 200;
% Notch Filter(s)
N = 60;
% Epochs
epochWin = [-0.3, 2];


%% EEGLab Pipeline

% Before running this pipeline, you must use the EEGLab GUI to make any
% rejections based on visual inspection. The expected filename will be of
% the form [exp '_' sub '_Vis.set']

% Load Visually Inspected Data
EEG = pop_loadset('filename',[exp '_' sub '_Vis.set'],'filepath',[dataRoot '/preproc-partial/']);
EEG = eeg_checkset( EEG );

% Load Channel Locations
EEG=pop_chanedit(EEG, 'load',{[dataRoot 'electrode_locs.xyz'] 'filetype' 'xyz'},'load',{[dataRoot 'electrode_locs.xyz'] 'filetype' 'xyz'},'save',[]);
EEG = eeg_checkset( EEG );
% Plot them to make sure they're correct
% figure; topoplot([],EEG.chanlocs, 'style', 'blank',  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo);
EEG = pop_saveset( EEG, 'filename',[exp '_' sub '_Vis_Locs.set'],'filepath',[dataRoot '/preproc-partial/']);
EEG = eeg_checkset( EEG );

% Band Pass Filter
EEG = pop_eegfiltnew(EEG, HP, LP, 846, 0, [], 1);
EEG.setname=[exp '_' sub '_Vis_Locs_BP' num2str(LP) '-' num2str(HP)];
EEG = eeg_checkset( EEG );
EEG = pop_saveset( EEG, 'filename',[EEG.setname '.set'],'filepath',[dataRoot '/preproc-partial/']);

% Notch Filter
EEG = pop_eegfiltnew(EEG, N - 5, N + 5, 846, 1, [], 1);
EEG.setname=[EEG.setname '_N' num2str(N)];
EEG = eeg_checkset( EEG );
EEG = pop_saveset( EEG, 'filename',[EEG.setname '.set'],'filepath',[dataRoot '/preproc-partial/']);

% Re-reference the electrodes to the group mean
EEG = pop_reref( EEG, [],'exclude',[63 65:72] );
EEG.setname=[EEG.setname '_Ref'];
EEG = eeg_checkset( EEG );
EEG = pop_saveset( EEG, 'filename',[EEG.setname '.set'],'filepath',[dataRoot '/preproc-partial/']);

% Parse into Epochs and remove Baseline
EEG = pop_epoch( EEG, {  }, epochWin, 'newname', [EEG.setname '_Epochs'], 'epochinfo', 'yes');
EEG = eeg_checkset( EEG );
EEG = pop_rmbase( EEG, [-300.7812             0]);
EEG.setname=[EEG.setname '_Base'];
oldSetName = EEG.setname;
EEG = pop_saveset( EEG, 'filename',[EEG.setname '.set'],'filepath',[dataRoot '/preproc-partial/']);

% Run ICA
EEG = pop_runica(EEG, 64, 'icatype', 'sobi');
EEG = eeg_checkset( EEG );
EEG = pop_saveset( EEG, 'savemode','resave');
EEG = eeg_checkset( EEG );

% Remove Blinks
EEG = pop_subcomp( EEG, [1 2], 0);
EEG.setname=[oldSetName '_ICA1-2'];
EEG = eeg_checkset( EEG );
EEG = pop_saveset( EEG, 'filename',[EEG.setname '.set'],'filepath',[dataRoot '/preproc-partial/']);

%% Create .mat file

if strcmp(exp, 'CompEEG')
    if sub < 'F'
        [data, labels, time] = getDataLabels_pilot(EEG);
    else
        [data, labels, time] = getDataLabels_Comp(EEG);
    end
else
    [data, labels, time] = getDataLabels_KR(EEG);
end

save([dataRoot '/preproc-final/' exp '_' sub '_Preprocessed.mat'], 'data', 'labels', 'time');
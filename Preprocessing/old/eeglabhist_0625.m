% EEGLAB history file generated on the 25-Jun-2015
% ------------------------------------------------

EEG = pop_biosig('/Users/nrafidi/Documents/MATLAB/CompEEG-data/CompEEG_C.bdf');
EEG.setname='CompEEG_C';
EEG = eeg_checkset( EEG );
EEG = pop_loadset('filename','CompEEG_C.set','filepath','/Users/nrafidi/Documents/MATLAB/CompEEG-data/');
EEG = eeg_checkset( EEG );
pop_eegplot( EEG, 1, 1, 1);
pop_eegplot( EEG, 1, 0, 1);
EEG.setname='CompEEG_C_Vis';
EEG = pop_loadset('filename','CompEEG_C_Vis.set','filepath','/Users/nrafidi/Documents/MATLAB/CompEEG-data/');
EEG = eeg_checkset( EEG );
EEG=pop_chanedit(EEG, 'lookup','/Users/nrafidi/Documents/MATLAB/Toolboxes/EEGLAB/plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp','load',[],'load',{'/Users/nrafidi/Documents/MATLAB/CompEEG-data/electrode_locs.xyz' 'filetype' 'xyz'},'save','/Users/nrafidi/Documents/MATLAB/CompEEG-data/electrode_locs_cartesian.ced');
EEG = eeg_checkset( EEG );
figure; topoplot([],EEG.chanlocs, 'style', 'blank',  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo);
EEG=pop_chanedit(EEG, 'load',{'/Users/nrafidi/Documents/MATLAB/CompEEG-data/electrode_locs.xyz' 'filetype' 'xyz'},'load',{'/Users/nrafidi/Documents/MATLAB/CompEEG-data/electrode_locs.xyz' 'filetype' 'xyz'},'save','/Users/nrafidi/Documents/MATLAB/CompEEG-data/eeg_locs_cartesian.ced');
EEG = eeg_checkset( EEG );
figure; pop_spectopo(EEG, 1, [0      1905900.4309], 'EEG' , 'percent', 15, 'freq', [6 10 22], 'freqrange',[2 25],'electrodes','off');
EEG = pop_saveset( EEG, 'filename','CompEEG_C_Vis_Locs.set','filepath','/Users/nrafidi/Documents/MATLAB/CompEEG-data/');
EEG = eeg_checkset( EEG );
EEG = pop_eegfiltnew(EEG, 2, 200, 846, 0, [], 1);
EEG.setname='CompEEG_C_Vis_Locs_BP2-200';
EEG = eeg_checkset( EEG );
EEG = pop_eegfiltnew(EEG, 55, 65, 846, 1, [], 1);
EEG.setname='CompEEG_C_Vis_Locs_BP2-200_N60';
EEG = eeg_checkset( EEG );
EEG = pop_reref( EEG, [],'exclude',[65:72] );
EEG.setname='CompEEG_C_Vis_Locs_BP2-200_N60_Ref';
EEG = eeg_checkset( EEG );
EEG = pop_epoch( EEG, {  }, [-0.3           3], 'newname', 'CompEEG_C_Vis_Locs_BP2-200_N60_Ref_Epochs', 'epochinfo', 'yes');
EEG = eeg_checkset( EEG );
EEG = pop_rmbase( EEG, [-300.7812             0]);
EEG.setname='CompEEG_C_Vis_Locs_BP2-200_N60_Ref_Epochs_Base';
EEG = eeg_checkset( EEG );

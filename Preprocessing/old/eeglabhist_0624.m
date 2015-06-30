% EEGLAB history file generated on the 24-Jun-2015
% ------------------------------------------------
% Load BDF
EEG = pop_biosig('/Users/nrafidi/Documents/MATLAB/CompEEG-data/CompEEG_C.bdf');
EEG.setname='CompEEG_C';
EEG = eeg_checkset( EEG );
% Save as .set and reload
EEG = pop_loadset('filename','CompEEG_C.set','filepath','/Users/nrafidi/Documents/MATLAB/CompEEG-data/');
EEG = eeg_checkset( EEG );
% Visually inspecting the data
pop_eegplot( EEG, 1, 1, 1);
pop_eegplot( EEG, 1, 0, 1);
EEG.setname='CompEEG_C_Vis';
EEG = pop_loadset('filename','CompEEG_C_Vis.set','filepath','/Users/nrafidi/Documents/MATLAB/CompEEG-data/');
EEG = eeg_checkset( EEG );
pop_eegplot( EEG, 1, 1, 1);
% Visually inspecting the frequency domain
figure; pop_spectopo(EEG, 1, [0      1905900.4309], 'EEG' , 'percent', 15, 'freqrange',[1 200],'electrodes','off');
% Band Pass Filtering the data between 2 and 200Hz
EEG = pop_eegfiltnew(EEG, 2, 200, 846, 0, [], 1);
EEG.setname='CompEEG_C_Vis_BP2-200';
EEG = eeg_checkset( EEG );
% Notch Filtering the data at 60Hz
figure; pop_spectopo(EEG, 1, [0      1905900.4309], 'EEG' , 'percent', 15, 'freqrange',[0 250],'electrodes','off');
EEG = pop_eegfiltnew(EEG, 55, 65, 846, 1, [], 1);
EEG.setname='CompEEG_C_Vis_BP2-200_notch60';
EEG = eeg_checkset( EEG );
% Inspecting Post Filtering
figure; pop_spectopo(EEG, 1, [0      1905900.4309], 'EEG' , 'percent', 15, 'freqrange',[2 250],'electrodes','off');
% Re-referencing the data to the mean of all SCALP electrodes
EEG = pop_reref( EEG, [],'exclude',[65:72] );
EEG.setname='CompEEG_C_Vis_BP2-200_Notch60_Ref';
EEG = eeg_checkset( EEG );
pop_eegplot( EEG, 1, 1, 1);
pop_eegplot( EEG, 1, 1, 1);
% Parsed Data into epochs
EEG = pop_epoch( EEG, {  }, [-0.3           3], 'newname', 'CompEEG_C_Vis_BP2-200_Notch60_Ref_epochs', 'epochinfo', 'yes');
EEG = eeg_checkset( EEG );
% Removed Baseline activity
EEG = pop_rmbase( EEG, [-300.7812             0]);
EEG.setname='CompEEG_C_Vis_BP2-200_Notch60_Ref_epochs_base';
EEG = eeg_checkset( EEG );
EEG = eeg_checkset( EEG );
% Generated ICA components using sobi
EEG = pop_runica(EEG, 72);
EEG = eeg_checkset( EEG );
% Visually inspected the components and rejected component 2
pop_eegplot( EEG, 0, 1, 1);
EEG = pop_subcomp( EEG, [2], 0);
EEG.setname='CompEEG_C_Vis_BP2-200_Notch60_Ref_epochs_base_ICA2';
EEG = eeg_checkset( EEG );
% Checking for abnormal ICA component values
EEG = pop_eegthresh(EEG,0,[1:71] ,-25,25,-0.30078,2.998,2,0);
EEG = eeg_checkset( EEG );
% Checking for abnormal trends in ICA component values
EEG = pop_rejtrend(EEG,0,[1:71] ,1690,50,0.3,0,0);
EEG = eeg_checkset( EEG );
% We couldn't reject the marks (bug) so we stored them
function [featData, labels, winTime, options, featOptions] = preprocPipeline(sub, experi, varargin)
% preprocPipeline: a preprocessing pipeline for CompEEG and CompEEG__KR
% data. Runs a user-selected set of preprocessing procedures (set with the
% optional struct options) and then extracts the features for use in
% classification. Currently uses feature extraction defaults. Note that all
% intermediate files are saved.
%
% Inputs:
%   sub: the subject to be processed, e.g. 'AA'
%   experi: the experiment to be processed, either 'CompEEG' or
%   'CompEEG__KR'
%   options: (optional) sets the parameters to use when preprocessing,
%   which include:
%       isVis: set to true if using visually inspected data
%       HP: the edge of the high-pass filter to apply (nan if none)
%       LP: the edge of the low-pass filter to apply (nan if none)
%       N: the center of the notch filter to apply (nan if none)
%       doRef: set to true to rereference electrodes to the group mean
%       (recommended)
%       doBase: baseline correct using prestimulus data (recommended)
%       runICA: compute ICA weights and subtract to remove blinks
%
% Outputs:
%   featData: the data in the form of samples x features
%   labels: the labels corresponding to these data
%   options: the basic preprocessing options used
%   featOptions: the feature extraction options used (currently default)

clear EEG;
eeglab;

% Path to data files
dataRoot = '/Users/nrafidi/Documents/MATLAB/compEEG-data-rep/';
% Intermediate files
saveInterDir = [dataRoot '/preproc-partial/' sub '/'];
% Final files
preProcFinalRoot = [dataRoot '/preproc-final/' sub '/'];
if ~exist(preProcFinalRoot, 'dir');
    mkdir(preProcFinalRoot);
end
saveFname = [preProcFinalRoot experi '_' sub];
% Epochs
epochWin = [-0.3, 2];

if nargin > 2
    options = varargin{1};
else
    options = struct;
end

% Visual Inspection
if ~isfield(options, 'isVis')
    options.isVis = true;
end
% High Pass Filter
if ~isfield(options, 'HP')
    options.HP = 2;
end
% Low Pass Filter
if ~isfield(options, 'LP')
    options.LP = 200;
end
% Notch Filter(s)
if ~isfield(options, 'N')
    options.N = 60;
end
% Re-reference electrodes to group mean
if ~isfield(options, 'doRef')
    options.doRef = true;
end
% Baseline Correct Epochs
if ~isfield(options, 'doBase')
    options.doBase = true;
end
% Do ICA from scratch
if ~isfield(options, 'runICA')
    options.runICA = true;
end
% Apply pre-computed ICA weights - haven't implemented yet
if ~isfield(options, 'useICA')
    options.useICA = false;
end

%% EEGLab Pipeline

% Before running this pipeline, you should use the EEGLab GUI to make any
% rejections based on visual inspection. The expected filename will be of
% the form [exp '_' sub  '_Vis.set']

fnameRoot = [experi '_' sub];

if options.isVis
    fnameRoot = [fnameRoot '_Vis'];
    saveFname = [saveFname '_Vis'];
end


% Load Visually Inspected Data
EEG = pop_loadset('filename',[fnameRoot '.set'],'filepath', saveInterDir);
EEG = eeg_checkset( EEG );
runningFname = [EEG.setname];

% % Load Channel Locations - currently not working
% EEG=pop_chanedit(EEG, 'load',{[dataRoot 'electrode_locs.xyz'] 'filetype' 'xyz'},'load',{[dataRoot 'electrode_locs.xyz'] 'filetype' 'xyz'},'save',[]);
% EEG = eeg_checkset( EEG );
% % Plot them to make sure they're correct
% % figure; topoplot([],EEG.chanlocs, 'style', 'blank',  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo);
% EEG.setname=[EEG.setname '_Locs'];
% EEG = pop_saveset( EEG, 'filename',[EEG.setname '.set'],'filepath',[dataRoot '/preproc-partial/']);
% EEG = eeg_checkset( EEG );

% Band Pass Filter
if ~isnan(options.HP) && ~isnan(options.LP) && ...
        ~exist([saveInterDir runningFname '_BP' num2str(options.HP) '-' num2str(options.LP) '.set'], 'file');
    fprintf('BP filtering\n')
    tic
    EEG = pop_eegfiltnew(EEG, options.HP, options.LP, 846, 0, [], 1);
    toc
    EEG.setname=[EEG.setname '_BP' num2str(options.HP) '-' num2str(options.LP)];
    saveFname = [saveFname '_BP' num2str(options.HP) '-' num2str(options.LP)];
    EEG = eeg_checkset( EEG );
%     EEG = pop_saveset( EEG, 'filename',[EEG.setname '.set'],'filepath', saveInterDir);
elseif ~isnan(options.HP) && ~isnan(options.LP)
    saveFname = [saveFname '_BP' num2str(options.HP) '-' num2str(options.LP)];
    runningFname = [runningFname '_BP' num2str(options.HP) '-' num2str(options.LP)];
    EEG = pop_loadset('filename',[runningFname '.set'],'filepath', saveInterDir);
end

% Notch Filter
if ~isnan(options.N) && ...
        ~exist([saveInterDir runningFname '_N' num2str(options.N) '.set'], 'file');
    fprintf('Notch filtering\n');
    tic
    EEG = pop_eegfiltnew(EEG, options.N - 5, options.N + 5, 846, 1, [], 1);
    toc
    EEG.setname=[EEG.setname '_N' num2str(options.N)];
    saveFname = [saveFname '_N' num2str(options.N)];
    EEG = eeg_checkset( EEG );
%     EEG = pop_saveset( EEG, 'filename',[EEG.setname '.set'],'filepath', saveInterDir);
elseif ~isnan(options.N)
    saveFname = [saveFname '_N' num2str(options.N)];
    runningFname = [runningFname '_N' num2str(options.N)];
    EEG = pop_loadset('filename',[runningFname '.set'],'filepath', saveInterDir);
end

% Re-reference the electrodes to the group mean
if options.doRef && ~exist([saveInterDir runningFname '_Ref.set'], 'file')
    fprintf('Re-referencing\n');
    tic
    EEG = pop_reref( EEG, [],'exclude',[63 65:72] );
    toc
    EEG.setname=[EEG.setname '_Ref'];
    saveFname = [saveFname '_Ref'];
    EEG = eeg_checkset( EEG );
%     EEG = pop_saveset( EEG, 'filename',[EEG.setname '.set'],'filepath', saveInterDir);
elseif options.doRef
    saveFname = [saveFname '_Ref'];
    runningFname = [runningFname '_Ref'];
    EEG = pop_loadset('filename',[runningFname '.set'],'filepath', saveInterDir);
end

% Parse into Epochs and remove Baseline
if options.doBase
    if ~exist([saveInterDir runningFname '_Epochs_Base.set'], 'file')
        EEG = pop_epoch( EEG, {  }, epochWin, 'newname', [EEG.setname '_Epochs'], 'epochinfo', 'yes');
        EEG = eeg_checkset( EEG );
        EEG = pop_rmbase( EEG, [-300.7812             0]);
        EEG.setname=[EEG.setname '_Epochs_Base'];
        saveFname = [saveFname '_Epochs_Base'];
        
        oldSetName = EEG.setname;
%         EEG = pop_saveset( EEG, 'filename',[EEG.setname '.set'],'filepath', saveInterDir);
    else
        saveFname = [saveFname '_Epochs_Base'];
        runningFname = [runningFname '_Epochs_Base'];
        EEG = pop_loadset('filename',[runningFname '.set'],'filepath', saveInterDir);
    end
else
    if ~exist([saveInterDir runningFname '_Epochs.set'], 'file')
        EEG = pop_epoch( EEG, {  }, epochWin, 'newname', [EEG.setname '_Epochs'], 'epochinfo', 'yes');
        EEG = eeg_checkset( EEG );
        saveFname = [saveFname '_Epochs'];
        oldSetName = EEG.setname;
%         EEG = pop_saveset( EEG, 'filename',[EEG.setname '.set'],'filepath', saveInterDir);
    else
        runningFname = [runningFname '_Epochs'];
        saveFname = [saveFname '_Epochs'];
        EEG = pop_loadset('filename',[runningFname '.set'],'filepath', saveInterDir);
    end
end


% Run ICA
if options.runICA && ~exist([saveInterDir runningFname '_ICA1-2.set'], 'file')
    %     keyboard;
    fprintf('Running ICA\n');
    tic
    EEG = pop_runica(EEG, 64, 'icatype', 'sobi');
    toc
    EEG = eeg_checkset( EEG );
%     EEG = pop_saveset( EEG, 'savemode','resave');
    EEG = eeg_checkset( EEG );
elseif options.useICA
    %***Need to figure out how to apply pre-computed weights
end

% Remove Blinks
if (options.runICA || options.useICA) && ...
        ~exist([saveInterDir runningFname '_ICA1-2.set'], 'file')
    fprintf('Removing Blinks\n');
    tic
    EEG = pop_subcomp( EEG, [1 2], 0);
    toc
    EEG.setname=[oldSetName '_ICA1-2'];
    saveFname = [saveFname '_ICA1-2'];
    EEG = eeg_checkset( EEG );
    EEG = pop_saveset( EEG, 'filename',[EEG.setname '.set'],'filepath',saveInterDir);
elseif options.runICA
    runningFname = [runningFname '_ICA1-2'];
    saveFname = [saveFname '_ICA1-2'];
    EEG = pop_loadset('filename',[runningFname '.set'],'filepath', saveInterDir);
end

% Create .mat file

if strcmp(experi, 'CompEEG')
    if any(sub < 'F') && length(sub) == 1
        [data, labels, time] = getDataLabels_pilot(EEG); %#ok<*ASGLU>
    else
        [data, labels, time] = getDataLabels_Comp(EEG);
    end
else
    [data, labels, time] = getDataLabels_KR(sub, EEG);
end
%
disp(sum(labels(:,1) == 1));
bar(labels(:,1));
title(sub);
% keyboard;

preProcOptions = options; %#ok<*NASGU>
save([saveFname '.mat'], 'data', 'labels', 'time', 'preProcOptions');

% Extract Relevant Features

[featData, labels, winTime, featOptions] = extractFeatures([saveFname '.mat']);

save([saveFname '_Features_Overlap_Time.mat'], 'featData', 'labels', 'featOptions', 'winTime');

end
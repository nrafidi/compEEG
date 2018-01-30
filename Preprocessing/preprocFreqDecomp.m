function [featData, labels, options, featOptions] = ...
    preprocFreqDecomp(sub, experi, isRep, varargin)
% preprocFreqDecomp: a preprocessing pipeline for CompEEG and CompEEG__KR
% data. Runs a user-selected set of preprocessing procedures (set with the
% optional struct options) and then extracts the features for use in
% classification. Runs a Hilbert Transform and selects the given frequency
% window (default theta band) for features.
%
% Inputs:
%   sub: the subject to be processed, e.g. 'AA'
%   experi: the experiment to be processed, either 'CompEEG' or
%   'CompEEG__KR'
%   isRep: is the subject from the original or the replication experiment
%   options: (optional) sets the parameters to use when preprocessing,
%   which include:
%       isVis: set to true if using visually inspected data
%       HP: the edge of the high-pass filter to apply (nan if none)
%       LP: the edge of the low-pass filter to apply (nan if none)
%       N: the center of the notch filter to apply (nan if none)
%       doRef: set to true to rereference electrodes to the group mean
%       (recommended)
%       bpfind: frequency boundaries to select
%       bpname: name of frequency band (e.g. 'theta')
%
% Outputs:
%   featData: the data in the form of samples x features
%   labels: the labels corresponding to these data
%   options: the basic preprocessing options used
%   featOptions: the feature extraction options used (currently default)

clear EEG;
eeglab;

% Path to data files
if isRep
    dataRoot = '/Users/nrafidi/Documents/MATLAB/compEEG-data-rep/';
else
    dataRoot = '/Users/nrafidi/Documents/MATLAB/compEEG-data/';
end
% Intermediate files
saveInterDir = [dataRoot '/preproc-partial/' sub '/'];
% Final files
saveFname = [dataRoot '/preproc-final/' sub '/' experi '_' sub];
% Epochs
epochWin = [-0.3, 2];

if nargin > 3
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
if ~isfield(options, 'bpFind')
    options.bpFind = [4, 8];
end
if ~isfield(options, 'bpName')
    options.bpName = 'theta';
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
    runningFname = [runningFname '_BP' num2str(options.HP) '-' num2str(options.LP)];
    saveFname = [saveFname '_BP' num2str(options.HP) '-' num2str(options.LP)];
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
    runningFname = [runningFname '_N' num2str(options.N)];
    saveFname = [saveFname '_N' num2str(options.N)];
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
    runningFname = [runningFname '_Ref'];
    saveFname = [saveFname '_Ref'];
    EEG = pop_loadset('filename',[runningFname '.set'],'filepath', saveInterDir);
end

% Hilbert transform
if ~any(isnan(options.bpFind))
    runningFname = [runningFname '_Hilbert-' options.bpName];
    hilpow = hilbertTransform(EEG, options.bpFind);
    EEG.data = squeeze(hilpow);
    EEG.setname=[EEG.setname '_Hilbert-' options.bpName];
    saveFname = [saveFname '_Hilbert-' options.bpName];
%     EEG = pop_saveset( EEG, 'filename',[EEG.setname '.set'],'filepath', saveInterDir);
end

% Parse into Epochs and remove Baseline
if ~exist([saveInterDir runningFname '_Epochs.set'], 'file')
    EEG = pop_epoch( EEG, {  }, epochWin, 'newname', [EEG.setname '_Epochs'], 'epochinfo', 'yes');
    EEG = eeg_checkset( EEG );
    EEG.setname=[EEG.setname '_Epochs'];
    saveFname = [saveFname '_Epochs'];
    oldSetName = EEG.setname;
%     EEG = pop_saveset( EEG, 'filename',[EEG.setname '.set'],'filepath', saveInterDir);
else
    runningFname = [runningFname '_Epochs'];
    saveFname = [saveFname '_Epochs'];
    EEG = pop_loadset('filename',[runningFname '.set'],'filepath', saveInterDir);
end

% Create .mat file

if strcmp(experi, 'CompEEG')
    [data, labels, time] = getDataLabels_Comp(EEG); %#ok<*ASGLU>
else
    [data, labels, time] = getDataLabels_KR(sub, EEG);
end
%%
figure;
bar(labels(:,1));
title(sprintf('%s = %d\n', sub, sum(labels(:,1) == 1)));


preProcOptions = options; %#ok<*NASGU>
save([saveFname '.mat'], 'data', 'labels', 'time', 'preProcOptions');

% Extract Relevant Features

[featData, labels, winTime, featOptions] = extractFeatures([saveFname '.mat']);

save([saveFname '_Features_Overlap_Time.mat'], 'featData', 'labels', 'winTime', 'featOptions');

end
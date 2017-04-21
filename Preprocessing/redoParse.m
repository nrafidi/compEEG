% Script to redo parsing/feature extraction for all subjects

subjects = {'YY', 'WW', 'TT', 'GGG', 'HHH'};
experi = 'CompEEG__KR';

% eeglab;

% Path to data files
dataRoot = '/Users/nrafidi/Documents/MATLAB/compEEG-data-rep/';
options = struct;
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
% Do ICA from scratch
if ~isfield(options, 'runICA')
    options.runICA = true;
end
% Apply pre-computed ICA weights - haven't implemented yet
if ~isfield(options, 'useICA')
    options.useICA = false;
end



for s = 4:length(subjects)
    sub = subjects{s};
    clear EEG;
    
    
    % Intermediate files
    saveInterDir = [dataRoot '/preproc-partial/' sub '/'];
    % File to Load
    loadFname = [experi '_' sub];
    % Final files
    saveFname = [dataRoot '/preproc-final/' sub '/' experi '_' sub];
    
    
    if options.isVis
        loadFname = [loadFname '_Vis'];
        saveFname = [saveFname '_Vis']; %#ok<*AGROW>
    end
    
    if ~isnan(options.HP) && ~isnan(options.LP)
        saveFname = [saveFname '_BP' num2str(options.HP) '-' num2str(options.LP)];
        loadFname = [loadFname '_BP' num2str(options.HP) '-' num2str(options.LP)];
    end
    
    if ~isnan(options.N)
        saveFname = [saveFname '_N' num2str(options.N)];
        loadFname = [loadFname '_N' num2str(options.N)];
    end
    
    if options.doRef
        saveFname = [saveFname '_Ref'];
        loadFname = [loadFname '_Ref'];
    end
    
    saveFname = [saveFname '_Epochs_Epochs_Base'];
    loadFname = [loadFname '_Epochs_Epochs_Base'];
    
    if options.runICA || options.useICA
        saveFname = [saveFname '_ICA1-2'];
        loadFname = [loadFname '_ICA1-2'];
    end
    
    
    EEG = pop_loadset( 'filename',[loadFname '.set'],'filepath',saveInterDir);
    
    % Create .mat file
    
    if strcmp(experi, 'CompEEG')
        [data, labels, time] = getDataLabels_Comp(sub, EEG);
    else
        [data, labels, time] = getDataLabels_KR(sub, EEG);
    end
    disp(sum(labels(:,1) == 1))
    bar(labels(:,1));
    keyboard;
    preProcOptions = options;  %#ok<*NASGU>
    save([saveFname '.mat'], 'data', 'labels', 'time', 'preProcOptions');
    
    % Extract Relevant Features
    
    [featData, labels, featOptions] = extractFeatures([saveFname '.mat']);
    
    save([saveFname '_Features.mat'], 'featData', 'labels', 'featOptions');
    
end
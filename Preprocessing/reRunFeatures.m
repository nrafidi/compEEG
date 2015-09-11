function [featData, labels, options, featOptions] = reRunFeatures(sub, experi, varargin)
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
%       runICA: compute ICA weights and subtract to remove blinks
%
% Outputs:
%   featData: the data in the form of samples x features
%   labels: the labels corresponding to these data
%   options: the basic preprocessing options used
%   featOptions: the feature extraction options used (currently default)

% Path to data files
dataRoot = '/Users/nrafidi/Documents/MATLAB/compEEG-data/';
% Final files
loadFname = [dataRoot '/preproc-final/' sub '/' experi '_' sub];


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
% Do ICA from scratch
if ~isfield(options, 'runICA')
    options.runICA = true;
end
% Apply pre-computed ICA weights - haven't implemented yet
if ~isfield(options, 'useICA')
    options.useICA = false;
end

if options.isVis
    loadFname = [loadFname '_Vis'];
end

% Band Pass Filter
if ~isnan(options.HP) && ~isnan(options.LP)
    loadFname = [loadFname '_BP' num2str(options.HP) '-' num2str(options.LP)];
end

% Notch Filter
if ~isnan(options.N)
    loadFname = [loadFname '_N' num2str(options.N)];
end

% Re-reference the electrodes to the group mean
if options.doRef
    loadFname = [loadFname '_Ref'];
end

% Parse into Epochs and remove Baseline
loadFname = [loadFname '_Epochs_Base'];

% Remove Blinks
if (options.runICA || options.useICA)
    loadFname = [loadFname '_ICA1-2'];
end

% Extract Relevant Features
featOptions = struct;
featOptions.overLap = true;
[featData, labels, featOptions] = extractFeatures([loadFname '.mat'], featOptions);


save([loadFname '_Features_Overlap.mat'], 'featData', 'labels', 'featOptions');

end
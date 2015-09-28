function [featData, labels, options] = extractFeatures(loadFname, varargin) %#ok<*STOUT>
% extractFeatures: a function for extracting features from EEG data
% recorded in the file loadFname, with optional input parameters. Note that
% this is intended to be called from preprocPipeline.m
%
% Inputs:
%   loadFname: the filename containing the preprocessed EEG data. It should
%   contain a variable data that is samples x channels x time, a variable
%   labels, and a variable time
%   options: (optional) sets the parameters to use when extracting, which
%   include:
%       erpWinSize: the size of window to average
%       specWinSize: the size of window to use for spectral decomposition
%       specStartFrequency: the start of each frequency bin to extract
%       specEndFrequency: the end of each frequency bin (should be same
%       length as previous)
%       maxTime: the maximum amount of time post-onset to consider for
%       extraction.
%
% Outputs:
%   featData: the data in the form of samples x features
%   labels: the labels corresponding to these data
%   options: the feature extraction options used

% Parameters
if nargin > 1
    options = varargin{1};
else
    options = struct;
end
% Window size to average voltages (ms)
if ~isfield(options, 'erpWinSize')
    options.erpWinSize = 50;
end
% Window size for spectral decomposition (ms)
if ~isfield(options, 'specWinSize')
    options.specWinSize = 200;
end
% Overlap windows?
if ~isfield(options, 'overLap')
    options.overLap = false;
end
% Frequencies to include (lower edge, Hz)
if ~isfield(options, 'specStartFreq')
    options.specStartFreq = [4, 9, 13];
end
% Frequencies to include (upper edge, Hz)
if ~isfield(options, 'specEndFreq')
    options.specEndFreq = [8, 12, 20];
end
% Maximum time post-stim to consider (ms)
if ~isfield(options, 'maxTime')
    options.maxTime = 800;
end

load(loadFname);

numSamp = size(data, 1); %#ok<*NODEF>

if ~isfield(options, 'minTime')
    minTime = min(time(time >= 0));
else
    minTime = options.minTime;
end

if options.overLap
    erpWin = minTime:20:options.maxTime;
    %     specWin = minTime:20:options.maxTime;
else
    erpWin = minTime:options.erpWinSize:options.maxTime;
    %     specWin = minTime:options.specWinSize:options.maxTime;
end

numErp = length(erpWin);
% numSpec = length(specWin);

featData = [];

for e = 2:numErp
    erp = (time >= erpWin(e-1)) & (time < (erpWin(e-1)+options.erpWinSize));
    newData = squeeze(mean(data(:,:,erp), 3));
    featData = cat(2, featData, newData);
end

% for p = 2:numSpec
%     for sp = 1:length(options.specStartFreq)
%         spec = (time >= specWin(p-1)) & (time < (specWin(p-1)+options.specWinSize));
%         newData = [];
%         for n = 1:numSamp
%             dataToTrans = squeeze(data(n,:,spec));
%             trans = abs(fft(dataToTrans'));
%             freqData = mean(trans(options.specStartFreq(sp):options.specEndFreq(sp), :));
%             newData = cat(1, newData, freqData);
%         end
%         featData = cat(2, featData, newData);
%     end
% end

end
function [data, labels, time] = getDataLabels_KR(EEG)
% getDataLabels_Comp: converts the eeglab data structure EEG into a more
% workable format. Designed for CompEEG__KR.
%
% Inputs:
%   EEG: the data structure to be converted
%
% Outputs:
%   data: the data as a double array of size samples x channels x time
%   labels: the labels for the data
%   time: a vector corresponding the time dimension of data with
%   stimulus-relative timepoints

% Extract the epochs
epochInfo = EEG.epoch;
numEpochs = length(epochInfo);

labelInds = false(numEpochs, 1);
stimLabels = nan(numEpochs, 2); %Study/Quiz binary, label identity

oddblock = false;
prevEvent = str2num(epochInfo(1).eventtype);%#ok<*ST2NM>
for i = 2:numEpochs
    currEvent = str2num(epochInfo(i).eventtype);
    
    if currEvent == 1 && prevEvent == 1
        oddblock = ~oddblock;
    end
    
    if currEvent > 1
        labelInds(i) = true;
        stimLabels(i, 2) = currEvent;
        stimLabels(i, 1) = double(~oddblock);
    end
    prevEvent = currEvent;
end

labels = stimLabels(labelInds, :);
data = EEG.data(:,:,labelInds);
data = permute(data, [3 1 2]);
data = data(:, 1:64, :);
time = EEG.times;

end

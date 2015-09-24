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
% epochInfo = [epochInfo(1:62) epochInfo(65:end)]; % for sub N
numEpochs = length(epochInfo);

labelInds = false(numEpochs, 1);
stimLabels = nan(numEpochs, 2); %Study/Quiz binary, label identity
% 
if ischar(epochInfo(1).eventtype)
    firstEvent = str2num(epochInfo(1).eventtype);%#ok<*ST2NM>
else
    firstEvent = epochInfo(1).eventtype;
end
% if firstEvent == 1 || firstEvent == 255
%     oddblock = false;
% else
%     oddblock = true;
% end
oddblock = true;
collectLabels = [];
prevprevprevEvent = nan;
prevprevEvent = nan;
prevEvent = firstEvent;
for i = 2:numEpochs
    if ischar(epochInfo(i).eventtype);
        currEvent = str2num(epochInfo(i).eventtype);
    else
        currEvent = epochInfo(i).eventtype;
    end
%     if i == 128
%         keyboard;
%     end
    if (currEvent == 255) && (prevprevEvent ~= 255) && (prevprevprevEvent ~=255)
        collectLabels = [];
        if stimLabels(i-1, 2) == 1;
            labelInds(i-1) = false;
        end
        oddblock = ~oddblock;
%         keyboard;
    else
        labelInds(i) = true;
        collectLabels = cat(1, collectLabels, currEvent);
        
%         if length(unique(collectLabels)) < (length(collectLabels) - 1)
%             oddblock = ~oddblock;
%             collectLabels = [];
% %             keyboard;
%         end
%         
        stimLabels(i, 2) = currEvent;
        stimLabels(i, 1) = double(oddblock);
    end
    
    prevprevprevEvent = prevprevEvent;
    prevprevEvent = prevEvent;
    prevEvent = currEvent;
    
end

labels = stimLabels(labelInds, :);
data = EEG.data(:,:,labelInds);
data = permute(data, [3 1 2]);
data = data(:, 1:64, :);
time = EEG.times;

end

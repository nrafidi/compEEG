function [data, labels, time] = getDataLabels_KR(EEG)

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

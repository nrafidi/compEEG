function [data, labels, time] = getDataLabels_Comp(EEG)

epochInfo = EEG.epoch;

numEpochs = length(epochInfo);

labelInds = false(numEpochs, 1);
stimLabels = nan(numEpochs, 2); %Comp/nonComp binary, category identity

labelStart = 19;

for i = 1:numEpochs
    if iscell(epochInfo(i).eventtype)
        epochInfo(i).eventtype = epochInfo(i).eventtype{1};
    end
    label = str2num(epochInfo(i).eventtype);%#ok<*ST2NM>
    
    if (label > labelStart) && (label ~= 255)
        labelInds(i) = true;
        stimLabels(i, 2) = label;
        if label > labelStart*2
            stimLabels(i, 1) = 0;
        else
            stimLabels(i, 1) = 1;
        end
    end
end

labels = stimLabels(labelInds, :);
data = EEG.data(:,:,labelInds);
data = permute(data, [3 1 2]);
data = data(:, 1:64, :);
time = EEG.times;

end

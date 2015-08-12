function [data, labels, time] = getDataLabels_pilot(EEG)

epochInfo = EEG.epoch;

numEpochs = length(epochInfo);

labelInds = false(numEpochs, 1);
stimLabels = nan(numEpochs, 2); %Comp/nonComp binary, label identity

% Start after the second set of instructions
firstInstSeen = false;
secondInstSeen = false;
thirdInstSeen = false;
fourthInstSeen = false;
for i = 1:numEpochs
    if iscell(epochInfo(i).eventtype)
        epochInfo(i).eventtype = epochInfo(i).eventtype{1};
    end
    if ~isstr(epochInfo(i).eventtype)
        epochInfo(i).eventtype = num2str(epochInfo(i).eventtype);
    end
    if fourthInstSeen
        label = str2num(epochInfo(i).eventtype);%#ok<*ST2NM>
        
%         keyboard;
        if label > 1
            labelInds(i) = true;
            stimLabels(i, 2) = label;
            if label > 97
                stimLabels(i, 1) = 0;
            else
                stimLabels(i, 1) = 1;
            end
        end
    elseif strcmp(epochInfo(i).eventtype, '1') && thirdInstSeen
        
        fourthInstSeen = true;
    elseif strcmp(epochInfo(i).eventtype, '1') && secondInstSeen
        thirdInstSeen = true;
    elseif strcmp(epochInfo(i).eventtype, '1') && firstInstSeen
        secondInstSeen = true;
    elseif strcmp(epochInfo(i).eventtype, '1') && ~firstInstSeen
        firstInstSeen = true;
    end
end

labels = stimLabels(labelInds, :);
data = EEG.data(:,:,labelInds);
data = permute(data, [3 1 2]);
data = data(:, 1:64, :);
time = EEG.times;

end

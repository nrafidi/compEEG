function alignData = pdtwEEG(oldData)
[numSamp, numE, numTime] = size(oldData);

Xs = cell(numSamp, 1);
for i = 1:numSamp
    Xs{i} = squeeze(oldData(i,:,:));
end

l = round(numTime*1.1);
ali0.P = repmat((1:l)', 1, numSamp);
aliT = [];
parDTW = struct;

ali = pdtw(Xs, ali0, aliT, parDTW);


alignData = nan(numSamp, numE, numTime);
for i = 1:numSamp
    holdData = oldData(i, :, ali.P(:,i));
    alignData(i, :, :) = holdData(:,:,1:numTime);
end

end
% Extracts features after the data have been processed

% Parameters
erpWinSize = 50;%ms
specWinSize = 200;%ms
specStartFreq = 4;%Hz, Theta
specEndFreq = 8;%Hz, Theta
subjects = 'C';
exp = 'CompEEG';
dataRoot = '/Users/nrafidi/Documents/MATLAB/compEEG-data/preproc-final/';

numSub = length(subjects);

for s = 1:numSub
    sub = subjects(s);
    
    load([dataRoot exp '_' sub '_Preprocessed.mat']);
    
    [numSamp, numChan, ~] = size(data);
    
    minTime = min(time(time >= 0));
    maxTime = max(time);
    
    erpWin = minTime:erpWinSize:maxTime;
    specWin = minTime:erpWinSize:maxTime;
    
    numErp = length(erpWin);
    numSpec = length(specWin);
    
    featData = [];
    
    for e = 2:numErp
        erp = (time >= erpWin(e-1)) & (time < erpWin(e));
        newData = squeeze(mean(data(:,:,erp), 3));
        featData = cat(2, featData, newData);
    end
    
    for p = 2:numSpec
        spec = (time >= specWin(e-1)) & (time < specWin(e));
        newData = [];
        for n = 1:numSamp
            dataToTrans = squeeze(data(n,:,spec));
            trans = abs(fft(dataToTrans'));
            freqData = mean(trans(specStartFreq:specEndFreq, :));
            newData = cat(1, newData, freqData);
        end
        featData = cat(2, featData, newData);
    end
    
    save([dataRoot exp '_' sub '_Features.mat'], 'featData', 'labels', 'erpWin', 'specWin');
    
end
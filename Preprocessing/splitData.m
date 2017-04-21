clear;
eeglab;
inputDir = '/Users/nrafidi/Documents/MATLAB/compEEG-data-rep/raw/';
outputDir = '/Users/nrafidi/Documents/MATLAB/compEEG-data-rep/preproc-partial/%s/';

subjects = {'III', 'JJJ', 'KKK'};%{'HHH'};
%{'BBB', 'GGG', 'HHH', 'AAA', 'CCC', 'DDD', 'EEE', 'FFF', 'MM', ...
 %    'NN', 'OO', 'PP', 'QQ', 'RR', 'SS', 'TT', 'UU', 'VV', 'WW',...
  %   'YY'};

numSub = length(subjects);

trueDiff255 = [154, 78, 78, 78, 78, 62, 62, 62, 62, 62, 62, 62]';
timeToAdd = 4000;
for s = 1:numSub
    sub = subjects{s}; 
    fprintf('Subject: %s\n', sub);
    
    if ismember(sub, {'CCC', 'FFF', 'OO', 'RR'})
        timeToAdd = 8000;
    elseif ismember(sub, {'DDD', 'EEE', 'PP', 'QQ', 'MM'})
        timeToAdd = 6000;
    elseif strcmp(sub, 'NN')
        timeToAdd = 10000;
    end
    
    EEG = pop_biosig([inputDir, sub, '.bdf']);
    
    subOutputDir = sprintf(outputDir, sub);
    if ~exist(subOutputDir, 'dir')
        mkdir(subOutputDir);
    end
    
    events = EEG.event;
    numEvent = length(events);
    occ255 = [];
    latency = nan(numEvent, 1);
    for iEvent = 1:numEvent
        latency(iEvent) = events(iEvent).latency;
        if events(iEvent).type == 255
            occ255 = cat(1, occ255, iEvent);
        end
    end
    
    if strcmp(sub, 'GGG')
        diffLate = diff(latency);
        occLarge = find(diffLate >= 8000);
        occ255 = [2; occLarge(1:(end-1)) + 2];
    end
    
    diff255 = diff(occ255);
    
    
    
    if length(diff255) ~= length(trueDiff255)
        fprintf('Subject %s is funky.\n', sub);
        keyboard;
        % Subjects with false starts
        if ismember(sub, {'AAA', 'SS', 'TT', 'GGG', 'HHH'})
            fprintf('False Start\n');
            
            falseStartInd = find(diff255 < 8);
            
            startFalse = (latency(occ255(falseStartInd)-2)*(EEG.times(end)/latency(end)))/1000;
            endFalse = (latency(occ255(falseStartInd+1)-2)*(EEG.times(end)/latency(end)))/1000;
            EEG = pop_select(EEG, 'notime', [startFalse, endFalse]);
            occ255 = occ255([1:(falseStartInd-1), (falseStartInd+1):end]);
            
            if strcmp(sub, 'TT')
                timeToAdd = -2000;
            elseif strcmp(sub, 'GGG')
                timeToAdd = -32000;
            elseif strcmp(sub, 'HHH')
                timeToAdd = -15000;
            end
        end
        
        %Subjects missing comp presentation
        if ismember(sub, {'MM', 'NN', 'YY'})
            occ255 = [NaN; occ255];
            keyboard;
        end
    elseif ~all(diff255 == trueDiff255)
        fprintf('Subject %s is funky.\n', sub);
        keyboard;
    end
    cutOff = ((latency(occ255(6)-2) + timeToAdd)*(EEG.times(end)/latency(end)))/1000;
    compEEG = pop_select(EEG, 'time', [0, cutOff]);
    krEEG = pop_select(EEG, 'notime', [0, cutOff]);
    
    compEEG.setname = ['CompEEG_' sub];
    krEEG.setname = ['CompEEG__KR_' sub];
    
    compEEG = eeg_checkset( compEEG );
    krEEG = eeg_checkset( krEEG );
    
    compEEG = pop_saveset( compEEG, 'filename',[compEEG.setname '.set'],'filepath',subOutputDir);
    krEEG = pop_saveset( krEEG, 'filename',[krEEG.setname '.set'],'filepath',subOutputDir);
    
    
end
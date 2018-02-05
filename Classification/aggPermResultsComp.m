%%
resultDir = '/Users/nrafidi/Documents/MATLAB/compEEG-data%s/results/';
addpath /Users/nrafidi/Documents/MATLAB/Toolboxes/export_fig/


repSubjects = {'III', 'JJJ', 'KKK', 'BBB', 'GGG', 'HHH', 'AAA', 'CCC', 'DDD', 'EEE', 'FFF', 'MM', ...
    'OO', 'PP', 'QQ', 'RR', 'SS', 'TT', 'WW',...
    'YY'};
subjects =  {'L', 'P', 'W', 'CC', 'FF', 'H', 'I', 'J', ...
    'AA', 'BB', 'DD', 'EE', 'GG', 'HH', 'JJ', 'M', 'N',...
    'F', 'K', 'O', 'R', 'S', 'T', 'U', 'V', 'X', 'Y', 'Z'};



resultDirRep = sprintf(resultDir, '-rep');

resultDir = sprintf(resultDir, '');

loadFileString_A = '%s/CompEEG_5FCV_win%d_permAccs100.mat';
loadFileString_B = 'CompEEG_%s_5FCV_win%d_permAccs100.mat';
saveFileString = 'CompEEG_5FCV_win%d_permAccs.mat';

%%
trueSubAccs = nan(length(repSubjects), 5);
permSubAccs = nan(length(repSubjects), 100, 5);
for w = [1:6, 41:54]
    for s = 1:length(repSubjects)
        sub = repSubjects{s};
        if strcmp(sub, 'OO') && w == 32
            loadFile = sprintf(loadFileString_A, sub, w);
        else
            loadFile = sprintf(loadFileString_B, sub, w);
        end
        load([resultDirRep loadFile]);
        trueSubAccs(s, :) = trueAcc;
        permSubAccs(s, :, :) = permAccs;
    end
    saveFile = sprintf(saveFileString, w);
    save([resultDirRep saveFile], 'trueSubAccs', 'permSubAccs');
end
%%
trueSubAccs = nan(length(subjects), 5);
permSubAccs = nan(length(subjects), 100, 5);
for w = [1:6, 41:54]
    for s = 1:length(subjects)
        sub = subjects{s};
        if w > 6 && w < 26
            loadFile = sprintf(loadFileString_A, sub, w);
        else
            loadFile = sprintf(loadFileString_B, sub, w);
        end
        load([resultDir loadFile]);
        trueSubAccs(s, :) = trueAcc;
        permSubAccs(s, :, :) = permAccs;
    end
    saveFile = sprintf(saveFileString, w);
    save([resultDir saveFile], 'trueSubAccs', 'permSubAccs');
end
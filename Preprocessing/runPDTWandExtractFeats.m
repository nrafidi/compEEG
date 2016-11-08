% Run PDTW and then extract features
cd /usr1/homes/nrafidi/MATLAB/groupRepo/shared/toolboxes/ctw
make
addPath
cd /usr1/homes/nrafidi/MATLAB/groupRepo/nrafidi/Preprocessing

subjects = {'L', 'P', 'W', 'CC', 'FF', 'H', 'I', 'J', ...
    'AA', 'BB', 'DD', 'F', 'EE', 'GG', 'HH', 'JJ', ...
    'K', 'M', 'O', 'R', 'S', 'T', 'U', 'V', 'X', 'Y', 'Z', 'N'};
numSub = length(subjects);

dataDir = '/usr1/homes/nrafidi/MATLAB/compEEG-data/';
fileToLoad = [dataDir, ...
    '%s/CompEEG_%s_Vis_BP2-200_N60_Ref_Epochs_Base_ICA1-2.mat'];
fileToSavePDTW = [dataDir, ...
    '%s/CompEEG_%s_Vis_BP2-200_N60_Ref_Epochs_Base_ICA1-2_PDTW.mat'];
fileToSaveFinal = [dataDir, ...
    '%s/CompEEG_%s_Vis_BP2-200_N60_Ref_Epochs_Base_ICA1-2_PDTW_Features_Overlap_Time.mat'];
h = waitbar(0, 'Subjects Completed');
for s = 1:numSub
    sub = subjects{s};
    waitbar(s/numSub, 'Subjects Completed');
    disp(sub)
    tic
    load(sprintf(fileToLoad, sub, sub));
    
    data = pdtwEEG(data);
    
    save(sprintf(fileToSavePDTW, sub, sub), 'data', 'labels', 'preProcOptions', 'time');
    
    [featData, labels, winTime, featOptions] = extractFeatures(sprintf(fileToSavePDTW, sub, sub));
    
    save(sprintf(fileToSaveFinal, sub, sub), 'featData', 'labels', 'featOptions', 'winTime');
    toc
end
subjects = {'L', 'P', 'W', 'CC', 'FF', 'H', 'I', 'J', ...
    'AA', 'BB', 'DD', 'F', 'EE', 'GG', 'HH', 'JJ', ...
    'K', 'M', 'O', 'R', 'S', 'T', 'U', 'V', 'N', 'X', 'Y', 'Z'};
numSub = length(subjects);

pOpt = {'', '_PDTW'};
numOpt = length(pOpt);

dataDir = ...
    '/Users/nrafidi/Documents/MATLAB/compEEG-data/preproc-final/';

fnameStub = [dataDir ...
    '%s/CompEEG_%s_Vis_BP2-200_N60_Ref_Epochs_Base_ICA1-2%s_Features_Overlap_Time.mat'];

for s = 1:numSub
    sub = subjects{s};
    
    for p = 1:numOpt
        load(sprintf(fnameStub, sub, sub, pOpt{p}));
        
        uniItems = unique(labels(:,2));
        numItems = length(uniItems);
        oldFeatData = featData;
        featData = nan(numItems, size(featData, 2));
        for i = 1:numItems
            featData(i,:) = mean(oldFeatData(labels(:,2) == uniItems(i), :));
        end
        
        oldLabels = labels;
        labels = [[ones(numItems/2, 1); zeros(numItems/2, 1)], uniItems];
        
        save(sprintf(fnameStub, sub, sub, [pOpt{p} '_ItemAvg']), 'featData', 'featOptions', 'labels', 'winTime');
    end
    fprintf('Subject %s Complete.\n', sub);
end
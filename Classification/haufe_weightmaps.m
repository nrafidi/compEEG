%%
fPrefixR = '/Users/nrafidi/Documents/MATLAB/compEEG-data-rep/preproc-final/';
fPrefixO = '/Users/nrafidi/Documents/MATLAB/compEEG-data/preproc-final/';
fSuffix = '_Vis_BP2-200_N60_Ref_Epochs_Base_ICA1-2_Features_Overlap_Time.mat';
subjects =  {'L', 'P', 'W', 'CC', 'FF', 'H', 'I', 'J', ...
        'AA', 'BB', 'DD', 'EE', 'GG', 'HH', 'JJ', 'M', 'N',...
        'F', 'K', 'O', 'R', 'S', 'T', 'U', 'V', 'X', 'Y', 'Z', ...
        'III', 'JJJ', 'KKK', 'BBB', 'GGG', 'HHH', 'AAA', 'CCC', 'DDD', 'EEE', 'FFF', 'MM', ...
        'OO', 'PP', 'QQ', 'RR', 'SS', 'TT', 'WW', 'YY'};
%%
load ../../compEEG-data/results/CompEEG_CV_Slide_Models_PLOS.mat
subModelsOrig = subModels;
load ../../compEEG-data-rep/results/CompEEG_CV_Slide_Models_PLOS.mat
subModels = cat(1, subModels, subModelsOrig);

[num_sub, num_win] = size(subModels);
%%
subMaps = nan(num_sub, num_win, 65);

for i_sub = 1:num_sub
    sub = subjects{i_sub};
    if i_sub < 29
        loadFname = [fPrefixO sub '/CompEEG_' sub fSuffix];
    else
        loadFname = [fPrefixR sub '/CompEEG_' sub fSuffix];
    end
    load(loadFname);
    for p = 1:num_win
        
        W = subModels{i_sub, p};
        
        startInd = (p-1)*64 + 1;
        endInd = p*64;
        X = [zscore(featData(:,startInd:endInd)), ones(size(featData, 1), 1)];
        
        S = X*W;
        
        cov_X = cov(X);
        
        cov_S = cov(S);
        
        subMaps(i_sub, p, :) = cov_X*W/cov_S;
        
    end
end

save('/Users/nrafidi/Documents/MATLAB/compEEG-data-rep/results/haufeMaps.mat', 'subMaps');

%%
load('/Users/nrafidi/Documents/MATLAB/compEEG-data-rep/results/haufeMaps.mat');
loc_file = '/Users/nrafidi/Documents/MATLAB/compEEG-data/electrode_locs.xyz';

block1 = 12:19;
block2 = 20:34;

block1Maps = squeeze(mean(mean(subMaps(:, block1, :), 1), 2));
block2Maps = squeeze(mean(mean(subMaps(:, block2, :), 1), 2));

absMax = max(abs([block1Maps; block2Maps]));


eeglab;

map = [block1Maps(1:64); zeros(8, 1)];
f = figure;
topoplot(map, loc_file, 'maplimits', [-absMax, absMax]);
colorbar
set(gca, 'fontsize', 18);
set(gcf, 'color', 'w');
export_fig(f, '../../compEEG-data-rep/results/figures/block1Topo.pdf');

map = [block2Maps(1:64); zeros(8, 1)];
f = figure;
topoplot(map, loc_file, 'maplimits', [-absMax, absMax]);
colorbar
set(gca, 'fontsize', 18);
set(gcf, 'color', 'w');
export_fig(f, '../../compEEG-data-rep/results/figures/block2Topo.pdf');
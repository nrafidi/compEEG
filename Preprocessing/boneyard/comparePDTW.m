% Compare PDTW to no PDTW (raw signal)

dataDir = '/Users/nrafidi/Documents/MATLAB/compEEG-data/preproc-final/%s/';

fileName = [dataDir 'CompEEG_%s_Vis_BP2-200_N60_Ref_Epochs_Base_ICA1-2%s.mat'];

sub = 'AA';

load(sprintf(fileName, sub, sub, ''));

comp = squeeze(mean(data(labels(:,1) == 1,:,:)));
ncomp = squeeze(mean(data(labels(:,1) == 0,:,:)));

load(sprintf(fileName, sub, sub, '_PDTW'));

comp_p = squeeze(mean(data(labels(:,1) == 1,:,:)));
ncomp_p = squeeze(mean(data(labels(:,1) == 0,:,:)));

figure
imagesc(comp)
figure
imagesc(ncomp)
figure
imagesc(comp_p)
figure
imagesc(ncomp_p)
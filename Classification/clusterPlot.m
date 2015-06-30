subjects = 'CDE';
numSub = length(subjects);

fPrefix = '/Users/nrafidi/Documents/MATLAB/compEEG-data/preproc-final/CompEEG_';
fSuffix = '_Preprocessed.mat';



for s = 1:length(subjects)
    
    sub = subjects(s);
    load([fPrefix sub fSuffix]);
    
       
    Y = labels(:,1);
    
    N = size(data, 1);
    
    figure;
    X = squeeze(mean(data, 3));
    P = size(X, 2);
    numP = sqrt(P);
    for p = 1:P
        subplot(numP, numP, p);
        clust1 = X(Y == 0, p);
        clust2 = X(Y == 1, p);
        
        m = min(length(clust1), length(clust2));
        
        hist([clust1(1:m), clust2(1:m) ]);
        if p <= P/2
            title(['A' num2str(p)]);
        else
            title(['B' num2str(p-P/2)]);
        end
    end
end
% here is a simple example showing how to train and test a Gaussian Naive Bayes classifier using this code. 
% Tom Mitchell, March 20, 2014
%

% define 
% 1. the examples (one per row, with any number of integer and real-valued features)
examples=[1 4.1 5;  
          2 5 7;
          3 3 5.2;
          8 2 6;
          9 4 9;
          10 4 8;
          8 2 6;
          9 4 9;
          10 4 8;]
% 2.  the labels (integer labels in a column vector)
labels=[3 3 3 1 1 1 1 1 1 ]'
% 3. classPriors (prior over distinct class labels, in same sequence as unique(labels) 
classPriors=[1/3 2/3]

% 4. train the GNB classifier
GNBmodel = nbayes_train(examples, labels, 1); % third arg says to use pooled variance estimate over
                                              % different class labels

% 5. apply the trained classifier to some examples (in this case, the training examples)
result = nbayes_apply(examples, GNBmodel);   
display(result);                             
% columns in result are labels sorted same as GNBmodel.labelVocab

% That's it -- enjoy!




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here's a simple k-fold train/test loop, with labels held out in the same proportion as they are
% represented in the entire set of examples (hence, class prior estimates over the training set will
% be appropriate over the test set).  You can cut and paste this into Matlab, if you have already
% defined examples, labels, and classPriors using the first part of this file.

% Should first calculate classPriors from labels
classLabels=[3,1];      % and the definition of which label is which in classPriors
nexamples=size(examples,1);
nfolds=3; % nfolds times we'll hold out 1/nfolds of the examples for test
ntest=floor(nexamples/nfolds);
ntrain=nexamples-ntest;
nperlabeltest=ntest*classPriors;  % there must be integers, of course
nperlabeltrain=ntrain*classPriors;

for k=1:length(classLabels) % create permedLabIdxs{i} as random perm of examples with ith label
    lab=classLabels(k);
    labIdxs=find(labels==lab);
    permedLabIdxs{k}=labIdxs(randperm(length(labIdxs)));
end

correct=0; incorrect=0;  % counters for test accuracy
for f=1:nfolds
    % 1. create a label-balanced set of train and test example indexes into the examples matrix
    testExampleIdxs=[]; % we'll build this up by randomly choosing the right number of 
    labelTestIdxsStarts= ((f-1)*nperlabeltest)+1;
    labelTestIdxsEnds= f*nperlabeltest;
    for k=1:size(labelTestIdxsStarts,2) % for each label, add its test examples to testExampleIdxs
       testExampleIdxs=[testExampleIdxs;  permedLabIdxs{k}(labelTestIdxsStarts(k):labelTestIdxsEnds(k))];
    end
    trainExampleIdxs=setdiff([1:nexamples],testExampleIdxs);
    
    % 2. train on training examples
    GNBmodel = nbayes_train(examples(trainExampleIdxs,:), labels(trainExampleIdxs), 1); 
    
    % 3. test on this fold's test examples
    result = nbayes_apply(examples(testExampleIdxs,:), GNBmodel);    % apply the learned  model
                                                                     % to the training examples
    [Y predictedLabelIndexes]=max(result,[],2);
    predictedLabels= GNBmodel.labelVocab(predictedLabelIndexes);
    correctLabels=labels(testExampleIdxs);
    for e=1:size(predictedLabels,1)
        if predictedLabels(e)==correctLabels(e) 
            correct=correct+1;
        else
            incorrect=incorrect+1;
        end
    end
    fprintf('true test labels: ');
    fprintf('%d ',correctLabels);
    fprintf('\npredicted test labels: ');
    fprintf('%d ', predictedLabels);
    fprintf('\n');
end

fprintf('Result of %d-fold train/test over %d examples:\n',nfolds,nexamples);
fprintf('Accuracy = %d/%d = %f\n',correct,correct+incorrect,correct/(correct+incorrect));


function [lambda] = ChooseLambda( mergeJobOutputFile )

    jobOutputStruct = load(mergeJobOutputFile, 'jobOutput');
    jobOutput = jobOutputStruct.jobOutput;
    
    Fold_X_Lambda_CVErr = zeros(size(jobOutput, 1), size(jobOutput, 2));
    Fold_X_Lambda_Lambda = zeros(size(jobOutput, 1), size(jobOutput, 2));
    
    for indexFold = 1:size(jobOutput, 1)
        for indexLambda = 1:size(jobOutput, 2)
            
            Fold_X_Lambda_CVErr(indexFold, indexLambda) = ...
                jobOutput{indexFold, indexLambda}.cvErr;
            
            Fold_X_Lambda_Lambda(indexFold, indexLambda) = ...
                jobOutput{indexFold, indexLambda}.lambda;
            
        end
    end
    
    lambda_CVErr = squeeze(mean(Fold_X_Lambda_CVErr, 1));
    [~, minIndex] = min(lambda_CVErr);
    lambda = Fold_X_Lambda_Lambda(1, minIndex);

end
function [jobOutput] = ...
    LogRegChooseLambdaAndTrainAll( ...
        dataLoader,...
        regBias, ...
        mergeJobOutputFile)
   
    lambda = ChooseLambda( mergeJobOutputFile );
    [ Features, Labels ] = dataLoader();
    Weights = GradDesc(Features, Labels, lambda, regBias);
    
    jobOutput = struct('lambda', lambda, 'Weights', Weights);
    
end
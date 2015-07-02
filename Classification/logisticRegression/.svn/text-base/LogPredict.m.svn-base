function [prediction] = LogPredict(Data, weights, isRegularizedBias)

    Input = Data;
    if (~isRegularizedBias)
        
        Input = [Input ones(size(Data, 1), 1)];
        
    end
    
    expProd = exp(Input*weights);
    sumExpProd = sum(expProd, 2);
    
    prediction = expProd ./ repmat(sumExpProd, 1, size(expProd, 2));
    depPred = ones(size(Input, 1), 1) ./ sumExpProd;
    prediction = [prediction depPred];
    
end
function [ Weights ] = StochasticAverageGradient(Obs_X_Features, Obs_X_Labels, lambda, regBias)

    eps = 1e-4;

    [numObservations, numFeatures] = size(Obs_X_Features);
    numLabelDims = size(Obs_X_Labels, 2);
    
    ApproxAverageGradient = zeros(numFeatures, numLabelDims - 1);
    LastDerivatives = zeros(numObservations, numFeatures, numLabelDims - 1);
    lastBiasDerivatives = zeros(numObservations, numLabelDims - 1);
    approxBiasAverageGradient = zeros(1, numLabelDims - 1);
    isCovered = zeros(numObservations, 1);
    Weights = zeros(numFeatures, numLabelDims - 1);
    weightsBias = zeros(1, numLabelDims - 1);
    
    % Lmax is maximum Lipschitz across observations
    Lmax = .25*max(sum(Obs_X_Features.^2,2)) + lambda;
    stepSize = 1/Lmax;
    
    tr = -1;
    iteration = 0;
    cumDeltaTr = 0;
    iterationsPerTrCheck = 20;
    
    while (true)
        
        iteration = iteration + 1;
        
        indexCurrentObservation = int32(ceil(rand(1) * numObservations));
        
        % we always pass True for isRegularizedBias to LogPredict because
        % we handle the bias weights separately in this function
        predError = Obs_X_Labels(indexCurrentObservation, :) ...
            - LogPredict(Obs_X_Features(indexCurrentObservation, :), Weights, true);
        
        tr_prev = tr;
        tr = norm(predError);
        
        if(iteration > 1)
            cumDeltaTr = cumDeltaTr + abs(tr - tr_prev);
        end
        
        if(mod(iteration,iterationsPerTrCheck) == 0)
            avgDeltaTr = cumDeltaTr / iterationsPerTrCheck;
            cumDeltaTr = 0;
            fprintf('average delta tr: %g\n', avgDeltaTr);
            if (avgDeltaTr < eps)
                if (~regBias)
                    Weights = vertcat(Weights, weightsBias); %#ok<AGROW>
                end
                break;
            end
        end
        
        derivative = zeros(numFeatures, numLabelDims - 1);
        for indexLabelDim = 1:numLabelDims - 1
            derivative(:, indexLabelDim) = Obs_X_Features(indexCurrentObservation, :) * predError(indexLabelDim);
        end
        
        derivativeDelta = derivative - squeeze(LastDerivatives(indexCurrentObservation, :, :));
        singleObs_Features_X_LabelsLess1 = repmat(Obs_X_Features(indexCurrentObservation, :)', 1, numLabelDims - 1);
        
        ApproxAverageGradient = ApproxAverageGradient + derivativeDelta .* singleObs_Features_X_LabelsLess1;
        
        LastDerivatives(indexCurrentObservation, :, :) = derivative;
        isCovered(indexCurrentObservation) = true;
        
        Weights = Weights .* (1 - stepSize * lambda);
        Weights = Weights - (stepSize/sum(isCovered)) * ApproxAverageGradient;
        
        if (~regBias)
            biasDerivative = predError(1:size(predError - 1));
            biasDerivativeDelta = biasDerivative - lastBiasDerivatives(indexCurrentObservation, :);
            approxBiasAverageGradient = approxBiasAverageGradient + biasDerivativeDelta;
            lastBiasDerivatives(indexCurrentObservation, :) = biasDerivative;
            weightsBias = weightsBias - (stepSize/sum(isCovered)) * approxBiasAverageGradient;
        end
        
    end

end
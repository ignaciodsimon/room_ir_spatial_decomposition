function [estimatedAngle, NUMBER_OF_GENERATIONS, covarianceMatrix, ...
                    POPULATION, fitnessEvaluation, initialAnglePop] = GASMLMultiSources(sensorData, ...
                                                  sensorPositions, ...
                                                  searchDistance, ...
                                                  angleRange, ...
                                                  sourceFrequency, ...
                                                  propagationSpeed, ...
                                                  sampleRate, ...
                                                  amountOfSources)
                                              
%This function will execute the a genetic algorithm to minimize the
%stochastic maximum likelihood function.

    % Constants used to define the population of angles and the number of
    % generations to be executed before stopping.
    
    POPULATION = 100;
    NUMBER_OF_GENERATIONS = 150;    
    
    % Initial values required (from input data)
    amountOfSensors = size(sensorData, 1);
    timeWindowWidth = size(sensorData, 2);
    
    % Form identity matrix
    idMatrix = eye(amountOfSensors);
        
    % Estimate the main frequency of the input signal (from the signals'
    % spectrum)
%     [~, detectedPeakFrequency] = max(abs(fft(sensorData(1,:), 24000)));
%     detectedPeakFrequency = (detectedPeakFrequency - 1) / 24000 * sampleRate;

    detectedPeakFrequency = sourceFrequency;

    % -- Compute covariance matrix
    covarianceMatrix = zeros(amountOfSensors);
    timeRange = [1, timeWindowWidth];
    for t = timeRange(1) : timeRange(2)
        covarianceMatrix = covarianceMatrix + ...
            sensorData(:, t) * sensorData(:, t)'; 
    end
    covarianceMatrix = covarianceMatrix / max(max(covarianceMatrix));



    % -- Initialize population randomly
    rng(0);
    upperAngleLimit = max(angleRange);
    lowerAngleLimit = min(angleRange);
    initialAnglePop = (upperAngleLimit-lowerAngleLimit).*rand(POPULATION,amountOfSources) + lowerAngleLimit;



    % -- Create pool of individuals from initial pool
    individualPool = initialAnglePop;

    % Form fitness evaluation matrix
    fitnessEvaluation = zeros(size(individualPool,1),1);

    % FIFO
    previous10Generations = zeros(10, amountOfSources);

    % -- Perform genetic search on solution domain
    for generation = 1 : NUMBER_OF_GENERATIONS
        for currentIndividual = 1 : size(individualPool,1)
            vandermonde = zeros(amountOfSensors, amountOfSources);
            for sourceIndex = 1 : amountOfSources
                for sensorIndex = 1 : amountOfSensors
                    currentAmplitude = 1;
                    currentAngle = individualPool(currentIndividual, sourceIndex);

                    % Position of the current sensor for the new Vandermonde input
                    currentSensorCoordinates = sensorPositions(sensorIndex, :);

                    % Of the point we're beamforming at (the radius given by the time and the angle, respect to the 0,0)
                    beamformingPosition = rms(searchDistance) * [cos(currentAngle) sin(currentAngle)];
                    currentDistanceVector = beamformingPosition - currentSensorCoordinates;

                    lambda = propagationSpeed / detectedPeakFrequency;
                    xsi = norm(currentDistanceVector) / lambda;
                    vandermonde(sensorIndex, sourceIndex) = currentAmplitude * exp(-1i * 2 * pi * xsi);% * (sensorIndex-1));
                end
            end

            % -- Form the UML cost function and calculate the fitness      

            pseudoInversVandermonde = pinv(vandermonde' * vandermonde) * vandermonde';
            projectorOnVandermonde = vandermonde * pseudoInversVandermonde;
            orthogonalProjectorOnVandermonde = idMatrix - projectorOnVandermonde;

            % -- Estimate noise power
            noisePower = 1 / (amountOfSensors - amountOfSources) * trace(orthogonalProjectorOnVandermonde * covarianceMatrix);

            % -- Estimate the signal covariance matrix
            estimatedSignalCovarianceMatrix = pseudoInversVandermonde * (covarianceMatrix - noisePower * idMatrix) * pseudoInversVandermonde';

            % -- Calculate the fitness (cost function) for given generation
            currentFitness = log( det( vandermonde * estimatedSignalCovarianceMatrix * vandermonde' + noisePower * idMatrix));

            fitnessEvaluation(currentIndividual) = real(currentFitness);
        end


        % -- Start mutation and sorting algorithms to evalutate the current
        % generation

        [sortedFitnessEvaluation, sortedIndividualPool] = rankAndSort(real(fitnessEvaluation), individualPool);

        parentPool = EMSelection(sortedIndividualPool);

        copiedPoolFull = copyPool(sortedIndividualPool);

        copiedPoolFull = bestFitNoised(copiedPoolFull, sortedIndividualPool);

        childrenPool = crossOver(parentPool);

        childrenPoolMutated = mutation(childrenPool);

        currentBestPool = pickTheBestPool(childrenPoolMutated, copiedPoolFull);

        finalPool = passBestParents(sortedIndividualPool, currentBestPool);

        individualPool = finalPool;


        % -- Save the last 10 generations and check if they have
        % changed significantly

        for fifoIndex = 1 : size(previous10Generations,1) - 1
            previous10Generations(size(previous10Generations,1) - fifoIndex + 1,:) = previous10Generations(size(previous10Generations,1) - fifoIndex,:);
        end
        previous10Generations(1,:) = sortedIndividualPool(1,:);

        if generation > 20
            if max(previous10Generations) - min(previous10Generations) < 0.00017453
                    disp(['Change is below 0.01 degree for 10 generations'])
                    disp(['Generation reached:       ' num2str(generation) ' of ' num2str(NUMBER_OF_GENERATIONS)])
                break
            end
        end
    end

    estimatedAngle = sortedIndividualPool(1,:);

end

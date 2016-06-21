function [estimatedAngles, outputPower, covarianceMatrix] = smlMultiSource(sensorData, ...
                                                  sensorPositions, ...
                                                  searchDistance, ...
                                                  scanAngleRange, ...
                                                  sourceFrequency, ...
                                                  propagationSpeed, ...
                                                  sampleRate, ...
                                                  amountOfSources)

    % This function will execute the a genetic algorithm to minimize the
    % stochastic maximum likelihood function.

    % Constants used to define the population of angles and the number of
    % generations to be executed before stopping. number of sources should
    % be estimated beforehand
    NUMBER_OF_SOURCES = 1;

    % Initial values required (from input data)
    amountOfSensors = size(sensorData, 1);
    timeWindowWidth = size(sensorData, 2);

    % Try to estimate the main frequency if no information is given about it
    if sourceFrequency(1) == 0
        [~, detectedPeakFrequency] = max(abs(fft(sensorData(1,:), sampleRate/2)));
        detectedPeakFrequency = (detectedPeakFrequency - 1) * 2;
        if detectedPeakFrequency == sampleRate/2
            detectedPeakFrequency = sampleRate / 16;
        end
        disp(sprintf('No information on frequency. Estimated: %.1f [Hz]', detectedPeakFrequency));
    else
        detectedPeakFrequency = sourceFrequency(1);
    end

    for index = 1 : size(sensorData, 1)
        sensorData(index, :) = hilbert(sensorData(index, :));
    end

    % Beamform on the time windows and obtain the estimated angle
    % -- Compute covariance matrix
    covarianceMatrix = zeros(amountOfSensors);
    timeRange = [1, timeWindowWidth];
    for t = timeRange(1) : timeRange(2)
        covarianceMatrix = covarianceMatrix + ...
            sensorData(:, t) * sensorData(:, t)'; 
    end
    covarianceMatrix = covarianceMatrix / max(max(abs(covarianceMatrix)));

    % To avoid recomputation
    identityMatrix = eye(amountOfSensors);

    % Output "power" vector
    outputPower = zeros(size(scanAngleRange));

    % -- Perform genetic search on solution domain
    for currentAngle = 1 : length(scanAngleRange)

        % Compute the steering matrix formed by Vandermonde vectors
        steeringMatrix = zeros(amountOfSensors, NUMBER_OF_SOURCES);
        for sourceIndex = 1 : NUMBER_OF_SOURCES
            for sensorIndex = 1 : amountOfSensors
                currentAmplitude = 1;
                tempAngle = scanAngleRange(currentAngle);

                % Position of the current sensor for the new Vandermonde input
                currentSensorCoordinates = sensorPositions(sensorIndex, :);

                % Of the point we're beamforming at (the radius given by the time and the angle, respect to the 0,0)
                beamformingPosition = rms(searchDistance) * [cos(tempAngle) sin(tempAngle)];
                currentDistanceVector = beamformingPosition - currentSensorCoordinates;

                lambda = propagationSpeed / detectedPeakFrequency;
                xsi = norm(currentDistanceVector) / lambda;
                steeringMatrix(sensorIndex, sourceIndex) = currentAmplitude * exp(-1i * 2 * pi * xsi);% * (sensorIndex-1));
            end
        end

        % -- Form the UML cost function and calculate the fitness
        q_hat = trace(( identityMatrix - steeringMatrix * inv(steeringMatrix' * steeringMatrix) * steeringMatrix') * covarianceMatrix)...
                / (amountOfSensors - NUMBER_OF_SOURCES);

        P_hat = inv(steeringMatrix' * steeringMatrix) * steeringMatrix' * covarianceMatrix * steeringMatrix * inv(steeringMatrix' * steeringMatrix)...
                - q_hat * inv(steeringMatrix' * steeringMatrix);

        % Save the current value of the cost function
        outputPower(currentAngle) = abs(log( det(steeringMatrix * P_hat * steeringMatrix' + q_hat * identityMatrix) ));
    end

    % Normalize the output power
    outputPower = (outputPower - min(outputPower)) / max(outputPower - min(outputPower));

    % Find needed peaks and return only those above a threshold
    [tempIRPeaksValues, tempIRPeaksPlaces] = findpeaks(outputPower, 'SortStr', 'descend','NPeaks', amountOfSources);

    estimatedAngles = [];
    for index = 1 : length(tempIRPeaksValues)
        if tempIRPeaksValues(index) > 0.2
            estimatedAngles(size(estimatedAngles, 1) + 1) = scanAngleRange(tempIRPeaksPlaces(index));
        end
    end

end

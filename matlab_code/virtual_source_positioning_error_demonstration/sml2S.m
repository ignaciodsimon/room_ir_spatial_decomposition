function [estimatedAngles, outputPower, covarianceMatrix] = sml2S(sensorData, ...
                                                  sensorPositions, ...
                                                  searchDistance, ...
                                                  scanAngleRange, ...
                                                  sourceFrequency, ...
                                                  propagationSpeed, ...
                                                  sampleRate)

    % This function will execute the a genetic algorithm to minimize the
    % stochastic maximum likelihood function.

    NUMBER_OF_SOURCES = 2;

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
    for currentAngle1Index = 1 : length(scanAngleRange)
        for currentAngle2Index = 1 : length(scanAngleRange)

            % Compute the steering matrix formed by Vandermonde vectors
            steeringMatrix = zeros(amountOfSensors, NUMBER_OF_SOURCES);
            for sensorIndex = 1 : amountOfSensors
                currentAmplitude = 1;
                tempAngle1 = scanAngleRange(currentAngle1Index);
                tempAngle2 = scanAngleRange(currentAngle2Index);

                % Position of the current sensor for the new Vandermonde input
                currentSensorCoordinates = sensorPositions(sensorIndex, :);

                % Of the point we're beamforming at (the radius given by the time and the angle, respect to the 0,0)
                beamformingPosition1 = rms(searchDistance) * [cos(tempAngle1) sin(tempAngle1)];
                beamformingPosition2 = rms(searchDistance) * [cos(tempAngle2) sin(tempAngle2)];

                currentDistanceVector1 = beamformingPosition1 - currentSensorCoordinates;
                currentDistanceVector2 = beamformingPosition2 - currentSensorCoordinates;

                lambda = propagationSpeed / detectedPeakFrequency;
                xsi1 = norm(currentDistanceVector1) / lambda;
                xsi2 = norm(currentDistanceVector2) / lambda;

                %
                steeringMatrix(sensorIndex, 1) = currentAmplitude * exp(-1i * 2 * pi * xsi1);% * (sensorIndex-1));
                steeringMatrix(sensorIndex, 2) = currentAmplitude * exp(-1i * 2 * pi * xsi2);% * (sensorIndex-1));
            end

            % -- Form the UML cost function and calculate the fitness
            q_hat = trace(( identityMatrix - steeringMatrix * inv(steeringMatrix' * steeringMatrix) * steeringMatrix') * covarianceMatrix)...
                    / (amountOfSensors - NUMBER_OF_SOURCES);

            P_hat = inv(steeringMatrix' * steeringMatrix) * steeringMatrix' * covarianceMatrix * steeringMatrix * inv(steeringMatrix' * steeringMatrix)...
                    - q_hat * inv(steeringMatrix' * steeringMatrix);

            % Save the current value of the cost function
            outputPower(currentAngle1Index, currentAngle2Index) = abs(log( det(steeringMatrix * P_hat * steeringMatrix' + q_hat * identityMatrix) ));
        end
    end

    % Normalize the output power
    outputPower = (outputPower - min(min(outputPower))) / max(max(outputPower - min(min(outputPower))));

    % Return the highest peak
    [~, estimatedAnglePosition] = max(max(outputPower));
    estimatedAngles = scanAngleRange(estimatedAnglePosition);
    disp(['Angle is: ' num2str(estimatedAngles*180/pi)])

end

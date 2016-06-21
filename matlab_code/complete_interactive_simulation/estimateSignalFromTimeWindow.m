function estimatedSignal = estimateSignalFromTimeWindow(sensorData, sensorPositions, searchDistance, ...
                                                        estimatedAngle, numberOfSources, propagationSpeed, sampleRate)
 
amountOfSensors = size(sensorData, 1);

% Estimate the main frequency of the input signal (from the signals'
% spectrum)
[~, detectedPeakFrequency] = max(abs(fft(sensorData(1,:), 24000)));
detectedPeakFrequency = (detectedPeakFrequency - 1) / 24000 * sampleRate;

% Add noise to estimated angle
% estimatedAngle = estimatedAngle + 40 * pi / 180;
% estimatedAngle = [30 * pi / 180 81 * pi / 180];
% estimatedAngle = [0 * pi / 180];


% -- Estimate signal
    for sourceIndex = 1 : numberOfSources
        for sensorIndex = 1 : amountOfSensors
            currentAmplitude = 1;
            currentAngle = estimatedAngle(sourceIndex);

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
    
    pseudoInversVandermonde = pinv(vandermonde' * vandermonde) * vandermonde';
    
    estimatedSignal = pseudoInversVandermonde * sensorData;
    
end
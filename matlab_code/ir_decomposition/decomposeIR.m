function decomposeIR(inputFilename, outputFilename, PROPAGATION_SPEED, SAMPLE_RATE, ARRAY_CENTER, originalSimulationData)

    impulseResponsesData = load(inputFilename);

%     IR_LENGTH = size(impulseResponsesData.microphonesData, 2);
%     AMOUNT_OF_SENSORS = size(impulseResponsesData.microphonesData, 1);
    INITIAL_WINDOW_SIZE = round(0.001 * SAMPLE_RATE);
%     WINDOW_SIZE_STEP = 3;
    VIRTUAL_SOURCE_POSITION_MARGIN = 0.75; % Minimum distance between sources [m]
    INITIAL_SOURCE_ALLOCATION_SIZE = 500;

    enableInteractiveMode = 0;
    showCloudEvolution = 0;
    useHammingWindow = 1;
    useConvolution = 1;
    forceSingleSource = 0;
    showDebugOutput = 0;
    displayCloudPlot = 0;

    disp(sprintf('<> IR decomposition program started.'))
    impulseResponsesData.microphonesData = impulseResponsesData.microphonesData(:, 300 : size(impulseResponsesData.microphonesData, 2));
    
    for index = 1 : size(impulseResponsesData.microphonePositions, 1)
        impulseResponsesData.microphonePositions(index,:) = impulseResponsesData.microphonePositions(index,:) - ARRAY_CENTER;
    end

    % Estimate the array center point
%     ARRAY_CENTER = [mean([max(impulseResponsesData.microphonePositions(:, 1)) min(impulseResponsesData.microphonePositions(:, 1))]) ...
%                     mean([max(impulseResponsesData.microphonePositions(:, 2)) min(impulseResponsesData.microphonePositions(:, 2))])];

    % Find the sensors separation
    sensorsDistance = norm(impulseResponsesData.microphonePositions(1,:) - impulseResponsesData.microphonePositions(2,:));
    maximumArrayFreq = PROPAGATION_SPEED / sensorsDistance /2;

    % Output data var
    foundVirtualSources = zeros(INITIAL_SOURCE_ALLOCATION_SIZE, 3);

    % Normalize IRs
    impulseResponsesData.microphonesData = impulseResponsesData.microphonesData / max(max(abs(impulseResponsesData.microphonesData)));

    % Compute a low pass average of all IRs
    filteredIRs = averageIRs(impulseResponsesData.microphonesData);

%     % Find the peaks in the average
%     [peaksValues, peaksPlaces] = findpeaks(filteredIRs);
%     peaksValues = peaksValues / max(peaksValues);
% 
%     % TODO: This has to come from some analysis of the peaks / signals ...
%     peakThreshold = 10^(-25/20);
% 
%     % Filter the peaks according to their amplitude
%     filteredPeaks = [];
%     peakVector = zeros(size(filteredIRs));
%     for index = 1 : length(peaksPlaces)
%         if peaksValues(index) > peakThreshold
%             peakVector(peaksPlaces(index)) = peaksValues(index);
%             filteredPeaks(size(filteredPeaks, 1) + 1, :) = [peaksPlaces(index) peaksValues(index)];
%         end
%     end
    
    filteredPeaks = findPeaksOnIrs(filteredIRs);
    
    
    
    

    disp(sprintf('    > Detected %d zones of interest in the IRs.', size(filteredPeaks, 1)));

    % Process detected peaks with an individual
    amountOfFoundVirtualSources = 0;
    disp(sprintf('    > Processing time window on zone %.4d of %.4d', 0, size(filteredPeaks, 1)));
    for currentPeakIndex = 1 : size(filteredPeaks, 1)

        disp(sprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b %.4d of %.4d', currentPeakIndex, size(filteredPeaks, 1)))
        currentWindowCentralTime = filteredPeaks(currentPeakIndex, 1);

        % Obtain the current window
        initPoint = round(currentWindowCentralTime - INITIAL_WINDOW_SIZE/2);
        if initPoint < 1
            initPoint = 1;
            disp('Warning: Limiting the beginning of the time window!')
        end
        endPoint = round(currentWindowCentralTime + INITIAL_WINDOW_SIZE/2);
        if endPoint > size(impulseResponsesData.microphonesData, 2);
            endPoint = size(impulseResponsesData.microphonesData, 2);
            disp('Warning: Limiting the ending of the time window!')
        end
        currentWindow = impulseResponsesData.microphonesData(:, initPoint : endPoint);

        % This needs to be tested, the interpolation increases the sampling
        % rate, so the algorithm needs to work consequently
        interpolationFactor = 1;
        interpolatedWindows = interpolateWindows(currentWindow, interpolationFactor);

        associatedDistance = currentWindowCentralTime / SAMPLE_RATE * PROPAGATION_SPEED;

        if useHammingWindow
            for index = 1 : size(interpolatedWindows, 1)
                interpolatedWindows(index, :) = interpolatedWindows(index, :) .* hamming(size(interpolatedWindows, 2))';
            end
        end

        arrayFrequency = maximumArrayFreq/2.1;
%         arrayFrequency = maximumArrayFreq/5;
        if useConvolution
            sine = sin(2 * pi * arrayFrequency / SAMPLE_RATE * [0 : 1000]);
            convolvedWindows = zeros(size(interpolatedWindows, 1), size(interpolatedWindows, 2) + length(sine) - 1);
            for index = 1 : size(interpolatedWindows, 1)
                convolvedWindows(index, :) = conv(interpolatedWindows(index, :), sine);
            end
%             interpolatedWindows = convolvedWindows;
        end


        detectedAmount = detectAmountOfSignals(interpolatedWindows);
        if showDebugOutput
            disp(sprintf('Signals detected: %d', detectedAmount));
        end

        if forceSingleSource && detectedAmount > 1
            detectedAmount = 1;
            disp('    > Forcing single source detection ...')
        end
        
        [foundAngles, outputPowers, covarianceMatrix] = smlMultiSource(convolvedWindows, ...
                                              impulseResponsesData.microphonePositions, ...
                                              associatedDistance, ...
                                              [0 : 0.005 : 2*pi], ...
                                              arrayFrequency * interpolationFactor, ...
                                              PROPAGATION_SPEED, ...
                                              SAMPLE_RATE * interpolationFactor, ...
                                              detectedAmount);

        if length(foundAngles) ~= detectedAmount
            if showDebugOutput
                disp(sprintf('Looking for %d sources, but only found %d.', detectedAmount, length(foundAngles)))
            end
            detectedAmount = length(foundAngles);
        end

        for newAngleIndex = 1 : detectedAmount

%             if ARRAY_CENTER ~= [0 0]
%                 % Finding the source position transporting the angle to the origin
%                 A = norm(ARRAY_CENTER);
%                 d = associatedDistance;
%                 beta = atan2(ARRAY_CENTER(2), ARRAY_CENTER(1));
%                 angleDiff = foundAngles(newAngleIndex) - beta;
%                 theta = asin(A * sin(angleDiff) / d);
%                 gamma = pi - angleDiff - theta;
%                 r = A * sin(gamma) / sin(theta);
%             else
%                 r = associatedDistance;
%             end
%             foundPosition = r * [cos(foundAngles(newAngleIndex)) sin(foundAngles(newAngleIndex))];
%             
            foundPosition = ARRAY_CENTER + associatedDistance * [cos(foundAngles(newAngleIndex)) sin(foundAngles(newAngleIndex))];

            % Estimate amplitude based on the distance to the main source
            if isempty(foundVirtualSources)
                currentAmplitude = 0.2;
            else
                currentAmplitude = 0.5 / (1 + 0.5*norm(foundPosition - foundVirtualSources(1, 1:2)));
            end

            if amountOfFoundVirtualSources == 0
                foundVirtualSources(1, :) = [foundPosition currentAmplitude];
                amountOfFoundVirtualSources = amountOfFoundVirtualSources + 1;
            else

                % Check that there are no existing sources within the allowed radius
                wasRepeated = 0;
                for index = 1 : amountOfFoundVirtualSources

                    if norm(foundPosition - foundVirtualSources(index, 1:2)) < VIRTUAL_SOURCE_POSITION_MARGIN

                        if showDebugOutput
                            disp('[Debug] Found source close to an existing one ...')
                        end
                        % Find the intermediate position and modify the existing one
                        foundVirtualSources(index, 1:2) = (foundPosition + foundVirtualSources(index, 1:2))/2;
                        wasRepeated = 1;
                        break;
                    end
                end

                % Just add it in case it was not repeated
                if ~wasRepeated
                    % Save the found virtual source in the total list
            %         foundVirtualSources(size(foundVirtualSources, 1) + 1, :) = [foundPosition filteredPeaks(currentPeakIndex, 2)];
                    foundVirtualSources(amountOfFoundVirtualSources + 1, :) = [foundPosition currentAmplitude];
                    amountOfFoundVirtualSources = amountOfFoundVirtualSources + 1;
                end

            end
        end


        if showCloudEvolution
            clf
            subplot(1,2,1)
            if ~isempty(originalSimulationData)
                plot(originalSimulationData.source_position(:,1), originalSimulationData.source_position(:,2), 'xg');
                hold on
            end
            plot(impulseResponsesData.sourcesPositions(:,1), impulseResponsesData.sourcesPositions(:,2), 'or');
            hold on
            plot(impulseResponsesData.microphonePositions(:, 1), impulseResponsesData.microphonePositions(:, 2), 'd');
%             plot(foundVirtualSources(1 : amountOfFoundVirtualSources, 1), foundVirtualSources(1 : amountOfFoundVirtualSources, 2), 'db');

            for index = 1 : detectedAmount
                foundPosition = r * [cos(foundAngles(index)) sin(foundAngles(index))];
                plot(foundPosition(1), foundPosition(2), 'x', 'Color', 'green', 'LineWidth', 2)
            end
            grid
            title(sprintf('Amount detected: %d', detectedAmount))
            subplot(1,2,2)
            plot([0 : 0.005 : 2*pi]/pi*180, outputPowers);
            grid
            pause
        end


        if enableInteractiveMode
            % Plotting the stuff
            clf
            displayRange = 1 : 7000;
            subplot(2,3,1)
            plot(filteredIRs(displayRange), 'LineWidth', 2)
            hold on
            plot(impulseResponsesData.microphonesData(:, displayRange)');
            plot(filteredPeaks(:,1), filteredPeaks(:,2) + 0.05, 'd', 'LineWidth', 2, 'Color', 'black');
            plot(filteredPeaks(currentPeakIndex, 1), filteredPeaks(currentPeakIndex, 2) + 0.5, 's');
            ylim([-0.5 1.5])
            grid on
            subplot(2,3,2)
            plot(interpolatedWindows')
            grid

            subplot(2,3,3)
            plot([0 : 0.005 : 2*pi]/pi*180, outputPowers);
%             hold on
%             plot([0 : 0.05 : 2*pi]/pi*180, outputPowers2);
            grid
            title(sprintf('Dist: %.1f [m]\n ML-Angle: %.1f [deg]', associatedDistance, foundAngles/pi*180));
%             title(sprintf('Dist: %.1f [m]\n ML-Angle: %.1f [deg] - Cap-Angle: %.1f [deg]', associatedDistance, foundAngle/pi*180, foundAngle2/pi*180));
            subplot(2,3,4)
            imagesc(abs(covarianceMatrix));

            subplot(2,3,5)
%             clf
            if ~isempty(originalSimulationData)
                plot(originalSimulationData.source_position(:,1), originalSimulationData.source_position(:,2), 'xg');
                hold on
            end
            plot(impulseResponsesData.sourcesPositions(:,1), impulseResponsesData.sourcesPositions(:,2), 'or');
            hold on
            plot(impulseResponsesData.microphonePositions(:, 1), impulseResponsesData.microphonePositions(:, 2), 'd');        
            plot(foundVirtualSources(1 : amountOfFoundVirtualSources, 1), foundVirtualSources(1 : amountOfFoundVirtualSources, 2), 'db');
            plot(foundPosition(1), foundPosition(2), 'x', 'Color', 'green', 'LineWidth', 2)
            grid
            title('ML')
        pause
        end
    end
    
    % Trim cloud to the amount found
    foundVirtualSources = foundVirtualSources(1 : amountOfFoundVirtualSources, :);

    % Normalize amplitudes
    foundVirtualSources(:,3) = foundVirtualSources(:,3) / max(abs(foundVirtualSources(:,3)));

    if isempty(foundVirtualSources)
        disp('No virtual sources found!')
    else
        if ~enableInteractiveMode && displayCloudPlot
            if ~isempty(originalSimulationData)
                plot(originalSimulationData.source_position(:,1), originalSimulationData.source_position(:,2), 'xg');
                hold on
            end
            plot(impulseResponsesData.sourcesPositions(:,1), impulseResponsesData.sourcesPositions(:,2), 'or');
            plot(impulseResponsesData.microphonePositions(:, 1), impulseResponsesData.microphonePositions(:, 2), 'd');        
            plot(foundVirtualSources(:, 1), foundVirtualSources(:, 2), 'db');
            grid on
        end
    end

    disp(sprintf(' > Saving results to "%s" ...', outputFilename));
    image_sources_list = foundVirtualSources;
    timeWindowSize = INITIAL_WINDOW_SIZE;
    save(outputFilename', 'image_sources_list', 'timeWindowSize');

end

function filteredIRs = averageIRs(inputData)

    filterLength = 10;
    filterTimes = 3;

    % Average all sensors inputs
    filteredIRs = abs(sum(inputData, 1));

    % Low pass filter the averaged signal
    for index = 1 : filterTimes
        filteredIRs = conv(filteredIRs, ones(filterLength, 1) / filterLength);
    end
    filteredIRs = filteredIRs(round(filterTimes * filterLength / 2) : length(filteredIRs));

end

function interpolatedWindows = interpolateWindows(windows, interpolationFactor)

%     interpolationFactor = 100;

    newXvalues = [1 : 1 / interpolationFactor : size(windows, 2)];
    interpolatedWindows = zeros(size(windows, 1), size(newXvalues, 2));
    for i = 1 : size(windows, 1)
        interpolatedWindows(i, :) = spline([1 : size(windows, 2)], windows(i, :), newXvalues);

    end

%     subplot(2,1,1)
%     plot(windows(3,:), 'x-')
%     grid
%     subplot(2,1,2)
%     plot(interpolatedWindows(3,:), 'x-')
%     grid
%     asd()

end

function detetedAmount = detectAmountOfSignals(inputSignals)

    % Compute a low pass average of all IRs
    windowPeaksCount = zeros(size(inputSignals, 1), 1);
    for index = 1 : size(inputSignals, 1)

        % Estimate the slope
        filteredIR = averageIRShortly(inputSignals(index,:));

        % Normalize and vertical stretch
        filteredIR = (filteredIR - min(filteredIR)) / max(filteredIR - min(filteredIR));

        % Find the peaks in the average
        [tempIRPeaksValues, tempIRPeaksPlaces] = findpeaks(filteredIR);
        tempIRPeaksValues = tempIRPeaksValues / max(tempIRPeaksValues);

        tempPeaksCount = 0;
        for index2 = 1 : length(tempIRPeaksPlaces)
            if tempIRPeaksValues(index2) > 0.02
                tempPeaksCount = tempPeaksCount + 1;
            end
        end
        windowPeaksCount(index) = tempPeaksCount;      
    end

    % Return the value that appears most
    detetedAmount = mode(windowPeaksCount);
end

function filteredIRs = averageIRShortly(inputSignal)

    filterLength = 5;
    filterTimes = 3;

    % Average all sensors inputs
    filteredIRs = abs(sum(inputSignal, 1));

    % Low pass filter the averaged signal
    for index = 1 : filterTimes
        filteredIRs = conv(filteredIRs, ones(filterLength, 1) / filterLength);
    end
    filteredIRs = filteredIRs(round(filterTimes * filterLength / 2) : length(filteredIRs));

end

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

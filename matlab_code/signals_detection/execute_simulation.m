function execute_simulation()

    % Load data from simulation
    impulseResponsesData = load('../array_simulation/simulated_array_irs.mat');
    originalSimulationData = load('../image_sources/simulation_results_no_diffusion_rectangular_room.mat');

    % Simulation parameters
    SAMPLE_RATE = 48000;
    IR_LENGTH = size(impulseResponsesData.microphonesData, 2);
    AMOUNT_OF_SENSORS = size(impulseResponsesData.microphonesData, 1);
    INITIAL_WINDOW_SIZE = round(0.001 * SAMPLE_RATE);
    WINDOW_SIZE_STEP = 3;
    PROPAGATION_SPEED = 343;
    VIRTUAL_SOURCE_POSITION_MARGIN = 0.50; % Minimum distance between sources [m]

    INITIAL_SOURCE_ALLOCATION_SIZE = 1000;
    showCloudEvolution = 0;
    enableInteractiveMode = 0;
    useHammingWindow = 1;
    useConvolution = 1;
    forceSingleSource = 0;


    % Steps:
    %
    % With all the IRs from the array, sliced on time windows, perform:
    %
    %   - beamforming to find the angle
    %   - estimate the associated time, based on the time window position
    %   - save the new found virtual source position
    %
    % In the end
    %
    %   - compare the found virtual source positions map with the known
    %     cloud of virtual sources
%     plot(impulseResponsesData.microphonesData')
%     grid
%     return
    impulseResponsesData.microphonesData = impulseResponsesData.microphonesData(:, 300 : size(impulseResponsesData.microphonesData, 2));

%     plot(impulseResponsesData.microphonesData')
%     grid
%     return

    % Estimate the array center point
    ARRAY_CENTER = [mean([max(impulseResponsesData.microphonePositions(:, 1)) min(impulseResponsesData.microphonePositions(:, 1))]) ...
                    mean([max(impulseResponsesData.microphonePositions(:, 2)) min(impulseResponsesData.microphonePositions(:, 2))])];

    % Find the sensors separation
    sensorsDistance = norm(impulseResponsesData.microphonePositions(1,:) - impulseResponsesData.microphonePositions(2,:));
    maximumArrayFreq = PROPAGATION_SPEED / sensorsDistance /2;

    % Output data var
    foundVirtualSources = zeros(INITIAL_SOURCE_ALLOCATION_SIZE, 3);
    foundVirtualSources2 = [];

    % Normalize IRs
    impulseResponsesData.microphonesData = impulseResponsesData.microphonesData / max(max(abs(impulseResponsesData.microphonesData)));

    % Compute a low pass average of all IRs
    filteredIRs = averageIRs(impulseResponsesData.microphonesData);

    % Find the peaks in the average
    [peaksValues, peaksPlaces] = findpeaks( filteredIRs);
    peaksValues = peaksValues / max(peaksValues);

    % TODO: This has to come from some analysis of the peaks / signals ...
    peakThreshold = 10^(-25/20);

    % Filter the peaks according to their amplitude
    filteredPeaks = [];
    peakVector = zeros(size(filteredIRs));
    for index = 1 : length(peaksPlaces)
        if peaksValues(index) > peakThreshold
            peakVector(peaksPlaces(index)) = peaksValues(index);
            filteredPeaks(size(filteredPeaks, 1) + 1, :) = [peaksPlaces(index) peaksValues(index)];
        end
    end


    % Process detected peaks with an individual
    amountOfFoundVirtualSources = 0;
    for currentPeakIndex = 1 : length(filteredPeaks)

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
        
%         for index = 1 : size(convolvedWindows, 1)
%             convolvedWindows(index, :) = convolvedWindows(index, :) .* hamming(size(convolvedWindows, 2))';
%         end


%         save('window_example.mat', 'interpolatedWindows', 'convolvedWindows')
% 
%         plot(interpolatedWindows')
%         grid
%         return

%         arrayFrequency = maximumArrayFreq /2;
%         interpolatedWindows = bandFilterWindows(interpolatedWindows, passBandCentralFreq, SAMPLE_RATE);

        % Maybe compare both, and if they're similar, take the Capon value,
        % otherwise, take the ML always

        detectedAmount = detectAmountOfSignals(interpolatedWindows);
        disp(sprintf('Signals detected: %d', detectedAmount));


%         plot(convolvedWindows')
%         grid
%         return

        if forceSingleSource && detectedAmount > 1
            detectedAmount = 1;
            disp(' > Forcing single source detection ...')
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
            disp(sprintf('Looking for %d sources, but only found %d.', detectedAmount, length(foundAngles)))
            detectedAmount = length(foundAngles);
        end


%         [foundAngles2] = GASMLDetectionIncluded(interpolatedWindows, ...
%                                               impulseResponsesData.microphonePositions, ...
%                                               associatedDistance, ...
%                                               [0 : 0.005 : 2*pi], ...
%                                               arrayFrequency * interpolationFactor, ...
%                                               PROPAGATION_SPEED, ...
%                                               SAMPLE_RATE * interpolationFactor, ...
%                                               detectedAmount);


        for newAngleIndex = 1 : detectedAmount

            % Finding the source position transporting the angle to the origin
            A = norm(ARRAY_CENTER);
            d = associatedDistance;
            beta = atan(ARRAY_CENTER(2) / ARRAY_CENTER(1));
            angleDiff = foundAngles(newAngleIndex) - beta;
            theta = asin(A * sin(angleDiff) / d);
            gamma = pi - angleDiff - theta;
            r = A * sin(gamma) / sin(theta);
            foundPosition = r * [cos(foundAngles(newAngleIndex)) sin(foundAngles(newAngleIndex))];

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

                        disp('[Debug] Found source close to an existing one ...')
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
            plot(originalSimulationData.source_position(:,1), originalSimulationData.source_position(:,2), 'xg');
            hold on
            plot(impulseResponsesData.sourcesPositions(:,1), impulseResponsesData.sourcesPositions(:,2), 'or');
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
            title(sprintf('Dist: %.1f [m]\n ML-Angle: %.1f [deg]', associatedDistance, num2str(foundAngles/pi*180)));
%             title(sprintf('Dist: %.1f [m]\n ML-Angle: %.1f [deg] - Cap-Angle: %.1f [deg]', associatedDistance, foundAngle/pi*180, foundAngle2/pi*180));
            subplot(2,3,4)
            imagesc(abs(covarianceMatrix));

            subplot(2,3,5)
%             clf
            plot(originalSimulationData.source_position(:,1), originalSimulationData.source_position(:,2), 'xg');
            hold on
            plot(impulseResponsesData.sourcesPositions(:,1), impulseResponsesData.sourcesPositions(:,2), 'or');
            plot(impulseResponsesData.microphonePositions(:, 1), impulseResponsesData.microphonePositions(:, 2), 'd');        
            plot(foundVirtualSources(1 : amountOfFoundVirtualSources, 1), foundVirtualSources(1 : amountOfFoundVirtualSources, 2), 'db');
            plot(foundPosition(1), foundPosition(2), 'x', 'Color', 'green', 'LineWidth', 2)
            grid
            title('ML')

%             subplot(2,3,6)
% %             clf
%             plot(originalSimulationData.source_position(:,1), originalSimulationData.source_position(:,2), 'xg');
%             hold on
%             plot(impulseResponsesData.sourcesPositions(:,1), impulseResponsesData.sourcesPositions(:,2), 'or');
%             plot(impulseResponsesData.microphonePositions(:, 1), impulseResponsesData.microphonePositions(:, 2), 'd');        
%             plot(foundVirtualSources2(:, 1), foundVirtualSources2(:, 2), 'db');
%             plot(foundPosition2(1), foundPosition2(2), 'x', 'Color', 'green', 'LineWidth', 2)
%             grid
%             title('Capon')

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
        if ~enableInteractiveMode
            plot(originalSimulationData.source_position(:,1), originalSimulationData.source_position(:,2), 'xg');
            hold on
            plot(impulseResponsesData.sourcesPositions(:,1), impulseResponsesData.sourcesPositions(:,2), 'or');
            plot(impulseResponsesData.microphonePositions(:, 1), impulseResponsesData.microphonePositions(:, 2), 'd');        
            plot(foundVirtualSources(:, 1), foundVirtualSources(:, 2), 'db');
            xlim([-60 60])
            grid on
        end
    end
%     return

    % TODO:
    %
    % - Implement the function to reconstruct the sound field map from the
    %   found virtual sources.
    % - Improve the SML to detect how many sources are in the window and
    %   estimate the angles of them all.
    % - Improve the positioning and timing, seems very offset.
    % - The timing might improve by estimating the group delay on each
    %   window, instead of just taking a fixed value -> This does not seem
    %   to be the problem, it's .34m after all ...
    %
    %
    %

%     % Filter virtual sources to avoid repetition
%     MINIMUM_SOURCE_DISTANCE = 0.5;
%     filteredFoundVirtualSources = [];
%     for index = 1 : size(foundVirtualSources, 1)
%         
%         if norm( ) > MINIMUM_SOURCE_DISTANCE
%         
%             % Add source directly
%             
%             % "Eliminate" source from list
%             
%         else
%             
%             % Find the intermediate position and save it
%             
%         end
%         
%     end






    disp(' > Saving results to "foundSourcesCloud.mat"');
    image_sources_list = foundVirtualSources;
    save('foundSourcesCloud_1.mat', 'image_sources_list');

    disp(' > Reconstructing the sound field power map ...');
    reconstructSoundField({'foundSourcesCloud_1.mat'}, PROPAGATION_SPEED);
    return



























%     dist = 595 / 48000 * 343;
%     
%     [foundAngle, outputPowers, covarianceMatrix] = sml(impulseResponsesData.microphonesData(:,595 : 595 + 48), ...
% %     [foundAngle, outputPowers, covarianceMatrix] = sml(impulseResponsesData.microphonesData(:,:), ...
%                                           impulseResponsesData.microphonePositions, ...
%                                           dist, ...
%                                           [0 : 0.05 : 2*pi], ...
%                                           2500, ...
%                                           PROPAGATION_SPEED, ...
%                                           SAMPLE_RATE);
% 
% 
% %     newPosition = [cos(foundAngle) sin(foundAngle)] * dist;
%     subplot(2,1,1)
%     plot(impulseResponsesData.sourcesPositions(:,1), impulseResponsesData.sourcesPositions(:,2), 'or');
%     hold on
% %     plot(originalSimulationData.source_position(:,1), originalSimulationData.source_position(:,2), 'xg');
%     plot(impulseResponsesData.microphonePositions(:, 1), impulseResponsesData.microphonePositions(:, 2), 'd');
% %     
% %     plot(newPosition(1), newPosition(2), '^');
% %     
% %     grid
% %     xlim([-1 3])
% %     ylim([0 4])
%     subplot(2,1,2)
%     polar([0 : 0.05 : 2*pi], outputPowers)
%     grid
%     return
% 
% 
% 
% 
%     plot(impulseResponsesData.microphonesData(:,595 : 595 + 48)')
%     grid
%     return
    
    
    % Estimate the noise floor to use it as threshold for peak detection
    estimatedNoiseFloor = findNoiseFloor(impulseResponsesData.microphonesData(:,:));    

    % Time window slicing
    timeWindowIndex = 1;
    timeWindowStart = 1;
    timeWindowWidth = INITIAL_WINDOW_SIZE;
    timeWindowEnd = timeWindowStart + timeWindowWidth;
    while timeWindowEnd <= IR_LENGTH
        % Modify this part to make the time window variable
        timeWindowWidth = 1 * timeWindowWidth;
        % Update value
        timeWindowEnd = timeWindowStart + timeWindowWidth;

        if timeWindowEnd > IR_LENGTH
            break;
        end
        % Current time windows
        windows = impulseResponsesData.microphonesData(:, timeWindowStart : timeWindowEnd);

%         plot(impulseResponsesData.microphonesData(1,:))
%         grid
%         return
        
        
%         if enableInteractiveMode
%             plot(windows(1,:));
%             grid on
%             ylim([-0.1 1])
%             pause(0.1)
%         end


        % Skip this time window if the energy is not enough
        if max(max(abs(windows(:, :)))) >= estimatedNoiseFloor

            
        clf
            
%             title('ok')
%             pause

            % Associated time window central time
            centralTime = (timeWindowStart + timeWindowEnd) / 2;
            associatedDistance = centralTime * PROPAGATION_SPEED / 48000;

            interpolatedWindows = interpolateWindows(windows);


            [foundAngle, outputPowers, covarianceMatrix] = sml(interpolatedWindows, ...
                                                  impulseResponsesData.microphonePositions, ...
                                                  associatedDistance, ...
                                                  [0 : 0.05 : 2*pi], ...
                                                  2000, ...
                                                  PROPAGATION_SPEED, ...
                                                  SAMPLE_RATE);

%             [foundAngle, foundAngle2, outputPowers, outputPowers2] = beamform(interpolatedWindows, impulseResponsesData.microphonePositions, associatedDistance, [0 : 0.05 : 2*pi]);
%             [foundAngle, outputPowers] = beamform(windows, impulseResponsesData.microphonePositions, associatedDistance, [0.8*foundAngle : 0.01 : 1.2*foundAngle]);













            if enableInteractiveMode

                clf
                subplot(4,2,1)
                plot(impulseResponsesData.microphonesData(1,:))
                hold on
                windowFunction = zeros(size(impulseResponsesData.microphonesData(1,:)));
                grid
                windowFunction(timeWindowStart : timeWindowEnd) = ones(1, timeWindowEnd - timeWindowStart +1);
                plot(windowFunction);
                title('Impulse responses with sliding time window')
                xlabel('Time [samples]')
                ylabel('Norm. ampl. [.]')
                subplot(4,2,3)
                startCut = timeWindowStart - 100;
                endCut = timeWindowEnd + 100;
                if startCut < 1
                    startCut = 1;
                end
                if endCut > IR_LENGTH
                    endCut = IR_LENGTH;
                end
                plot(impulseResponsesData.microphonesData(:, startCut : endCut)')
                hold on
                plot([zeros(1, 100) ones(1, timeWindowWidth) zeros(1, 100)], 'k', 'LineWidth', 2);
                grid
                title('Impulse responses with sliding time window (zoom)')
                xlabel('Time [samples]')
                ylabel('Norm. ampl. [.]')
                subplot(4,2,5)
                for kk = 1 : AMOUNT_OF_SENSORS
                    plot(interpolatedWindows(kk, :))
                    hold on
                    grid on
                end
                title('Content of time window (data interpolated)')
                xlabel('Rel. time [samples]')
                ylabel('Norm. ampl. [.]')
                subplot(4,2,7)
                plot([0 : 0.05 : 2*pi], outputPowers)
%                 hold on
%                 plot([1 : 360/length(outputPowers2): 360], outputPowers2)
                title(sprintf('Beamformers output power. Detected peak at angles: %.1f [deg]', foundAngle/pi*180));%, foundAngle2/pi*180));
                xlabel('Incoming angle [deg]')
                ylabel('Norm. ampl. [.]')
                grid on
                
                
                subplot(4,2,2)
                plot(impulseResponsesData.sourcesPositions(:,1), impulseResponsesData.sourcesPositions(:,2), 'or');
                hold on
                plot(impulseResponsesData.microphonePositions(:, 1), impulseResponsesData.microphonePositions(:, 2), 'd');
                if size(foundVirtualSourcesPositions, 1) >= 1
                    plot(foundVirtualSourcesPositions(:,1), foundVirtualSourcesPositions(:,2), 'xb');
                end
                grid
                
                
                
                pause
            end




            % Checks that there is no other found virtual source close to
            % the new one found, such that they are not added twice.
            newSourcePosition = associatedDistance * [cos(foundAngle) sin(foundAngle)];
            sourceExistsAlready = 0;
            for i = 1 : size(foundVirtualSourcesPositions, 1)
                if norm(foundVirtualSourcesPositions(i, :) - newSourcePosition) < VIRTUAL_SOURCE_POSITION_MARGIN
                    sourceExistsAlready = 1;
                    disp('repeated');
                    break
                end
            end
            if ~sourceExistsAlready

%                 imagesc(covarianceMatrix)
%                 pause

                disp(sprintf('> New virtual source found on (%.2f, %.2f)', newSourcePosition(1), newSourcePosition(2)));
                foundVirtualSourcesPositions(size(foundVirtualSourcesPositions, 1) + 1, :) = newSourcePosition;
            end

%         else
% 
%             title('pass')
%             pause

        end

        % Update values
        timeWindowIndex = timeWindowIndex + 1;
        timeWindowStart = timeWindowStart + WINDOW_SIZE_STEP;
    end

    foundVirtualSourcesPositions

    plot(foundVirtualSourcesPositions(:,1), foundVirtualSourcesPositions(:,2), 'xb');
    hold on
    plot(impulseResponsesData.sourcesPositions(:,1), impulseResponsesData.sourcesPositions(:,2), 'or');
    plot(originalSimulationData.source_position(:,1), originalSimulationData.source_position(:,2), 'xg');
    grid

end

function [estimatedAngle1, estimatedAngle2, outputPowers1, outputPowers2] = beamform(sensorData, sensorPositions, searchDistance, angleRange)

    AMOUNT_OF_SENSORS = size(sensorData, 1);
    timeWindowWidth = size(sensorData, 2);

    % Beamform on the time windows and obtain the estimated angle
    % -- Compute covariance matrix
    covarianceMatrix = zeros(AMOUNT_OF_SENSORS);
    timeRange = [1, timeWindowWidth];
    for t = timeRange(1) : timeRange(2)
        covarianceMatrix = covarianceMatrix + ...
            sensorData(:, t) * sensorData(:, t)'; 
    end
%             covarianceMatrix = covarianceMatrix / (timeRange(2) - timeRange(1) + 1);
    covarianceMatrix = covarianceMatrix / max(max(covarianceMatrix));

    % -- Scan all angles
    outputPowers1 = zeros(1, length(angleRange));
    for currentAngleIndex = 1 : length(angleRange)

        currentAngle = angleRange(currentAngleIndex);
        vandermonde = zeros(AMOUNT_OF_SENSORS, 1);
        for sensorIndex = 1 : AMOUNT_OF_SENSORS
            currentAmplitude = 1;

            % Position of the current sensor for the new Vandermonde input
            currentSensorCoordinates = sensorPositions(sensorIndex, :);

            % Of the point we're beamforming at (the radius given by the time and the angle, respect to the 0,0)
            beamformingPosition = searchDistance * [cos(currentAngle) sin(currentAngle)];
            currentDistanceVector = beamformingPosition - currentSensorCoordinates;

            lambda = 343 / 400;
            xsi = norm(currentDistanceVector) / lambda;
            vandermonde(sensorIndex) = currentAmplitude * exp(1i * 2 * pi * xsi);% * (sensorIndex-1));
        end

        % Bartlett beamformer
        outputPowers1(currentAngleIndex) = abs(ctranspose(vandermonde) * covarianceMatrix * vandermonde) / (norm(vandermonde)^2);

        % Capon beamformer
        outputPowers2(currentAngleIndex) = abs(1 / ( ctranspose(vandermonde) * inv(covarianceMatrix) * vandermonde ));

    end

    % Normalize output power
    outputPowers1 = outputPowers1 / max(abs(outputPowers1));
    outputPowers2 = outputPowers2 / max(abs(outputPowers2));
%     outputPowers2 = 1 - outputPowers2;
%     outputPowers2 = outputPowers2 / max(abs(outputPowers2));

    % Return angle corresponding to peak on the output power
    [~, estimatedAngleIndex1] = max(outputPowers1);
    estimatedAngle1 = angleRange(estimatedAngleIndex1);

    [~, estimatedAngleIndex2] = max(outputPowers2);
    estimatedAngle2 = angleRange(estimatedAngleIndex2);
    
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

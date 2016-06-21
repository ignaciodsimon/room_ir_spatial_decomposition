function execute_demonstration()

    

    % Global constants
    SAMPLE_RATE = 48000;                        % [S/s]
    SNR_VALUES = 35;%[-20 : 2 : 40];                          % [dB]
    PROPAGATION_SPEED = 343;                    % [m/s]
    AMOUNT_OF_SENSORS = [10];

    % Source characteristics
    initialSourcesDistances = [3 3];          % [m]
    sourcesInitialAngles = [+000 35];     % [deg]
    sourcesAngularSpeed = [0.2 -0.2];              % [deg/iter]
    signalDuration = 0.0003;                     % [s]
    signalsFreqs = [6000 6000];                 % [Hz]
    signalsAmplitudes = [1 1 1];                % [.]
    interpolationFactor = 1;
    initialPhase = pi/2;
    showAnimatedPlot = 0;

    % TODO: make it for several frequencies / sources at the same time
    % TODO: put a message about the maximum frequency

    % ----------------- Computations -----------------
    
    % Generates source signal with hamming window and additive white noise
    sourcesSignals = zeros(length(signalsFreqs), round(signalDuration * SAMPLE_RATE));
    for index = 1 : length(signalsFreqs)
        sourcesSignals(index, :) = (signalsAmplitudes(index) * sin(initialPhase + 2 * pi * signalsFreqs(index) / SAMPLE_RATE * [0 : size(sourcesSignals, 2)-1])) + ...
                                   (signalsAmplitudes(index) * sin(initialPhase + 2 * pi * signalsFreqs(index)*2 / SAMPLE_RATE * [0 : size(sourcesSignals, 2)-1])) + ...
                                   (signalsAmplitudes(index) * sin(initialPhase + 2 * pi * signalsFreqs(index)*3 / SAMPLE_RATE * [0 : size(sourcesSignals, 2)-1]));
    sourcesSignals(index, :) = sourcesSignals(index, :) / max(abs(sourcesSignals(index, :)));
    sourcesSignals(index, :) = sourcesSignals(index, :).^4 .* sign(sourcesSignals(index,:));
    sourcesSignals(index, :) = (sourcesSignals(index, :) .* hamming(size(sourcesSignals, 2)).^10');

    end
    sourcesSignals = interpolateWindows(sourcesSignals, interpolationFactor);
    sourcesSignals = sourcesSignals(:, 3 * interpolationFactor : size(sourcesSignals, 2));

%     plot(sourcesSignals')
%     grid
%     return

%     results = zeros(length(AMOUNT_OF_SENSORS), length(SNR_VALUES));

    currentAmountOfSensors = AMOUNT_OF_SENSORS;
    results = [];

    for currentSNRIndex = 1 : length(SNR_VALUES)
        currentSNR = SNR_VALUES(currentSNRIndex);

        disp(sprintf('>> Running for SNR: %.1f [dB]', currentSNR));

        sourcesDistances = initialSourcesDistances;
        goodCount = 0;
        for currentDistanceDiff = [0 : 0.001 : 0.20]

            if goodCount >= 5
                disp(' ** The resolution is maximum and stable. Stopping iteration.')
                break
            end

            sourcesDistances = [sourcesDistances(1) sourcesDistances(2)+currentDistanceDiff];

            % Construct the spiral array
            arrayRadius = 0.075;

            sensorsPositions = zeros(currentAmountOfSensors, 2);
            sensorsAngles = [0 : 1 / currentAmountOfSensors : 1]' * 4 * pi;
            for index = 1 : length(sensorsAngles)
                sensorsPositions(index, :) = arrayRadius * [cos(sensorsAngles(index)) sin(sensorsAngles(index)) ];
            end

            sourcesAngles = sourcesInitialAngles;
            sourcesPosition = zeros(length(signalsFreqs), 2);

            sourcesAnglesHistory = [];

            beamformerErrors = [];
            while max(sourcesAngles) <= 180
                rng(0);

                % Rotate sources and update their positions
                sourcesAngles = sourcesAngles + sourcesAngularSpeed;
                for index = 1 : length(signalsFreqs)
                    sourcesPosition(index, :) = sourcesDistances(index) * [cos(sourcesAngles(index)/180*pi) sin(sourcesAngles(index)/180*pi)];
                end
                sourcesAnglesHistory(size(sourcesAnglesHistory, 1) + 1, :) = sourcesAngles;

                % Find the largest distance between the sources and any sensor
                maxDistance = 0;
                for index1 = 1 : length(signalsFreqs)
                    for index2 = 1 : size(sensorsPositions, 1)
                        newDistance = norm(sensorsPositions(index2, :) - sourcesPosition(index1, :));
                        if newDistance > maxDistance
                            maxDistance = newDistance;
                        end
                    end
                end

                % Estimates the length needed to store the delayed signals
                receivedSignalEstimatedLength = round((signalDuration * interpolationFactor * SAMPLE_RATE) + (maxDistance / PROPAGATION_SPEED * SAMPLE_RATE));

                % Estimate the "received" signals
                sensorsSignals = zeros(size(sensorsPositions, 1), receivedSignalEstimatedLength);
                for sensorIndex = 1 : size(sensorsPositions, 1)

                    % Initialize the sensed signal vector
                    sensorsSignals(sensorIndex, :) = zeros(size(sensorsSignals(sensorIndex, :)));

                    for sourceIndex = 1 : length(signalsFreqs)

                        % Estimate the traveling signal from the source to the sensor
                        newDistance = norm(sourcesPosition(sourceIndex, :) - sensorsPositions(sensorIndex,:));
                        newDelay = newDistance / PROPAGATION_SPEED * SAMPLE_RATE;
                        newDelayedSignal = 1 / newDistance * conv(fractionalDelta(newDelay), sourcesSignals(sourceIndex, :) + (randn(1, size(sourcesSignals, 2)) * 10^(-currentSNR / 20)));

                        newDelayedSignal_padded = pad_with_zeros(newDelayedSignal, size(sensorsSignals, 2));

                        % Add the new source contribution
                        sensorsSignals(sensorIndex, :) = sensorsSignals(sensorIndex, :) + (newDelayedSignal_padded(1 : size(sensorsSignals, 2)));
                    end
                end

                sensorsSignals = sensorsSignals(:, round(0.95 * mean(sourcesDistances) / PROPAGATION_SPEED * SAMPLE_RATE) : size(sensorsSignals, 2));

                detectedAmount = detectAmountOfSignals(sensorsSignals);
% 
%                 plot(sensorsSignals(1,:)')
%                 grid
%                 pause
%                 continue
                
                
                if showAnimatedPlot
                    clf
                    subplot(2,1,1)
                    plot(sensorsPositions(:,1), sensorsPositions(:,2), 'o')
                    hold on
                    grid on
                    for indexTest = 1 : size(sourcesPosition, 1)
                        plot(sourcesPosition(indexTest, 1), sourcesPosition(indexTest, 2), 'x')
                    end
                    xlim([-10 10])
                    ylim([-4 4])

                    subplot(2,1,2)
                    for index = 1 : size(sensorsSignals, 1)
                        plot(sensorsSignals(index, :)' + 0.2*index)
                        hold on
                    end
                    grid on
                    xlim([0 size(sensorsSignals, 2)])
                    title(sprintf('Amount detected: %d', detectedAmount));
                    pause(0.01)
                end

                if (detectedAmount ~= 2) || (max(sourcesAngles) - min(sourcesAngles) < 0.5)
                    disp(sprintf('Limit for dist-diff: %.3f [m] is %.1f, %.1f [deg]', currentDistanceDiff, sourcesAngles(1), sourcesAngles(2)))
                    results(size(results, 1) + 1, :) = [currentSNR currentDistanceDiff (max(sourcesAngles) - min(sourcesAngles))];
                    
                    if  (max(sourcesAngles) - min(sourcesAngles)) < 0.5
                        goodCount = goodCount + 1;
                    else
                        goodCount = 0;
                    end

                    break
                else
                    continue
                end

            end

        end
    end
    
    save('results_data.mat', 'results')
    close all
    plot(results(:,2), results(:,3))
    grid on
    
    

end

function padded_vector = pad_with_zeros(input_vector, desired_length)
    % Returns the input vector with the desired length. It pads with zeros
    % at the end if the input is too short, or trims it if it's the other
    % way around.
    %
    % Joe.

    if length(input_vector) == desired_length
        padded_vector = input_vector;
        return
    end

    if length(input_vector) > desired_length
        padded_vector = input_vector(1 : desired_length);
    end
    
    if length(input_vector) < desired_length
        padded_vector = zeros(1, desired_length);
        padded_vector(1 : length(input_vector)) = input_vector;
    end

end

function output_ir = fractionalDelta(fractionalDelay)
    % Returns a vector with a delta on time, corresponding to the delay in
    % samples (decimal values accepted). If the delay value is not an
    % integer, a sinc function is returned.
    %
    % Joe.

    output_ir = zeros(1, round(fractionalDelay * 2));
    for i = 0 : length(output_ir)-1

        arg = i - (fractionalDelay);
        if arg ~= 0
            output_ir(i+1) = sin(pi* arg) / (pi *arg);
        else
            disp(sprintf('[Debug] The fractional-delay is actually an integer value.        \n'))
            output_ir(i+1) = 1;
        end
    end
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

function detectedAmount = detectAmountOfSignals(inputSignals)



    for index = 1 : size(inputSignals, 1)
        inputSignals(index, :) = inputSignals(index, :) / max(abs(inputSignals(index, :)));
        plot(inputSignals(index, :) + index*1.2, 'LineWidth', 1.5)
        hold on
    end
    grid on
    xlabel(sprintf('Time [samples @ f_S = %.0f kHz]', 48000/1000))
    ylabel('Microphone signals [.]')
    % Save plot to file
    set(gcf, 'PaperPosition', [-0.5 -0.2 9.4 6.05]);
    set(gcf, 'PaperSize', [8.5 5.5]);
    saveas(gcf, 'microphone_signals.pdf', 'pdf');
    close all

%     return

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
            if tempIRPeaksValues(index2) > 0.25
                tempPeaksCount = tempPeaksCount + 1;
            end
        end
        windowPeaksCount(index) = tempPeaksCount;      


        plot(filteredIR' + (index - 1)*1.5, 'LineWidth', 1.5)
        hold on
        for index3 = 1 : length(tempIRPeaksPlaces)
            if tempIRPeaksValues(index3) > 0.25
                plot(tempIRPeaksPlaces(index3), filteredIR(tempIRPeaksPlaces(index3)) + (index-1) * 1.5, 'dk')
            end
        end

    end

    % Return the value that appears most
    detectedAmount = mode(windowPeaksCount);


    xlabel(sprintf('Time [samples @ f_S = %.0f kHz]', 48000/1000))
    ylabel('Filtered microphone signals [.]')
    grid on
    % Save plot to file
    set(gcf, 'PaperPosition', [-0.5 -0.2 9.4 6.05]);
    set(gcf, 'PaperSize', [8.5 5.5]);
    saveas(gcf, 'filtered_microphone_signals.pdf', 'pdf');
    close all
    asd()
    return





end


function filteredIRs = averageIRShortly(inputSignal)

    filterLength = 3;
    filterTimes = 3;

    % Average all sensors inputs
    filteredIRs = abs(sum(inputSignal, 1));

    % Low pass filter the averaged signal
    for index = 1 : filterTimes
        filteredIRs = conv(filteredIRs, ones(filterLength, 1) / filterLength);
    end
    filteredIRs = filteredIRs(round(filterTimes * filterLength / 2) : length(filteredIRs));

end

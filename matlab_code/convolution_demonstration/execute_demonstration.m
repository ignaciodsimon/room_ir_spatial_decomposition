function execute_demonstration()

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
    VIRTUAL_SOURCE_POSITION_MARGIN = 0.75; % Minimum distance between sources [m]

    scanRange = [0 : 0.005 : 2*pi];

    INITIAL_SOURCE_ALLOCATION_SIZE = 1000;
    enableInteractiveMode = 0;
    useHammingWindow = 1;
    useConvolution = 1;

    % I know for a fact that the real source was detected on the ray
    % tracing as:
    realSourceAngle = atan(5.934/3.069);

    impulseResponsesData.microphonesData = impulseResponsesData.microphonesData(:, 300 : size(impulseResponsesData.microphonesData, 2));

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

        arrayFrequency = maximumArrayFreq/4;

        [foundAngle1, outputPowers1, covarianceMatrix_withoutConvolution] = sml(interpolatedWindows, ...
                                              impulseResponsesData.microphonePositions, ...
                                              associatedDistance, ...
                                              scanRange, ...
                                              arrayFrequency * interpolationFactor, ...
                                              PROPAGATION_SPEED, ...
                                              SAMPLE_RATE * interpolationFactor);

        % Finding the source position transporting the angle to the origin
        A = norm(ARRAY_CENTER);
        d = associatedDistance;
        beta = atan(ARRAY_CENTER(2) / ARRAY_CENTER(1));
        angleDiff = foundAngle1 - beta;
        theta = asin(A * sin(angleDiff) / d);
        gamma = pi - angleDiff - theta;
        r = A * sin(gamma) / sin(theta);
        foundPosition = r * [cos(foundAngle1) sin(foundAngle1)];
        correctedAngle1 = atan(foundPosition(2) / foundPosition(1));

        sine = sin(2 * pi * arrayFrequency / SAMPLE_RATE * [0 : 100]);
        convolvedWindows = zeros(size(interpolatedWindows, 1), size(interpolatedWindows, 2) + length(sine) - 1);
        for index = 1 : size(interpolatedWindows, 1)
            convolvedWindows(index, :) = conv(interpolatedWindows(index, :), sine);
        end

        for index = 1 : size(convolvedWindows, 1)
            convolvedWindows(index, :) = convolvedWindows(index, :) .* hamming(size(convolvedWindows, 2))';
        end

        [foundAngle2, outputPowers2, covarianceMatrix_withConvolution] = sml(convolvedWindows, ...
                                              impulseResponsesData.microphonePositions, ...
                                              associatedDistance, ...
                                              scanRange, ...
                                              arrayFrequency * interpolationFactor, ...
                                              PROPAGATION_SPEED, ...
                                              SAMPLE_RATE * interpolationFactor);

        % Finding the source position transporting the angle to the origin
        A = norm(ARRAY_CENTER);
        d = associatedDistance;
        beta = atan(ARRAY_CENTER(2) / ARRAY_CENTER(1));
        angleDiff = foundAngle2 - beta;
        theta = asin(A * sin(angleDiff) / d);
        gamma = pi - angleDiff - theta;
        r = A * sin(gamma) / sin(theta);
        foundPosition = r * [cos(foundAngle2) sin(foundAngle2)];
        correctedAngle2 = atan(foundPosition(2) / foundPosition(1));

        
        subplot(1,3,1)
        plot(interpolatedWindows')
        grid
        title('Sensed signals on time window')
        xlabel('Time [samples at 48.0 kHz]')
        ylabel('Amplitude [.]')
        subplotHandler = subplot(1,3,2);
        newPosition = get(subplotHandler, 'pos') + [-0.027 -0.025 0 0];
        set(subplotHandler, 'pos', newPosition);
        surf(abs(covarianceMatrix_withoutConvolution));
        title(sprintf('Covariance matrix estimated from sensed\nsignals convolved with a sine'))
        view([35 68])
        xlabel('Sensor number')
        ylabel('Sensor number')
        zlabel('Covariance [.]')
        subplotHandler = subplot(1,3,3);
        newPosition = get(subplotHandler, 'pos') + [-0.1 0 0 0];
        set(subplotHandler, 'pos', newPosition);
        polar(scanRange, outputPowers1)
        grid
        title(sprintf('Main peak found at angle: %.1f [deg].\nReal source at: %.1f [deg]', correctedAngle1/pi*180, realSourceAngle/pi*180));
        xlabel('Impinging angle [degrees]')
        ylabel('Amplitude [.]')
        % Save plot to file
        set(gcf, 'PaperPosition', [-2.25 +0.10 22.5 3.25]);
        set(gcf, 'PaperSize', [15 3.5]);
        saveas(gcf, 'convolution_demonstration_without.pdf', 'pdf');
        close all

        
        
        subplot(1,3,1)
        plot(convolvedWindows')
        grid
        title('Sensed signals on time window convolved with sine')
        xlabel('Time [samples at 48.0 kHz]')
        ylabel('Amplitude [.]')
        subplotHandler = subplot(1,3,2);
        newPosition = get(subplotHandler, 'pos') + [-0.027 -0.025 0 0];
        set(subplotHandler, 'pos', newPosition);
        surf(abs(covarianceMatrix_withConvolution))
        title('Covariance matrix estimated from sensed signals')
        view([35 68])
        xlabel('Sensor number')
        ylabel('Sensor number')
        zlabel('Covariance [.]')
        subplotHandler = subplot(1,3,3);
        newPosition = get(subplotHandler, 'pos') + [-0.1 0 0 0];
        set(subplotHandler, 'pos', newPosition);
        polar(scanRange, outputPowers2)
        grid
        title(sprintf('Main peak found at angle: %.1f [deg].\nReal source at: %.1f [deg]', correctedAngle2/pi*180, realSourceAngle/pi*180));
        xlabel('Impinging angle [degrees]')
        ylabel('Amplitude [.]')
        % Save plot to file
        set(gcf, 'PaperPosition', [-2.25 +0.10 22.5 3.25]);
        set(gcf, 'PaperSize', [15 3.5]);
        saveas(gcf, 'convolution_demonstration_with.pdf', 'pdf');
        close all

        return

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
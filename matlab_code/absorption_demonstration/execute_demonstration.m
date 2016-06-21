function execute_demonstration()

    PROPAGATION_SPEED = 343;
    SAMPLE_RATE = 48000;
    SOURCE_LOCATION = [1 4];
    BOUNDARY_LOCATION = [0 1 0 10];
    NOISE_LEVEL = -60;
    ARRAY_LOCATION = [2 1];
    ARRAY_RADIUS = 0.075;
    ARRAY_ELEMENTS_COUNT = 10;
    BOUNDARY_ABSORPTION = [ 125  250  500 1000 2000 4000 8000 16000 24000
%                            0.99 0.99 0.99 0.99 0.99 0.99];
                   fliplr([0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.9 0.99])];
    WINDOW_SIZE = 0.001;

    % Obtain image source location (and check that it's possible actually)
    addpath('../complete_interactive_simulation/');
    BOUNDARY_LOCATION = BOUNDARY_LOCATION + randn(size(BOUNDARY_LOCATION))/1000;
    imagePosition = findImagePosition(SOURCE_LOCATION, BOUNDARY_LOCATION);
    if ~testPointIsWithinCone(imagePosition, BOUNDARY_LOCATION, ARRAY_LOCATION)
        disp(sprintf('ERROR: the combination of boundary, source and array\n       positions do not produce a valid reflection.'))
        return
    end

    % Construct array geometry
    arrayPositions = zeros(ARRAY_ELEMENTS_COUNT, 2);
    arrayElementsAngles = [0 : 2 * pi / (ARRAY_ELEMENTS_COUNT - 1) : 2 * pi];
    for index = 1 : length(arrayElementsAngles)
        arrayPositions(index, :) = ARRAY_LOCATION + ARRAY_RADIUS * [cos(arrayElementsAngles(index)) sin(arrayElementsAngles(index))];
    end

    % Create geometry plot
    plot(arrayPositions(:,1), arrayPositions(:,2), 'o', 'MarkerSize', 3, 'Color', [.8 .4 .2])
    hold on
    plot(SOURCE_LOCATION(1), SOURCE_LOCATION(2), 'd', 'MarkerSize', 10, 'Color', [.2 .4 .8])
    plot([BOUNDARY_LOCATION(1) BOUNDARY_LOCATION(3)], [BOUNDARY_LOCATION(2) BOUNDARY_LOCATION(4)], 'LineWidth', 5)
    plot(imagePosition(1), imagePosition(2), 'd', 'MarkerSize', 10, 'Color', [.6 .7 1]);
    plot([ARRAY_LOCATION(1) SOURCE_LOCATION(1)], [ARRAY_LOCATION(2) SOURCE_LOCATION(2)], 'Color', [.4 .6 .8 .4])
    plot([SOURCE_LOCATION(1) 0 ARRAY_LOCATION(1)],[SOURCE_LOCATION(2) 3 ARRAY_LOCATION(2)], 'Color', [.8 .6 .4 .4])
    plot([0 imagePosition(1)], [3 imagePosition(2)], '--', 'Color', [.8 .6 .4 .4])
    xlim([-2.25 2.25])
    ylim([0.5 4.5])
    xlabel('X-Coordinate [m]')
    ylabel('Y-Coordinate [m]')
    grid on
    % Save plot to file
    set(gcf, 'PaperPosition', [-0.5 -0.25 9.7 7.5]);
    set(gcf, 'PaperSize', [8.5 7.0]);
    saveas(gcf, 'geometry_absorption_example.pdf', 'pdf');
    close all

    % Obtain source signal
    loadedIR = load('genelec_ir.mat');
    sourceSignal = loadedIR.ir;

    % Obtain "recordings" from direct source
    maximumRecordingLength = round(norm(imagePosition - ARRAY_LOCATION) / PROPAGATION_SPEED * SAMPLE_RATE) + length(sourceSignal);
    microphoneRecordings = zeros(ARRAY_ELEMENTS_COUNT, maximumRecordingLength);
    for currentMicrophoneIndex = 1 : ARRAY_ELEMENTS_COUNT
        % Obtain the fractional delay IR for the current source - mic pair
        currentDistance = norm(SOURCE_LOCATION - arrayPositions(currentMicrophoneIndex, :));
        currentTimeLagSamples = currentDistance / PROPAGATION_SPEED * SAMPLE_RATE;
        currentFractionalDelayIR = place_fractional_delta(currentTimeLagSamples);

        % Apply the spherical law
        currentAmplitude = 1 / currentDistance;

        % Obtain new "recording"
        newRecording = conv(sourceSignal, currentFractionalDelayIR);
        newRecording = pad_with_zeros(newRecording, maximumRecordingLength);
        microphoneRecordings(currentMicrophoneIndex, :) = newRecording * currentAmplitude;
    end
    directSignal = newRecording * currentAmplitude;

    % Obtain "recordings" from image source
    addpath('../image_sources/');
    for currentMicrophoneIndex = 1 : ARRAY_ELEMENTS_COUNT

        % Obtain the fractional delay IR for the current source - mic pair
        currentDistance = norm(imagePosition - arrayPositions(currentMicrophoneIndex, :));
        currentTimeLagSamples = currentDistance / PROPAGATION_SPEED * SAMPLE_RATE;
        currentFractionalDelayIR = place_fractional_delta(currentTimeLagSamples);

        % Apply the spherical law
        currentAmplitude = 1 / currentDistance;

        % Apply the corresponding filter to the generated delta
        delayedSignal = conv(sourceSignal, currentFractionalDelayIR);
        [filteredSignal, filterGroupDelay] = generic_filter_design(delayedSignal, ...
                                                                     BOUNDARY_ABSORPTION(1,:), ...
                                                                     BOUNDARY_ABSORPTION(2,:));

        % Correct for the filter group delay
        filteredSignal = filteredSignal(filterGroupDelay : length(filteredSignal));

        % Pad and accumulate to the existing recording
        newRecording = pad_with_zeros(filteredSignal, maximumRecordingLength);
        microphoneRecordings(currentMicrophoneIndex, :) = microphoneRecordings(currentMicrophoneIndex, :) + ...
                                                          newRecording * currentAmplitude;
    end
    reflectedSignal = newRecording * currentAmplitude;

    % Add the microphone / preamplifier noise
    rng(0);
    for index = 1 : size(microphoneRecordings, 1)
        microphoneRecordings(index, :) = microphoneRecordings(index, :) + ...
                                                          randn(size(microphoneRecordings(index, :)))*10^(NOISE_LEVEL/20);
    end

    plot([1 : 1 : length(microphoneRecordings)] / SAMPLE_RATE * 1000, microphoneRecordings')
    grid on
    xlim([0 40])
    xlabel('Time [ms]')
    ylabel('Amplitude [.]')
    % Save plot to file
    set(gcf, 'PaperPosition', [-0.5 0.1 9.7 3.0]);
    set(gcf, 'PaperSize', [8.5 3.0]);
    saveas(gcf, 'irs_absorption_example.pdf', 'pdf');
    close all
    return
    
    WINDOW_SIZE = 0.0015;

    % Extract window for direct signal
    windowWidth = round(WINDOW_SIZE * SAMPLE_RATE);
    directSignalDistance = norm(SOURCE_LOCATION - ARRAY_LOCATION);
    directSignalWindowTimeLag = directSignalDistance / PROPAGATION_SPEED * SAMPLE_RATE + (WINDOW_SIZE * SAMPLE_RATE / 5);
    directSignalStartPoint = round(directSignalWindowTimeLag - (windowWidth/2));
    directSignalEndPoint = round(directSignalWindowTimeLag + (windowWidth/2));
    directSignalWindow = microphoneRecordings(:, directSignalStartPoint : directSignalEndPoint);

    % Extract window for reflected signal
    reflectedSignalDistance = norm(imagePosition - ARRAY_LOCATION);
    reflectedSignalWindowTimeLag = reflectedSignalDistance / PROPAGATION_SPEED * SAMPLE_RATE + (WINDOW_SIZE * SAMPLE_RATE / 5);
    reflectedSignalStartPoint = round(reflectedSignalWindowTimeLag - (windowWidth/2));
    reflectedSignalEndPoint = round(reflectedSignalWindowTimeLag + (windowWidth/2));
    reflectedSignalWindow = microphoneRecordings(:, reflectedSignalStartPoint : reflectedSignalEndPoint);

    % Beamform on each window to extract the original source signal
    directSignalAngle = atan((SOURCE_LOCATION(2) - ARRAY_LOCATION(2)) / (SOURCE_LOCATION(1) - ARRAY_LOCATION(1)));
    estimatedDirectSignal = estimateSignal(directSignalWindow, arrayPositions, directSignalDistance, directSignalAngle, 1, PROPAGATION_SPEED, SAMPLE_RATE);
    estimatedDirectSignal = abs(estimatedDirectSignal) .* sign(real(estimatedDirectSignal));
    reflectedSignalAngle = atan((imagePosition(2) - ARRAY_LOCATION(2)) / (imagePosition(1) - ARRAY_LOCATION(1)));
    estimatedReflectedSignal = estimateSignal(reflectedSignalWindow, arrayPositions, reflectedSignalDistance, reflectedSignalAngle, 1, PROPAGATION_SPEED, SAMPLE_RATE);
    estimatedReflectedSignal = abs(estimatedReflectedSignal) .* sign(real(estimatedReflectedSignal));

    
    % Compare their spectrums
    directSignalSpectrum = 20*log10(abs(fft(estimatedDirectSignal)));
    directSignalSpectrum = directSignalSpectrum(1 : round(length(directSignalSpectrum)/2));
    reflectedSignalSpectrum = 20*log10(abs(fft(estimatedReflectedSignal)));
    reflectedSignalSpectrum = reflectedSignalSpectrum(1 : round(length(reflectedSignalSpectrum)/2));



    subplot(2,1,1)
    plot(directSignalWindow')
    ylim([-.5 .5])
    grid on
    title('Direct signal window')
    xlabel('Rel. time [samples @ 48 kHz]')
    ylabel('Amplitude [.]')

    subplot(2,1,2)
    plot(reflectedSignalWindow')
    ylim([-.3 .3])
    grid on
    title('Reflected signal window')
    xlabel('Rel. time [samples @ 48 kHz]')
    ylabel('Amplitude [.]')
    % Save plot to file
    set(gcf, 'PaperPosition', [-0.5 -0.25 9.7 7.5]);
    set(gcf, 'PaperSize', [8.5 7.0]);
    saveas(gcf, sprintf('time_windows_%.4f.pdf', WINDOW_SIZE), 'pdf');
    close all

%     return








    % Code to generate a plot comparing the spectrums of direct / reflected
    semilogx(nan, ':b', 'LineWidth', 3)
    hold on
    semilogx(nan, ':r', 'LineWidth', 3)
    accumulatedDirectSpectrum = [];
    for index = 1 : size(directSignalWindow, 1)
        newSpectrum = averageSpectrum(20*log10(abs(fft(directSignalWindow(index, :)))));
        
        if isempty(accumulatedDirectSpectrum)
            accumulatedDirectSpectrum = newSpectrum;
        else
            accumulatedDirectSpectrum = accumulatedDirectSpectrum + newSpectrum;
        end
        semilogx([1 : 24000/length(newSpectrum) : 24000], newSpectrum, 'LineWidth', 1.5, 'Color', [.2 .4 .6 .3])
        hold on
    end
    grid on
    accumulatedDirectSpectrum = accumulatedDirectSpectrum / size(directSignalWindow, 1);
    accumulatedReflectedSpectrum = [];
    for index = 1 : size(reflectedSignalWindow, 1)
        newSpectrum = averageSpectrum(20*log10(abs(fft(reflectedSignalWindow(index, :)))));
        if isempty(accumulatedReflectedSpectrum)
            accumulatedReflectedSpectrum = newSpectrum;
        else
            accumulatedReflectedSpectrum = accumulatedReflectedSpectrum + newSpectrum;
        end
        semilogx([1 : 24000/length(newSpectrum) : 24000], newSpectrum, 'LineWidth', 1.5, 'Color', [.6 .4 .2 .3])
        hold on
    end
    grid on
    accumulatedReflectedSpectrum = accumulatedReflectedSpectrum / size(reflectedSignalWindow, 1);
    semilogx([1 : 24000/length(newSpectrum) : 24000], accumulatedDirectSpectrum, ':b', 'LineWidth', 3)
    semilogx([1 : 24000/length(newSpectrum) : 24000], accumulatedReflectedSpectrum, ':r', 'LineWidth', 3)
    legend({'Received direct signal spectrum (average all mics)', 'Received reflected signal spectrum (average all mics)'}, 'Location', 'best')
    xlim([20 24000])
    ylim([-30 10])
    xlabel('Frequency [Hz]')
    ylabel('Power [dB]')
    % Save plot to file
    set(gcf, 'PaperPosition', [-0.5 0.05 9.7 3.05]);
    set(gcf, 'PaperSize', [8.5 3.0]);
    saveas(gcf, sprintf('received_spectrums_comparison_%.4f.pdf', WINDOW_SIZE), 'pdf');
    close all


    % SECOND PLOT - Comparison to original spectrums
    originalSpectrum = averageSpectrum(20*log10(abs(fft(directSignal))));
    imageSourceSpectrum = averageSpectrum(20*log10(abs(fft(reflectedSignal))));

    % Compare them to the estimated ones from the direct / reflected windows
    semilogx([1 : 24000/length(newSpectrum) : 24000], accumulatedDirectSpectrum, 'b', 'LineWidth', 3)
    hold on
    semilogx([1 : 24000/length(newSpectrum) : 24000], accumulatedReflectedSpectrum, 'r', 'LineWidth', 3)
    semilogx([1 : 24000/length(originalSpectrum) : 24000], originalSpectrum, ':', 'Color', [0.2 0.8 .4], 'LineWidth', 3)
    semilogx([1 : 24000/length(imageSourceSpectrum) : 24000], imageSourceSpectrum, ':k', 'LineWidth', 3)
    xlim([40 12000])
    ylim([-30 10])
    grid
    legend({'Received direct spectrum (average all mics)', 'Received reflected spectrum (average all mics)', 'Original direct spectrum', 'Original reflected spectrum'}, 'Location', 'best')
    xlabel('Frequency [Hz]')
    ylabel('Power [dB]')
    % Save plot to file
    set(gcf, 'PaperPosition', [-0.5 0.05 9.7 3.05]);
    set(gcf, 'PaperSize', [8.5 3.0]);
    saveas(gcf, sprintf('received_vs_known_spectrums_comparison_%.4f.pdf', WINDOW_SIZE), 'pdf');
    close all

    % Show the absorption curves
    
    materialCharacteristics = ones(1, length(originalSpectrum)) - 10.^((originalSpectrum - imageSourceSpectrum)/20).^(-1);
    estimatedAbsorption = ones(1, length(accumulatedDirectSpectrum)) - 10.^((accumulatedDirectSpectrum - accumulatedReflectedSpectrum)/20).^(-1);
    save(sprintf('estimated_absorption_window_%.4f.mat', WINDOW_SIZE), 'materialCharacteristics', 'estimatedAbsorption');

    semilogx([1 : 24000/length(originalSpectrum) : 24000], materialCharacteristics, 'LineWidth', 2)
    hold on
    semilogx([1 : 24000/length(newSpectrum) : 24000], estimatedAbsorption, 'LineWidth', 2)
    grid on
    xlim([63 12000])
    ylim([0.5 1])
    xlabel('Frequency [Hz]')
    ylabel('Absorption [.]')
    legend({'Material characteristics', 'Estimated absorption values'}, 'Location', 'best')
    % Save plot to file
    set(gcf, 'PaperPosition', [-0.5 0.05 9.7 3.05]);
    set(gcf, 'PaperSize', [8.5 3.0]);
    saveas(gcf, sprintf('estimated_absorption_demonstration_%.4f.pdf', WINDOW_SIZE), 'pdf');
    close all


    return

    % This is the spectrum from the beamforming (not good)
    averagedDirectSignalSpectrum = averageSpectrum(directSignalSpectrum);
    averagedReflectedSignalSpectrum = averageSpectrum(reflectedSignalSpectrum);
    semilogx([1 : 24000/length(averagedDirectSignalSpectrum) : 24000], averagedDirectSignalSpectrum, 'd-');
    semilogx([1 : 24000/length(averagedReflectedSignalSpectrum) : 24000], averagedReflectedSignalSpectrum, 's-');
    return



    for index = 1 : size(microphoneRecordings, 1)
        plot(microphoneRecordings(index, :) + index * 1.5)
        hold on
    end
    grid on


end

function padded_vector = pad_with_zeros(input_vector, desired_length)

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

function output_ir = place_fractional_delta(fractionalDelay)

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

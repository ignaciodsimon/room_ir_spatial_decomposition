function corrected_signal = correct_air_absorption(input_signal, distance)
    % Filters an input signal to simulate the absorption produced by sound
    % traveling long distances. It is based on the octave band values
    % provided by the ISO9613-2. They have been interpolated and extended
    % for a full band from 20-20000 Hz. The input data is assumed to be
    % sampled at 48000 Hz.
    %
    % Usage:
    %   corrected_signal = correct_air_absorption(input_signal, distance)
    %
    % Being the parameter "distance" in metres.
    %
    % Joe.

    SAMPLE_FREQ = 48000;
    GENERATE_INTERPOLATED_PLOT = 0;

    % Table of "atmospheric attenuation coefficient" from ISO9613-2, for a 
    % temperature of 20 degrees and 70% relative humidity.
    alpha = [  63,  0.1
              125,  0.3
              250,  1.1
              500,  2.8
             1000,  5.0
             2000,  9.0
             4000, 22.9
             8000, 76.6];

    % Code to generate a plot of interpolated values, also used to show
    % that above ~10 kHz there's no energy left (-428 dB/km at 20kHz).
    freqBins = [1 : SAMPLE_FREQ/2];
    interpolated_alpha = zeros(length(freqBins), 2);
    interpolated_alpha(:, 1) = freqBins;
    interpolated_alpha(:, 2) = spline(alpha(:,1), alpha(:,2), freqBins);

    if GENERATE_INTERPOLATED_PLOT
        % Generates the interpolated values on frequency bins following a
        % geometric progression, to have an equally distributed set on a
        % logarithmic scale
        minFreq = 20;
        maxFreq = SAMPLE_FREQ/2;
        binStep = 2^0.01;
        currentFreq = minFreq;
        counter = 0;
        while currentFreq <= maxFreq
            counter = counter + 1;
            currentFreq = currentFreq * binStep;
        end
        freqBins = zeros(1, counter);
        freqBins(1) = minFreq;
        for i = 2 : counter
            freqBins(i) = freqBins(i-1) * binStep;
        end
        interpolated_alpha = zeros(length(freqBins), 2);
        interpolated_alpha(:, 1) = freqBins;
        interpolated_alpha(:, 2) = spline(alpha(:,1), alpha(:,2), freqBins);

        % Generate a plot with the interpolated values
        figureHandler = figure();
        semilogx(alpha(:,1), alpha(:,2), 'x-', 'LineWidth', 2);
        hold on
        semilogx(interpolated_alpha(:,1), interpolated_alpha(:,2), ':', 'LineWidth', 2);
        ylim([-5 120]);
        xlim([10 20000]);
        set(gca, 'FontSize', 12);
        grid
        legend({'Values from ISO9613-2', 'Cubic-interpolation values'}, 'Location', 'best', 'FontSize', 12)
        xlabel('Frequency [Hz]', 'FontSize', 12);
        ylabel(sprintf('Atmospheric attenuation\ncoefficient [dB / km]'), 'FontSize', 12);

        % Saves plot to PDF
        set(gcf, 'PaperPosition', [-0.42 +0.05 9.3 3.6]);
        set(gcf, 'PaperSize', [8.1 3.5]);
        saveas(gcf, 'interpolated_att_coefficient.pdf', 'pdf');
        close(figureHandler);

        return
    end


% 
%     % Filter input IR with the atmospheric attenuation, for the given distance
%     filter_response = [interpolated_alpha(1,2)              % Consider the first entry also the value for DC
%                        interpolated_alpha(:,2)              % First part of the spectrum
%                        flipud(interpolated_alpha(1 : length(interpolated_alpha)-1, 2))];    % Mirrored part of the spectrum
%     filter_response = 10.^(-filter_response * distance/1000 / 20);
%     filter_ir = real(ifft(filter_response));
% 
%     % Correct the obtained IR to reduce the group delay to the minimum
%     cut_point = round(0.997 * length(filter_ir));
%     length_cut = round(0.01 * length(filter_ir));
%     filter_ir = [filter_ir(cut_point : length(filter_ir))
%                  filter_ir(1 : cut_point -1)];
%     filter_ir = filter_ir(1:length_cut);


    % Filter input IR with the atmospheric attenuation, for the given distance
    filter_module = [interpolated_alpha(1,2)                                              % Consider the first entry also the value for DC
                     interpolated_alpha(:,2)                                              % First part of the spectrum
                     flipud(interpolated_alpha(1 : length(interpolated_alpha)-1, 2))];    % Mirrored part of the spectrum
    filter_module = 10.^(-filter_module * distance/1000 / 20);

    group_delay_samples = 40;

    filter_phase = [1 : length(filter_module)/2]' / (length(filter_module)/2) * group_delay_samples * -pi;
    filter_phase = [filter_phase
                    filter_phase];

    filter_response = filter_module .* (cos(filter_phase) + 1i*sin(filter_phase));
    filter_ir = real(ifft(filter_response));
    filter_ir = filter_ir(1 : 200);
    
    % Return the input data filtered
    corrected_signal = conv(input_signal, filter_ir);
    return

end

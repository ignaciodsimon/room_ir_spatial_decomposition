function [amplitudeMap, x_coordinates, y_coordinates] = generatePowerMapFast(propagationSpeed, inputSignals_coordinates, ...
                                                                             inputSignals_amplitudes, inputSignals_freqs, ...
                                                                             x_limits, y_limits, gridResolution, saveFilename)

    if length(inputSignals_freqs) == 1
        disp(sprintf('    > Computing power map for frequency %.1f [Hz].', inputSignals_freqs));
    else
        disp(sprintf('    > Computing power map for frequencies %.1f to %.1f [Hz].', min(inputSignals_freqs), max(inputSignals_freqs)));
    end

    x_coordinates = [x_limits(1) : gridResolution : x_limits(2)];
    y_coordinates = [y_limits(1) : gridResolution : y_limits(2)];

    amplitudeMap = zeros(length(y_coordinates), length(x_coordinates));

    disp('    > Processing row 0000 of 0000');

    % For each point on the defined grid
    for x_index = 1 : length(x_coordinates)
        x_coordinate = x_coordinates(x_index);
        disp(sprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b%.4d of %.4d', x_index, length(x_coordinates)));

        for y_index = 1 : length(y_coordinates)
            y_coordinate = y_coordinates(y_index);

            % Initialize the current signal contributions
            current_point_signal = zeros(1, length(inputSignals_freqs));

            % Add the contribution of each virtual source
            for current_source_index = 1 : size(inputSignals_coordinates, 1)                
                current_source_coordinates = inputSignals_coordinates(current_source_index, :);

                % Use a complex number to represent the amplitude and phase due to the traveling
                current_distance = norm(current_source_coordinates - [x_coordinate y_coordinate]);
                newContribution = (inputSignals_amplitudes(current_source_index) / current_distance) * ...
                             exp(1i * 2 * pi .* inputSignals_freqs / propagationSpeed * current_distance);

                current_point_signal = current_point_signal + newContribution;
            end

            % Save power of current point
            amplitudeMap(y_index, x_index) = 20*log10(sum(abs(current_point_signal)));
        end
    end
    amplitudeMap = amplitudeMap - max(max(amplitudeMap));

    if ~isempty(saveFilename)
        disp(sprintf('    > Saving power map to file "%s" ...', saveFilename));
        save(saveFilename, 'amplitudeMap', 'x_coordinates', ...
            'y_coordinates', 'gridResolution', 'inputSignals_coordinates', ...
            'inputSignals_freq', 'propagationSpeed');
    end
end

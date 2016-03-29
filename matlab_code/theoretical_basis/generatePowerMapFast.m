function generatePowerMapFast(propagationSpeed, inputSignals_coordinates, ...
                              inputSignals_amplitudes, inputSignals_freq, ...
                              x_limits, y_limits, gridResolution, saveFilename)

    x_coordinates = [x_limits(1) : gridResolution : x_limits(2)];
    y_coordinates = [y_limits(1) : gridResolution : y_limits(2)];

    amplitudeMap = zeros(length(y_coordinates), length(x_coordinates));

    disp(' > Processing row 0000 of 0000');

    % For each point on the defined grid
    for x_index = 1 : length(x_coordinates)
        x_coordinate = x_coordinates(x_index);
        disp(sprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b%.4d of %.4d', x_index, length(x_coordinates)));

        for y_index = 1 : length(y_coordinates)
            y_coordinate = y_coordinates(y_index);

            % Initialize the current signal
            current_point_signal = 0;

            % Add the contribution of each virtual source
            for current_source_index = 1 : size(inputSignals_coordinates, 1)                
                current_source_coordinates = inputSignals_coordinates(current_source_index, :);

                % Use a complex number to represent the amplitude and phase due to the traveling
                current_distance = norm(current_source_coordinates - [x_coordinate y_coordinate]);
                new_signal = (inputSignals_amplitudes(current_source_index) / current_distance) * ...
                             exp(1i * 2 * pi * inputSignals_freq / propagationSpeed * current_distance);

                current_point_signal = current_point_signal + new_signal;
            end

            % Save power of current point
            amplitudeMap(y_index, x_index) = 20*log10(abs(current_point_signal));
        end
    end

    save(saveFilename, 'amplitudeMap', 'x_coordinates', ...
        'y_coordinates', 'gridResolution', 'inputSignals_coordinates', ...
        'inputSignals_freq', 'propagationSpeed');
end

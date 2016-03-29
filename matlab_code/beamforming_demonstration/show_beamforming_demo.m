function show_beamforming_demo()
    % Performs a 1D Capon beamformer simulation on an arbitrary array
    % design for different values of the noise background. It displays a
    % plot with the resulting energy as a function of the angle for each
    % noise value.
    %
    % Joe.

    % ----------------------- SIMULATION PARAMS ----------------------- 

    SOURCE_POSITION = [2, 5];
    SENSOR_NOISE_VALUES = [-20 -40 -60 -200];
    INPUT_SIGNAL_LENGTH = 1000;
    PROPAGATION_SPEED = 343;
    INPUT_SIGNAL_FREQ = 500;
    SAMPLE_RATE = 48000;
    SCAN_ANGLES = [-180 : 180];
    SENSORS_COORDINATES = [-0.20 -0.05
                          -0.15 -0.00
                          -0.10  0.05
                          -0.05  0.10
                          -0.00  0.05
                           0.05  0.00
                           0.10 -0.05
                           0.15 -0.10
                           0.20  0.10];

    % ------------------------ SIMULATION CODE ------------------------

    inputSignal = sin(2*pi*INPUT_SIGNAL_FREQ/SAMPLE_RATE * [1 : INPUT_SIGNAL_LENGTH]);
    for current_noise_level = SENSOR_NOISE_VALUES

        % Seed the random generator
        rng(0);

        % Generate the signal of the source on each sensor (complex number)
        lambda = INPUT_SIGNAL_FREQ / PROPAGATION_SPEED;
        sensorSignals = zeros(size(SENSORS_COORDINATES, 1), INPUT_SIGNAL_LENGTH);
        for i = 1 : size(SENSORS_COORDINATES, 1)
            % Obtain distance from current sensor to each source
            currentSensorCoordinates = SENSORS_COORDINATES(i, :);
            currentDistanceVector = SOURCE_POSITION - currentSensorCoordinates;

            % Update the xsi value for each angle tested
            xsi = norm(currentDistanceVector) / lambda;

            currentAmplitude = 1 / (1 + norm(currentDistanceVector));

            % "Delay" source signal as necessary
            vandermonde = currentAmplitude * exp(1i * 2 * pi * xsi);% * (i-1));
            sensorSignals(i, :) = inputSignal * vandermonde + ...
                randn(1, INPUT_SIGNAL_LENGTH) * 10^(current_noise_level/20);

            % Save the coordinates of the sensor
            SENSORS_COORDINATES(i, :) = currentSensorCoordinates;
        end

        % Generate covariance matrix of input signals
        covarianceMatrix = zeros(size(SENSORS_COORDINATES, 1));
        timeRange = [1, INPUT_SIGNAL_LENGTH];
        for t = timeRange(1) : timeRange(2)
            covarianceMatrix = covarianceMatrix + ...
                sensorSignals(:, t) * sensorSignals(:, t)';     
        end

        % Scan all angles in the defined range
        outputPower = zeros(size(SCAN_ANGLES));
        for angle_index = 1 : length(SCAN_ANGLES);
            current_angle = SCAN_ANGLES(angle_index);

            % Recreate the Vandermonde vector with associated delays
            vandermonde = zeros(size(SENSORS_COORDINATES, 1), 1);
            for sensorIndex = 1 : size(SENSORS_COORDINATES, 1)

                % Obtain distance from current sensor to each source
                currentSensorCoordinates = SENSORS_COORDINATES(sensorIndex, :);
                currentDistanceVector = 5*[cos(current_angle/180*pi) sin(current_angle/180*pi)] - currentSensorCoordinates;

                % Update the xsi value for each angle tested
                xsi = norm(currentDistanceVector) / lambda;

                % Use distance to current scan point to compensate amplitude
                currentAmplitude = (1 + norm(currentDistanceVector));

                % Obtain current entry for the Vandermonde vector
                vandermonde(sensorIndex) = currentAmplitude * exp(1i * 2 * pi * xsi);% * (sensorIndex-1));
            end

            % Bartlett beamformer
    %         outputPower(angle_index) = real((ctranspose(vandermonde) * covarianceMatrix * vandermonde) / (norm(vandermonde)^2));

            % Compute the output power for current tested angle using
            % the Capon beamformer
            outputPower(angle_index) = real(1 / ( ctranspose(vandermonde) * inv(covarianceMatrix) * vandermonde ));
        end
        outputPower = 20*log10(outputPower / max(abs(outputPower)));

        plot(SCAN_ANGLES, outputPower, 'DisplayName', sprintf('Noise level: %.0f [dBFS]', current_noise_level))
    %     polar(scan_angles/180*pi, outputPower)
        hold on
    end

    grid on;
    legend_handler = legend(gca,'show');
    set(legend_handler, 'Location', 'best');
    set(legend_handler, 'FontSize', 12);
    set(gca, 'FontSize', 12);
    xlim([min(SCAN_ANGLES) max(SCAN_ANGLES)]);
    xlabel('Angle [degrees]', 'FontSize', 12);
    ylabel('Beamformed power [dBFS]', 'FontSize', 12);
    title(sprintf('Capon beamformer output example. Frequency: %.0f[Hz]', INPUT_SIGNAL_FREQ), 'FontSize', 12);

    % Save plot to file
    set(gcf, 'PaperPosition', [-0.5 -0.35 9.4 8.05]);
    set(gcf, 'PaperSize', [8.5 7.5]);
    saveas(gcf, 'capon_beamformer_1D_example.pdf', 'pdf');
end

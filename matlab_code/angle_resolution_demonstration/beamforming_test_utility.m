function beamforming_test_utility()
    % Some notes about the geometries:
    %
    %   - The linear is obviously the worse, with a very variable
    %   resolution that depends on the impinging angle. It's also unsuable
    %   due to the total symmetry.
    %
    %   - The circular array seems to provide the same resolution as the
    %   spiral, but due to the perfect symmetry, for angles multiples of 45
    %   degrees, the resulting covariance matrix is not invertible. It can
    %   be solved adding a small error to the positioning. This way, the
    %   symmetry is not perfect and the defect does not occur. In a real
    %   case, the geometry of the array will never be 100% perfect, so this
    %   problem should not appear.
    %
    %   - The spiral array works fine in all cases. The side lobe seems to
    %   be much lower than in the case of a circular array. This is an
    %   advantage, probably when noise increases the size of the side lobes
    %   (or the covariance matrix is not very well estimated). The side
    %   lobe is most of the time half the amplitude than for the circular
    %   (or even lower). 
    %
    %   Seems that the best would be to use the circular array but
    %   inserting a bit of noise in the placement of the sensors (a small
    %   tolerance of a couple of centimeters is enough). It's the only that
    %   that has exactly the same resolution in all angles (because of the
    %   "vertical symmetry"). If the side lobe is an issue, then switch to
    %   the spiral beamformer, as it introduces a reduction of the side
    %   lobe of ~6dB.
    %
    %
    %   Time window minimum size:
    %
    %   - For the circular: seems to work free of errors mostly for windows
    %   of at least 4.0 ms.
    %
    %   - For the spiral: seems to work free of errors mostly for windows
    %   of at least 0.6 ms.
    %
    %
    %   Amount of sensors and conditions:
    %
    %   - These simulations have been conducted using 20 sensors, at a
    %   frequency of 300 Hz. The noise level was set to -25 dBFS. The
    %   source was moved in steps of 5 degrees and the beamforming
    %   algorithms were executed for steps of 1 degree.
    %
    %   In many cases, if the algorithm performs well, the precision limit
    %   is mostly due to the scan step size.
    %
    %
    %   General note:
    %
    %   - It seems that the side lobe(s) get bigger when the signal is
    %   short (few samples), so when working on critically short signals,
    %   the side lobes must be "kept under supervision", as they can easily
    %   make the algorithm give a wrong estimation. This is going to be a
    %   limitation, since we want to use short time windows if possible, so
    %   that not many reflections are contained on each window.
    %
    %
    %   Joe.


    % Global constants
    SAMPLE_RATE = 48000;                        % [S/s]
    SNR_VALUES = [0 : 5 : 85];                  % [dB]
    SNR_STEP_SIZE = 5;                          % [dB]
    PROPAGATION_SPEED = 343;                    % [m/s]
    AMOUNT_OF_SENSORS = [3 : 3 : 40];

    % Source characteristics
    sourcesDistances = [20];          % [m]
    sourcesInitialAngles = [+000];     % [deg]
    sourcesAngularSpeed = [1];              % [deg/iter]
    signalDuration = 0.001;                     % [s]
    signalsFreqs = [3500];                 % [Hz]
    signalsAmplitudes = [1 1 1];                % [.]

    % Scan setup (beamforming)
    scanAngleRange = [-180 : 1 : 180]/180*pi;

    % Array geometry selection
    useCircular = 1;
    useSpiral = 0;
    useLinear = 0;

    showInteractivePlot = 0;

    % TODO: make it for several frequencies / sources at the same time
    % TODO: put a message about the maximum frequency


    % ----------------- Computations -----------------

    % Generates source signal with hamming window and additive white noise
    sourcesSignals = zeros(length(signalsFreqs), round(signalDuration * SAMPLE_RATE));
    for index = 1 : length(signalsFreqs)
        sourcesSignals(index, :) = signalsAmplitudes(index) * sin(2 * pi * signalsFreqs(index) / SAMPLE_RATE * [1 : size(sourcesSignals, 2)]);
        sourcesSignals(index, :) = (sourcesSignals(index, :) .* hamming(size(sourcesSignals, 2))');
    end

    results = zeros(length(AMOUNT_OF_SENSORS), length(SNR_VALUES));
    
    for currentAmountOfSensorsIndex = 1 : length(AMOUNT_OF_SENSORS)
        currentAmountOfSensors = AMOUNT_OF_SENSORS(currentAmountOfSensorsIndex);
        
        disp(sprintf('*- Iteration for %d sensors:', currentAmountOfSensors))

        % Construct the spiral array
        sensorsPositions = zeros(currentAmountOfSensors, 2);
        sensorsAngles = [0 : 1 / currentAmountOfSensors : 1]' * 4 * pi;
        for index = 1 : currentAmountOfSensors
            sensorsPositions(index, :) = (0.20 + index/100) * (index - 1)/5 * [cos(sensorsAngles(index)) sin(sensorsAngles(index)) ];
        end

        for currentSNRIndex = 1 : length(SNR_VALUES)
            currentSNR = SNR_VALUES(currentSNRIndex);

            disp(sprintf('   Trying with SNR: %.1f [dB]', currentSNR));
            
            sourcesAngles = sourcesInitialAngles;
            sourcesPosition = zeros(length(signalsFreqs), 2);

            sourcesAnglesHistory = [];

            beamformerErrors = [];
            while max(sourcesAngles) < 360
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
                receivedSignalEstimatedLength = round((signalDuration * SAMPLE_RATE) + (maxDistance / PROPAGATION_SPEED * SAMPLE_RATE));

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

                sensorsSignals = sensorsSignals(:, 2500 : size(sensorsSignals, 2));

                [estimatedAngles, outputPower, covarianceMatrix] = sml(...
                                                                    sensorsSignals, ...
                                                                    sensorsPositions, ...
                                                                    sourcesDistances, ...
                                                                    scanAngleRange, ...
                                                                    signalsFreqs, ...
                                                                    PROPAGATION_SPEED, ...
                                                                    SAMPLE_RATE);

                % Errors in the angle estimations
                beamformerErrors(size(beamformerErrors, 1) + 1, :) = estimatedAngles/pi*180 - sourcesAngles;

            end
            % Count how many peaks are within the +-0.5 degrees margin
            amountOfErrorsBelowHalfDegree = 0;
            for index = 1 : length(beamformerErrors)
                if abs(wrapTo180(beamformerErrors(index))) <= 0.5
                    amountOfErrorsBelowHalfDegree = amountOfErrorsBelowHalfDegree + 1;
                end
            end
            amountOfErrorsBelowHalfDegree = amountOfErrorsBelowHalfDegree / length(beamformerErrors) * 100;

            disp(sprintf('     Percentage of fit solutions: %.1f', amountOfErrorsBelowHalfDegree));



            % Save new value in the results table
            results(currentAmountOfSensorsIndex, currentSNRIndex) = amountOfErrorsBelowHalfDegree;
            disp('** Saving results ...');
            save('results.mat', 'results', 'AMOUNT_OF_SENSORS', 'SNR_VALUES');
        end

    end

    disp('** Saving results ...');
    save('results.mat', 'results', 'AMOUNT_OF_SENSORS', 'SNR_VALUES');
    disp('** All done ***');
    

% 
%     subplot(2,1,1)
%     plot(sensorsSignals')
%     grid on
%     subplot(2,1,2)
%     plot(sourcesAnglesHistory, wrapTo180(beamformerErrors), 'LineWidth', 2);
%     ylim([-5 5])
%     grid on
% 
%     impingingAngles = sourcesAnglesHistory;
%     amountOfSensors = AMOUNT_OF_SENSORS;
%     filename = sprintf('angle_error_%d_sensors_noise_level_%.1f.mat', amountOfSensors, NOISE_LEVELS);
%     save(filename, 'impingingAngles', 'beamformerErrors', 'amountOfSensors');




















return
    
    
    
    
    

% 
%     for sourceAngleIndex = [1 : length(sourceAngles)]
%         disp(sprintf('Step: %d of %d', sourceAngleIndex, length(sourceAngles)))
% 
%         % Update the source position for the current angle
%         realAngle = sourceAngles(sourceAngleIndex);
%         sourcesPosition = sourcesDistances * [cos(realAngle/180*pi) sin(realAngle/180*pi)];
% 
%         % Finds the largest distance between the source and any sensor
%         maxDistance = sqrt(max(sum((sensorsPositions - repmat(sourcesPosition, size(sensorsPositions, 1), 1))'.^2)));
% 
%         % Estimates the length needed to store the delayed signals
%         receivedSignalEstimatedLength = round((signalDuration * SAMPLE_RATE) + (maxDistance / PROPAGATION_SPEED * SAMPLE_RATE));
% 
%         % Variable for estimated received signals
%         sensorsSignals = zeros(size(sensorsPositions, 1), receivedSignalEstimatedLength);
% 
%         % Estimate the received signals
%         for index = 1 : size(sensorsPositions, 1)
%             newDistance = norm(sourcesPosition - sensorsPositions(index,:));
%             newDelay = newDistance / PROPAGATION_SPEED * SAMPLE_RATE;
%             newDelayedSignal = 1 / newDistance * conv(fractionalDelta(newDelay), sourcesSignals + (randn(size(sourcesSignals)) * 10^(NOISE_LEVEL / 20)));
%             sensorsSignals(index, :) = (newDelayedSignal(1 : size(sensorsSignals, 2)));
%         end
% 
%         
% 
%     plot(sensorsSignals')
%     grid
%     return
% 
% 
% 
% 
% 
%         % Beamforming execution
%         searchDistance = norm(sourcePosition);
%         [estimatedBartlett, estimatedCapon, powerBartlett, powerCapon, covarianceMatrix] = beamformers(sensorsSignals, sensorsPositions, searchDistance, scanAngleRange, signalsFreqs, PROPAGATION_SPEED);
% 
%         % Errors in the estimations
%         newErrorBartlett = (realAngle/180*pi) - estimatedBartlett;
%         newErrorCapon = (realAngle/180*pi) - estimatedCapon;
% %         newErrorBartlett = wrapToPi(newErrorBartlett);
% %         newErrorCapon = wrapToPi(newErrorCapon);
%         bartlettErrors(sourceAngleIndex) = newErrorBartlett;
%         caponErrors(sourceAngleIndex) = newErrorCapon;
% %         bartlettSNR(sourceAngleIndex) = 
% %         caponSNR(sourceAngleIndex) = 
% 
%         if showInteractivePlot
%             clf
%             subplot(2,2,2)
%             plot(sensorsPositions(:,1), sensorsPositions(:, 2), 'x');
%             hold on
%             plot(sourcePosition(1), sourcePosition(2), 'o')
%             grid
%             title(sprintf('Sensors (x) and source (o)\nReal angle: %.1f [deg]', realAngle))
%             xlabel('X-axis [m]')
%             ylabel('Y-axis [m]')
%             xlim([-sourcesDistances*1.2 sourcesDistances*1.2])
%             ylim([-sourcesDistances*1.2 sourcesDistances*1.2])
% 
%             subplot(2,2,1)
%             plot(scanAngleRange/pi*180, powerBartlett)
%             grid
%             title(sprintf('Bartlett beamformer power. Peak at: %.1f [deg]', estimatedBartlett/pi*180));
%             xlabel('Angle [deg]')
%             ylabel('Power [.]')
% 
%             subplot(2,2,3)
%             plot(scanAngleRange/pi*180, powerCapon)
%             grid
%             title(sprintf('Capon beamformer power. Peak at: %.1f [deg]', estimatedCapon/pi*180));
%             xlabel('Angle [deg]')
%             ylabel('Power [.]')
% 
%             subplot(2,2,4)
%             imagesc(abs(covarianceMatrix))
%             title('Estimated covariance matrix')
% %             plot(sourceSignal)
% %             xlim([1 length(sourceSignal)])
% %             title(sprintf('Source signal. Freq: %.1f [Hz]', signalFreq))
% %             xlabel('Time [samples]')
% %             ylabel('Amplitude [.]')
% %             grid
%             pause(0.5)
%         end
% 
%     end
% 
%     figure
%     subplot(2,1,1)
%     plot(sourceAngles, bartlettErrors/pi*180)
%     title('Bartlett beamformer estimated angle error')
%     xlabel('Impinging angle [deg]')
%     ylabel('Error [deg]')
% %     ylim([-10 10])
%     grid
%     subplot(2,1,2)
%     plot(sourceAngles, caponErrors/pi*180)
%     title('Capon beamformer estimated angle  error')
%     xlabel('Impinging angle [deg]')
%     ylabel('Error [deg]')
%     grid
%     return

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

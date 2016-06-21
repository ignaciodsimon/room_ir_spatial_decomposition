function execute_array_simulation(varargin)
    % Uses the cloud of virtual sources obtained with the ray-tracing
    % simulation to estimate the IRs on the points defined by the array
    % geometry.
    %
    % Joe.

    if nargin > 0
        MAX_VIRTUAL_SOURCE_ORDER = cell2mat(varargin(3));
        SAMPLE_RATE = cell2mat(varargin(4));
        PROPAGATION_SPEED = cell2mat(varargin(5));
        AMOUNT_OF_SENSORS = cell2mat(varargin(6));
        ARRAY_ORIGIN = cell2mat(varargin(7));
        ARRAY_RADIUS = cell2mat(varargin(8));
    else
        disp('[WARNING] Array simulation is being executed with default parameters since non were given!')
        SAMPLE_RATE = 48000;
        PROPAGATION_SPEED = 343;
        MAX_VIRTUAL_SOURCE_ORDER = 10;
        AMOUNT_OF_SENSORS = 25;
        ARRAY_ORIGIN = [1 1];
        ARRAY_RADIUS = 0.075;
    end

    addNoiseToPositioning = 0;

%     % Positions of the sensors in the array
%     microphonePositions = [0.50 0.5
%                            0.51 0.5
%                            0.52 0.5
%                            0.53 0.5
%                            0.54 0.5];

% ARRAY_ORIGIN = [0 0];

    % Create circular array
    microphonePositions = zeros(AMOUNT_OF_SENSORS, 2);
    for index = 1 : AMOUNT_OF_SENSORS
        tempAngle = index / AMOUNT_OF_SENSORS * 2 * pi;
        microphonePositions(index, :) = ARRAY_ORIGIN + (ARRAY_RADIUS * [cos(tempAngle) sin(tempAngle)]);
    end
    
    % Adds noise to the positions
    if addNoiseToPositioning
        microphonePositions = microphonePositions + randn(size(microphonePositions))*10^(-40/20);
    end

%     % Plot array elements positions
%     for i = 1 : size(microphonePositions)
%         plot(microphonePositions(i,1), microphonePositions(i,2), 'xb');
%         hold on
%     end
%     grid on
% %     xlim([-0.5 2])
% %     xlim([-0.5 2])
%     return




    % Load data from simulation (position and impulse response of every
    % virtual source)
    if nargin == 0
        inputFilename = '../image_sources/simulation_results_no_diffusion_rectangular_room.mat';
    else
        inputFilename = char(varargin(1));
    end
    disp(sprintf('<> IR synthesizer started ...\n    Loading data from file "%s" ...', inputFilename));
    loadedData = load(inputFilename);

%     for index = 1 : size(loadedData.image_sources_list(:,1:4))
%         disp(index)
%         disp(loadedData.image_sources_list(index,1:4))
%     end
%     return
%     loadedData.image_sources_list = [loadedData.image_sources_list(10, :)
%                                      loadedData.image_sources_list(17, :)];
%     return
%     loadedData.image_sources_list = loadedData.image_sources_list([99 86 49], :);


    % Output data
    microphonesData = zeros(size(microphonePositions, 1), size(loadedData.image_sources_list, 2)-4);

    % Compute for each microphone position and find the total IR
    disp('    > Computing IR for microphone position: 000');
    for micIndex = 1 : size(microphonePositions, 1)
        disp(sprintf('\b\b\b\b\b %.3d', micIndex));
        for sourceIndex = 1 : size(loadedData.image_sources_list, 1)

            if loadedData.image_sources_list(sourceIndex, 4) <= MAX_VIRTUAL_SOURCE_ORDER

                % New component
                newComponent = loadedData.image_sources_list(sourceIndex, 5 : size(loadedData.image_sources_list, 2));

                % Delay new component
                distance = sqrt( (loadedData.image_sources_list(sourceIndex, 1) - microphonePositions(micIndex, 1))^2 + ...
                                 (loadedData.image_sources_list(sourceIndex, 2) - microphonePositions(micIndex, 2))^2 );
                delayValue = distance / PROPAGATION_SPEED * SAMPLE_RATE;
                newComponent = conv(place_fractional_delta(delayValue), newComponent);
                newComponent = newComponent(1 : size(loadedData.image_sources_list, 2) - 4);
                newComponent = newComponent / distance;

                % Add it to the total IR
                microphonesData(micIndex, :) = microphonesData(micIndex, :) + newComponent;
                
            end
        end
    end

    % Count the amount of sources that fulfil the maximum order criteria
    sourcesPositions = [];
    for sourceIndex = 1 : size(loadedData.image_sources_list, 1)
        if loadedData.image_sources_list(sourceIndex, 4) <= MAX_VIRTUAL_SOURCE_ORDER
            sourcesPositions(size(sourcesPositions, 1) + 1, :) = loadedData.image_sources_list(sourceIndex, 1 : 2);
        end
    end

    
%     % Plot obtained IRs
%     plot(microphonesData(1, :))
%     hold on
%     plot(microphonesData(2, :))
%     plot(microphonesData(3, :))
%     plot(microphonesData(4, :))
%     plot(microphonesData(5, :))
%     grid
%     pause
%     figure
%     plot(microphonesData')
%     title(sprintf('%d ', ARRAY_ORIGIN));
%     grid on
%     pause
%     figure
%     for ii = 1 : size(microphonesData, 1)
%         plot(microphonesData(ii,:) + 0.2*ii)
%         hold on
%     end
%     grid on
%     pause
    

    % Save results to file
    disp('    > Saving synthesised IRs to file "simulated_array_irs.mat" ...');
    
    if nargin > 0
        outputFilename = char(varargin(2));
    else
        outputFilename = 'simulated_array_irs.mat';
    end
    save(outputFilename, 'microphonesData', 'microphonePositions', 'sourcesPositions');
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

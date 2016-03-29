function execute_simulation()

    % Simulation constants
    PROPAGATION_SPEED = 343;
%     TEST_FREQUENCIES = [107 200 249 333 503 999];
%     TEST_FREQUENCIES = [60 : 10 : 3000];
    TEST_FREQUENCIES = [500];
    X_LIMITS = [0 7];
    Y_LIMITS = [0 8];
    GRID_RESOLUTION = 0.0125;
    MAXIMUM_VIRTUAL_SOURCE_ORDER = 10;

    % Load list of virtual sources with location, order and amplitude. It
    % also contains the boundaries of the room used when simulating the
    % virtual sources
%     loaded_data = load('room_simulation_results.mat');
    loaded_data = load('simulation_results_no_diffusion_order_10.mat');


    % Filter the sources to a maximum order and make sure there's only one
    % direct ray
    filteredSourceList = [];
    foundDirectSource = 0;
    keepSingleDirectSource = 1;
    counter = 1;
    for i = 1 : size(loaded_data.image_sources_list, 1)
        if loaded_data.image_sources_list(i, 4) <= MAXIMUM_VIRTUAL_SOURCE_ORDER

            % It's not first order
            if loaded_data.image_sources_list(i, 4) ~= 0
                add = 1;
            else
                % It's first order
                if ~foundDirectSource
                    add = 1;
                    if keepSingleDirectSource
                        foundDirectSource = 1;
                    end
                else
                    add = 0;
                end
            end

            if add
                filteredSourceList(counter, :) = loaded_data.image_sources_list(i, :);
                counter = counter + 1;
            end
        end
    end

    for currentFrequency_index = 1 : length(TEST_FREQUENCIES)
        currentFrequency = TEST_FREQUENCIES(currentFrequency_index);
        disp(sprintf(' > Frequency: %d [Hz] (iteration: %d/%d)', TEST_FREQUENCIES(currentFrequency_index), currentFrequency_index, length(TEST_FREQUENCIES)));
        saveFilename = sprintf('power_map_maxorder_%d_freq_%d.mat', MAXIMUM_VIRTUAL_SOURCE_ORDER, currentFrequency);

        % Execute simulation
        inputSignals_coordinates = filteredSourceList(:, 1 : 2);
        inputSignals_amplitudes = filteredSourceList(:, 3);
        generatePowerMapFast(PROPAGATION_SPEED, inputSignals_coordinates, ...
                                  inputSignals_amplitudes, currentFrequency, ...
                                  X_LIMITS, Y_LIMITS, GRID_RESOLUTION, saveFilename);
    end
    
    saveFilenames = {};
    for i = 1 : length(TEST_FREQUENCIES)
        saveFilenames{i} = sprintf('power_map_maxorder_%d_freq_%d.mat', MAXIMUM_VIRTUAL_SOURCE_ORDER, TEST_FREQUENCIES(i));
    end

    % Plot results from simulation
    room_boundaries = [0 0 7 0
                       7 0 7 8
                       7 8 0 8
                       0 8 0 0];
    room_boundaries = [];
    plotAmplitudeMap_image_fast(saveFilenames, room_boundaries, ...
                                TEST_FREQUENCIES, PROPAGATION_SPEED, MAXIMUM_VIRTUAL_SOURCE_ORDER);

end

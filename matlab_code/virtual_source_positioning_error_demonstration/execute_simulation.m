function execute_simulation()
    % This script is a copy of the one used for the theoretical
    % demonstrations but includes the posibility of adding gaussian noise
    % to the positions of the virtual sources

    % Simulation constants
    PROPAGATION_SPEED = 343;
%     TEST_FREQUENCIES = [107 200 249 333 503 999];
%     TEST_FREQUENCIES = [60 : 10 : 3000];
    TEST_FREQUENCIES = [150 : 7 : 4000];
    X_LIMITS = [-0.5 7.5];%[6 7.2];
    Y_LIMITS = [-0.5 8.5];%[-0.5 8.5];
    GRID_RESOLUTION = 0.05;%0.0125;
    MAXIMUM_VIRTUAL_SOURCE_ORDER = 5;
    noise_level = -12.0412;

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

    % Adds noise to the position of the sources
    rng(0);
    filteredSourceList(:, 1:2) = filteredSourceList(:, 1:2) + (10^(noise_level/20) * randn(size(filteredSourceList(:, 1:2))));

    for currentFrequency_index = 1 : length(TEST_FREQUENCIES)
        currentFrequency = TEST_FREQUENCIES(currentFrequency_index);
        disp(sprintf(' > Frequency: %d [Hz] (iteration: %d/%d)', TEST_FREQUENCIES(currentFrequency_index), currentFrequency_index, length(TEST_FREQUENCIES)));
        saveFilename = sprintf('power_map_maxorder_%d_freq_%d_noise_%.3f.mat', MAXIMUM_VIRTUAL_SOURCE_ORDER, currentFrequency, 10^(noise_level/20));

        % Execute simulation
        inputSignals_coordinates = filteredSourceList(:, 1 : 2);
        inputSignals_amplitudes = filteredSourceList(:, 3);
        generatePowerMapFast(PROPAGATION_SPEED, inputSignals_coordinates, ...
                                  inputSignals_amplitudes, currentFrequency, ...
                                  X_LIMITS, Y_LIMITS, GRID_RESOLUTION, saveFilename);
    end

    saveFilenames = {};
    for i = 1 : length(TEST_FREQUENCIES)
        saveFilenames{i} = sprintf('power_map_maxorder_%d_freq_%d_noise_%.3f.mat', MAXIMUM_VIRTUAL_SOURCE_ORDER, TEST_FREQUENCIES(i), 10^(noise_level/20));
    end

    % Plot results from simulation
    room_boundaries = [0 0.25 7 0.25
                       7 0.25 7 8
                       7 8 0 8
                       0 8 0 0.25];
%     room_boundaries = [];
    plotAmplitudeMap_image_fast(saveFilenames, room_boundaries, ...
                                TEST_FREQUENCIES, PROPAGATION_SPEED, MAXIMUM_VIRTUAL_SOURCE_ORDER, 10^(noise_level/20));

end

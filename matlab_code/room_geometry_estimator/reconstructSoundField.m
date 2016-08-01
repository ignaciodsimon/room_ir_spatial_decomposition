function reconstructSoundField(foundVirtualSourcesFilenames, propagationSpeed, powerMapFrequencies, outputFilename, x_limits, y_limits, gridResolution)
    % Reconstructs the sound field power map from a file with the cloud of
    % estimated virtual sources.
    %
    % Joe

%     originalSimulationData = load();


    % TODO: This should be calculated from the sources list, taking into acount
    % the ones that are closest to the main source.
%     x_limits = [-1 7];
%     y_limits = [-1 9];
%     gridResolution = 0.05;

%     primeNumbers = primes(3000);

%     powerMapFrequencies = primeNumbers(length(primeNumbers)-350 : length(primeNumbers));%[200 : 7 : 3000];
%     powerMapFrequencies = 2000;
    MAXIMUM_VIRTUAL_SOURCE_ORDER = -1;

    disp(sprintf('<> Sound field reconstruction utility started.'))

    newAmplitudeMap = [];
    for filenameIndex = 1 : length(foundVirtualSourcesFilenames)

        loadedData = load(char(foundVirtualSourcesFilenames(filenameIndex)));
%         foundVirtualSources = loadedData.foundVirtualSources;
            foundVirtualSources = loadedData.image_sources_list;

        inputSignals_coordinates = foundVirtualSources(:,1:2);
        inputSignals_freq = powerMapFrequencies;

        % Generate the new power map
        [amplitudeMap, x_coordinates, y_coordinates] = ...
            generatePowerMapFast(propagationSpeed, inputSignals_coordinates, ...
                                 foundVirtualSources(:,3), powerMapFrequencies, ...
                                 x_limits, y_limits, gridResolution, '');

        if isempty(newAmplitudeMap)
            newAmplitudeMap = amplitudeMap;
        else
            newAmplitudeMap = newAmplitudeMap + amplitudeMap;
        end
    end
    amplitudeMap = newAmplitudeMap;

    save(outputFilename, 'amplitudeMap', 'x_coordinates', ...
            'y_coordinates', 'gridResolution', 'inputSignals_coordinates', ...
            'inputSignals_freq', 'propagationSpeed');

    return
    % Plot results from simulation
    room_boundaries = [0 0 6 0
                       6 0 6 8
                       6 8 0 8
                       0 8 0 0];
    room_boundaries = [];
    plotAmplitudeMap_image_fast({outputFilename}, room_boundaries, ...
                                powerMapFrequencies, propagationSpeed, MAXIMUM_VIRTUAL_SOURCE_ORDER);













    return

    amplitudeMap = [];
    for currentFrequency = [200 : 50 : 3000];

        [newAmplitudeMap, x_coordinates, y_coordinates] = ...
            generatePowerMapFast(propagationSpeed, foundVirtualSources(:, 1 : 2), ...
                                 foundVirtualSources(:, 3), currentFrequency, ...
                                 x_limits, y_limits, gridResolution, '');
%         [newAmplitudeMap, x_coordinates, y_coordinates] = ...
%             generatePowerMapFast(PROPAGATION_SPEED, originalSimulationData.image_sources_list(:,1:2), ...
%                                  originalSimulationData.image_sources_list(:,1:2), currentFrequency, ...
%                                  x_limits, y_limits, gridResolution, '');

        if isempty(amplitudeMap)
            amplitudeMap = newAmplitudeMap;
        else
            amplitudeMap = amplitudeMap + newAmplitudeMap;
        end
    end

    % Normalize map
    amplitudeMap = amplitudeMap / max(max(abs(amplitudeMap)));

    disp('Saving map ...')
    save('total_map_found.mat', 'amplitudeMap', 'x_coordinates', 'y_coordinates');

%     disp('Displaying map')
%     figure
%     imagesc(x_coordinates, y_coordinates, 20*log10(abs(amplitudeMap)))

end

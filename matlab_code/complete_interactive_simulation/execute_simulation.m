function execute_simulation()
    clc
    disp(sprintf('------------------------------------------\n *** GENERAL INTERACTIVE DEMO STARTED ***\n   Jose I. D. Simon - Christian S. Jøns\n  AAU - Acoustics - Master Thesis - 2016\n------------------------------------------\n\n'))

    roomGeometryFilename = 'room_geometry_definition.mat';

    % Execute the room geometry capture utility
    disp(' > Launching geometry capture utility ...')
    import_room_geometry(roomGeometryFilename);

    disp(sprintf('\b\b [DONE]'))
    disp(' > Loading captured geometry ...')
    loadedGeometry = load(roomGeometryFilename);
    disp(sprintf('\b\b [DONE]'))

    % Create global simulation parameters
    ROOM_BOUNDARIES = loadedGeometry.room_boundaries;
    ARRAY_CENTER = loadedGeometry.mic_position;
    SOURCE_POSITIONS = loadedGeometry.source_positions;
    BOUNDARY_COLLISION_ERROR_MARGIN = 0.05;
    INCLUDE_REAL_SOURCE = 1;
    MAX_REFLECTION_ORDER = 5;
    MIC_COLLISION_ERROR_MARGIN = 0.05;
    PROPAGATION_SPEED = 343;
    SAMPLE_RATE = 48000;
    SCAN_ANGLE_MAX = 180;
    SCAN_ANGLE_MIN = -180;
    SCAN_ANGLE_RES = 0.25;
    AMOUNT_OF_SENSORS = 25;
    ARRAY_RADIUS = 0.075;
    POWER_MAP_FREQUENCIES = primes(2000);
    x_limits = [-4 10];
    y_limits = [-4 10];
    POWER_MAP_GRID_RESOLUTION = 0.05;

    % Temporal override
%     ROOM_BOUNDARIES = [-1 -1 -1 10
%                        -1 10 6 10
%                        6 10 6 -1
%                        6 -1 -1 -1];
%     ROOM_BOUNDARIES = [-1 -1 -2 10
%                        -2 10 6 9
%                        6 9 5.5 -1.5
%                        5.5 -1.5 -1 -1];
%     ARRAY_CENTER = [0 0];
%     SOURCE_POSITIONS = [2.5 4.5];
%                         4 6];
%                         [2 3
%                         2 8
%                         5 7
%                         4 1];
%     x_limits = [-2 7];
%     y_limits = [-2 11];
    POWER_MAP_GRID_RESOLUTION = 0.1;



%     ARRAY_CENTER = [5 9];
    show_animated_plot = 0;
    animated_plot_pause_duration = 0.1;

    recomputeEverything = 1;


    % Generate synthetic materials with the same absorption for all frequencies
    disp(' > Generating synthetic materials ...')
    ABSORPTION_FREQUENCIES = [125,250,500,1000,2000,4000,8000,16000,24000];
    ABSORPTION_COEFFICIENTS = ones(length(loadedGeometry.boundaries_reflection_coefficients), length(ABSORPTION_FREQUENCIES));
    for index = 1 : size(ABSORPTION_COEFFICIENTS, 1)
        ABSORPTION_COEFFICIENTS(index, :) = ABSORPTION_COEFFICIENTS(index, :) * loadedGeometry.boundaries_reflection_coefficients(index);
    end
    disp(sprintf('\b\b [DONE]'))

    if recomputeEverything
        disp(sprintf('\n[>] Main computation part started.\n'))
        for currentSourcePositionIndex = 1 : size(SOURCE_POSITIONS, 1)

            % Update the new source position
            currentSourcePosition = SOURCE_POSITIONS(currentSourcePositionIndex, :);
            disp(sprintf(' > Source is now placed at (%.1f, %.1f) [m]. Position %d of %d.', currentSourcePosition(1), currentSourcePosition(2), currentSourcePositionIndex, size(SOURCE_POSITIONS, 1)))

            % Save current simulation vars to file (for ray-tracing)
            disp(' > Generating temp file for current simulation ...')
            tempRaytracingFilename = sprintf('ray_tracing_%d.mat', currentSourcePositionIndex);
            generateTempFileForRayTracing(ROOM_BOUNDARIES, ARRAY_CENTER, currentSourcePosition, BOUNDARY_COLLISION_ERROR_MARGIN, ...
                                          INCLUDE_REAL_SOURCE, MAX_REFLECTION_ORDER, MIC_COLLISION_ERROR_MARGIN, PROPAGATION_SPEED, ...
                                          SAMPLE_RATE, SCAN_ANGLE_MAX, SCAN_ANGLE_MIN, SCAN_ANGLE_RES, ABSORPTION_COEFFICIENTS, ABSORPTION_FREQUENCIES, ...
                                          tempRaytracingFilename, show_animated_plot, animated_plot_pause_duration);

            % Perform ray-tracing
            disp(' > Executing ray-tracing ...')
            addpath('../image_sources/', '../image_sources/air_absorption/');
            [image_sources_list, reflectogram] = mirror_images(tempRaytracingFilename);

            % Save clouds to file
            tempVirtualSourcesCloudFile = sprintf('temp_cloud_%d.mat', currentSourcePositionIndex);
            saveCloudToFile(image_sources_list, tempVirtualSourcesCloudFile);

            % Generate the IRs for current cloud -> array
            tempIRsFile = sprintf('temp_irs_%d.mat', currentSourcePositionIndex);
            addpath('../array_simulation/');
            execute_array_simulation(tempVirtualSourcesCloudFile, tempIRsFile, MAX_REFLECTION_ORDER, SAMPLE_RATE, PROPAGATION_SPEED, AMOUNT_OF_SENSORS, ...
                                     ARRAY_CENTER, ARRAY_RADIUS);

            % Perform IR decomposition
            decomposedCloudFilename = sprintf('temp_decomposed_cloud_%d.mat', currentSourcePositionIndex);
            addpath('../ir_decomposition/');
            decomposeIR(tempIRsFile, decomposedCloudFilename, PROPAGATION_SPEED, SAMPLE_RATE, ARRAY_CENTER, []);

            % Reconstruct sound field image and save data to enumerated file
            powerMapFilename = sprintf('temp_map_%d.mat', currentSourcePositionIndex);
            addpath('../room_geometry_estimator/');
            reconstructSoundField({decomposedCloudFilename}, PROPAGATION_SPEED, POWER_MAP_FREQUENCIES, powerMapFilename, x_limits, y_limits, POWER_MAP_GRID_RESOLUTION);

%             figure
%             decomposedCloud = load(decomposedCloudFilename);
%             rayTracingCloud = load(tempVirtualSourcesCloudFile);
%             subplot(2,1,1)
%             plot(rayTracingCloud.image_sources_list(:,1), rayTracingCloud.image_sources_list(:,2), 'o')
%             hold on
%             grid on
%             plot(decomposedCloud.image_sources_list(:,1), decomposedCloud.image_sources_list(:,2), 'x')
%             legend({'Ray tracing sources', 'Decomposed IR sources'})
%             subplot(2,1,2)
%             reconstructedSoundField = load(powerMapFilename);
%             imagesc(reconstructedSoundField.x_coordinates, reconstructedSoundField.y_coordinates, reconstructedSoundField.amplitudeMap)
%             set(gca,'YDir','normal')
%             xlabel('X-coordinate [m]')
%             ylabel('Y-coordinate [m]')
%             pause


        end
    end
% return

    % Average all data and create animation
    disp(sprintf(' > Generating animation from all computed power maps ...'));
    [averagedPowerMap, x_coordinates, y_coordinates] = createAnimation(SOURCE_POSITIONS, ARRAY_CENTER, ROOM_BOUNDARIES);
    disp(sprintf('\b\b [DONE]'));
    save('averaged_power_map.mat', 'averagedPowerMap', 'x_coordinates', 'y_coordinates', 'SOURCE_POSITIONS', 'ARRAY_CENTER', 'ROOM_BOUNDARIES');


        imagesc(x_coordinates, y_coordinates, averagedPowerMap);
        hold on
        plot(ROOM_BOUNDARIES(:,1), ROOM_BOUNDARIES(:,2), 'xr', 'LineWidth', 2);
        plot(ARRAY_CENTER(1), ARRAY_CENTER(2), 'o', 'MarkerSize', 10, 'Color', 'red', 'LineWidth', 3);
        for index = 1 : size(SOURCE_POSITIONS, 1)
            plot(SOURCE_POSITIONS(index,1), SOURCE_POSITIONS(index,2), 'xb', 'LineWidth', 3)
        end
        set(gca,'YDir','normal')
        xlabel('X-coordinate [m]')
        ylabel('Y-coordinate [m]')

    return


    % Capture the corners using the image
    enteredCorners = captureCorners(x_coordinates, y_coordinates, averagedPowerMap);

    estimatedBoundaries = zeros(size(enteredCorners, 1), 4);
    for index = 1 : size(enteredCorners, 1)
        if index ~= size(enteredCorners, 1)
            estimatedBoundaries(index, :) = [enteredCorners(index, :) enteredCorners(index + 1, :)];
        else
            estimatedBoundaries(index, :) = [enteredCorners(index, :) enteredCorners(1, :)];
        end
    end
    save('clicked_points.mat', 'estimatedBoundaries', 'enteredCorners');
    return

    % Find out where is the real source for each source position
    estimatedSourcePositions = zeros(size(SOURCE_POSITIONS, 1), 2);
    for currentCloudIndex = 1 : size(SOURCE_POSITIONS, 1)

        % Load it's associated cloud
        loadedDecomposedCloudData = load(sprintf('temp_decomposed_cloud_%d.mat', currentCloudIndex));
        currentCloud = loadedDecomposedCloudData.image_sources_list;

        % Find out which one is the original source
        bestCandidateDistance = inf;
        for tempCloudIndex = 1 : size(currentCloud, 1)
            if norm(currentCloud(tempCloudIndex, 1 : 2) - ARRAY_CENTER) < bestCandidateDistance
                bestCandidateDistance = norm(currentCloud(tempCloudIndex, 1 : 2) - ARRAY_CENTER);
                bestCandidateIndex = tempCloudIndex;
            end
        end

        % Save real source position
        estimatedSourcePositions(currentCloudIndex, :) = currentCloud(bestCandidateIndex, 1 : 2);
    end
%     plot(realSourcePositions(:, 1), realSourcePositions(:, 2), 'o')
%     grid
%     hold on
%     plot(ARRAY_CENTER(1), ARRAY_CENTER(2), 'x');
%     for index = 1 : size(ROOM_BOUNDARIES, 1)
%         plot([ROOM_BOUNDARIES(index, 1) ROOM_BOUNDARIES(index, 3)], [ROOM_BOUNDARIES(index, 2) ROOM_BOUNDARIES(index, 4)]);
%     end

    for currentBoundaryIndex = 1 : size(estimatedBoundaries, 1)
        
        for currentSourcePositionIndex = 1 : size(estimatedSourcePositions, 1)
        
            % Check if it fulfils the "geometric criterium" for first order
            % reflection

            % - Find image position            
            imagePosition = findImagePosition(estimatedSourcePositions(currentSourcePositionIndex, :), estimatedBoundaries(currentBoundaryIndex, :));

            % - Test cone condition
            if ~testPointIsWithinCone(imagePosition, estimatedBoundaries(currentBoundaryIndex, :), ARRAY_CENTER)
                % There's no first order reflection in this case
                break;
            end

            % If it does, get the time window that is closest to the
            % associated time window of the event

            % - Find time windows of current IRs
            
            % 1 - For the image source
            associatedDistance_ImageSource = norm(imagePosition - ARRAY_CENTER);
            associatedTime_ImageSource = associatedDistance_ImageSource / PROPAGATION_SPEED * SAMPLE_RATE;

            % 2 - For the real source
            associatedDistance_RealSource = norm(estimatedSourcePositions(currentSourcePositionIndex) - ARRAY_CENTER);
            associatedTime_RealSource = associatedDistance_RealSource / PROPAGATION_SPEED * SAMPLE_RATE;

            loadedIRs = load(sprintf('temp_irs_%d.mat', currentBoundaryIndex));            
            timeWindowsCenterTime = findTimeWindows(loadedIRs.microphonesData);

            % 1 - Find the closest time window to the image source time
            closestTimeWindowIndex_ImageSource = 1;
            closestTimeWindowDiff_ImageSource = inf;
            for index = 1 : length(timeWindowsCenterTime)
                if abs((timeWindowsCenterTime(index) - associatedTime_ImageSource)) < closestTimeWindowDiff_ImageSource
                    closestTimeWindowDiff_ImageSource = abs((timeWindowsCenterTime(index) - associatedTime_ImageSource));
                    closestTimeWindowIndex_ImageSource = index;
                end
            end

            % 2 - Find the closest time window to the real source time
            closestTimeWindowIndex_RealSource = 1;
            closestTimeWindowDiff_RealSource = inf;
            for index = 1 : length(timeWindowsCenterTime)
                if abs((timeWindowsCenterTime(index) - associatedTime_RealSource)) < closestTimeWindowDiff_RealSource
                    closestTimeWindowDiff_RealSource = abs((timeWindowsCenterTime(index) - associatedTime_RealSource));
                    closestTimeWindowIndex_RealSource = index;
                end
            end

            % Load it's associated cloud (to know the used time-window size)
            loadedDecomposedCloudData = load(sprintf('temp_decomposed_cloud_%d.mat', currentSourcePositionIndex));

            % Perform the detection and estimate the spectrum of the
            % virtual source by its angle (respect to the array position)

            % 1 - Get the time window associated to the detected image source
            timeWindow_ImageSource = getTimeWindowContent(loadedIRs.microphonesData, ...
                                                   timeWindowsCenterTime(closestTimeWindowIndex_ImageSource), ...
                                                   loadedDecomposedCloudData.timeWindowSize);

            % 1 - Detect the amount of signals
            detectedAmountOfSignals_ImageSource = detectAmountOfSignals(timeWindow_ImageSource);

            detectedAmountOfSignals_ImageSource = 1;
            
            % 1 - Perform an estimation of the image source signal
            associatedDOA_Image = atan( (imagePosition(2) - ARRAY_CENTER(2)) / (imagePosition(1) - ARRAY_CENTER(1)) );
            estimatedSignal_Image = estimateSignalFromTimeWindow(timeWindow_ImageSource, ...
                                                                 loadedIRs.microphonePositions, ...
                                                                 associatedDistance_ImageSource, ...
                                                                 associatedDOA_Image, ...
                                                                 detectedAmountOfSignals_ImageSource, ...
                                                                 PROPAGATION_SPEED, ...
                                                                 SAMPLE_RATE);




            % 2 - Get the time window associated to the detected real source
            timeWindow_RealSource = getTimeWindowContent(loadedIRs.microphonesData, ...
                                                   timeWindowsCenterTime(closestTimeWindowIndex_RealSource), ...
                                                   loadedDecomposedCloudData.timeWindowSize);

            % 2 - Detect the amount of signals (should be only one in most cases)
            detectedAmountOfSignals_RealSource = detectAmountOfSignals(timeWindow_RealSource);

            associatedDOA_Real = atan( (estimatedSourcePositions(currentSourcePositionIndex, 1) - ARRAY_CENTER(2)) / ...
                                       (estimatedSourcePositions(currentSourcePositionIndex, 1) - ARRAY_CENTER(1)) );
            estimatedSignal_Real = estimateSignalFromTimeWindow(timeWindow_RealSource, ...
                                                                 loadedIRs.microphonePositions, ...
                                                                 associatedDistance_RealSource, ...
                                                                 associatedDOA_Real, ...
                                                                 detectedAmountOfSignals_RealSource, ...
                                                                 PROPAGATION_SPEED, ...
                                                                 SAMPLE_RATE);


%             estimatedSignal_Real = estimatedSignal_Real / max(abs(estimatedSignal_Real));
%             estimatedSignal_Image = estimatedSignal_Image / max(abs(estimatedSignal_Image));
% 
%             subplot(2,1,1)
%             plot(real(estimatedSignal_Real))
%             grid on
%             subplot(2,1,2)
%             plot(1 - real((estimatedSignal_Image)))
%             grid on
%             return








            % Compensate estimated espectrums for distance attenuation
            % (both direct and reflected) so the differences are only due
            % to the boundary's absorption
%             distanceCompensationRatio = norm(imagePosition - ARRAY_CENTER) / ...
%                                         norm(estimatedSourcePositions(currentSourcePositionIndex, :) - ARRAY_CENTER);
%             estimatedSignal_Image = estimatedSignal_Image * distanceCompensationRatio;
%             normalizationFactor = max(abs(estimatedSignal_Real));
%             estimatedSignal_Real = estimatedSignal_Real / normalizationFactor;
%             estimatedSignal_Image = estimatedSignal_Image / normalizationFactor;

%             randomNoise = randn(1, SAMPLE_RATE);
%             estimatedSignal_Real = conv(estimatedSignal_Real, randomNoise);
%             estimatedSignal_Image = conv(estimatedSignal_Image, randomNoise);
            
            
            
            % Compare spectrums to estimate the absorption characteristics
            estimatedSpectrum_Image = 20*log10(abs(fft(estimatedSignal_Image)));
            estimatedSpectrum_Real = 20*log10(abs(fft(estimatedSignal_Real)));
            spectrumDifference = estimatedSpectrum_Real - estimatedSpectrum_Image;
            estimatedSpectrum_Real = estimatedSpectrum_Real(1 : round(length(spectrumDifference)/2));
            estimatedSpectrum_Image = estimatedSpectrum_Image(1 : round(length(spectrumDifference)/2));
            spectrumDifference = spectrumDifference(1 : round(length(spectrumDifference)/2));

            semilogx(estimatedSpectrum_Real)
            hold on
            semilogx(estimatedSpectrum_Image)
            semilogx(spectrumDifference)
            grid on
            legend({'Real', 'Image', 'Difference '})
            return
            
            % Save estimation for current boundary ...
        end
        
        
    end







    return
    
    % For each boundary, find the position where the virtual source should
    % be found

    for currentBoundaryIndex = 1 : size(estimatedBoundaries, 1)
        
        % For each source
        for currentSourcePositionIndex = 1 : size(SOURCE_POSITIONS, 1)
            SOURCE_POSITIONS(index,:)
        end
    end
    

    % Finding the places of the 1st order image sources for each boundary

    
    
    
    
    
    
    
    
    
    
    
    
    
    


    return














% 
%     filenames = {};
%     for index = 1 : size(SOURCE_POSITIONS, 1)
%         filenames(length(filenames) +1) = {sprintf('temp_map_%d.mat', index)};
%     end
%     plotAmplitudeMap_image_fast(filenames, [], ...
%                             POWER_MAP_FREQUENCIES, PROPAGATION_SPEED, -1);
%     return
%     
%     
%     
%     
%     
%     
%     
%     propagationSpeed = PROPAGATION_SPEED;
%     gridResolution = loadedMapData.gridResolution;
%     inputSignals_coordinates = [];
%     inputSignals_freq = loadedMapData.inputSignals_freq;
%     amplitudeMap = averagedPowerMap;
%     x_coordinates = loadedMapData.x_coordinates;
%     y_coordinates = loadedMapData.y_coordinates;
%     save('averagedMap.mat', 'propagationSpeed', 'gridResolution', 'inputSignals_coordinates', ...
%          'inputSignals_freq', 'amplitudeMap', 'x_coordinates', 'y_coordinates');
%     plotAmplitudeMap_image_fast({'averagedMap.mat'}, [], ...
%                                 POWER_MAP_FREQUENCIES, PROPAGATION_SPEED, -1);
% 

end

function enteredCorners = captureCorners(x_coordinates, y_coordinates, averagedPowerMap)

    enteredCorners = [];
    while 1

        clf();
        imagesc(x_coordinates, y_coordinates, averagedPowerMap);
        set(gca,'YDir','normal')
        xlabel('X-coordinate [m]');
        ylabel('Y-coordinate [m]');

        if ~isempty(enteredCorners)
            hold on
            plot(enteredCorners(:, 1), enteredCorners(:, 2), 'x', 'Color', 'red', 'LineWidth', 2);
        end

        title(sprintf('Left click -> Enter new corner.\nRight click -> Finish.'))

        [x, y, button] = ginput(1);
        if button == 3
            break;
        end
        enteredCorners(size(enteredCorners, 1) + 1, :) = [x, y];
    end
    close all

end

function [averagedPowerMap, x_coordinates, y_coordinates] = createAnimation(sourcePositions, arrayCenter, roomBoundaries)
    averagedPowerMap = [];
    frameCounter = 1;
    figureHandler = figure();
    for index = 1 : size(sourcePositions, 1)
        powerMapFilename = sprintf('temp_map_%d.mat', index);
        loadedMapData = load(powerMapFilename);
        
        loadedMapData.amplitudeMap = limitDynamicRange(loadedMapData.amplitudeMap, -55, -10);
        
        if isempty(averagedPowerMap)
            averagedPowerMap = loadedMapData.amplitudeMap;
        else
            averagedPowerMap = averagedPowerMap * ((index-1)/index) + (loadedMapData.amplitudeMap/index);
        end
        averagedPowerMap = averagedPowerMap - max(max(averagedPowerMap));
        minValue = min(min(averagedPowerMap));
        averagedPowerMap = limitDynamicRange(averagedPowerMap, minValue * 18.9 / 20, 0);

        % TODO: May need to limit the dynamic range here
        clf
        imagesc(loadedMapData.x_coordinates, loadedMapData.y_coordinates, averagedPowerMap);
        hold on
        plot(roomBoundaries(:,1), roomBoundaries(:,2), 'xr', 'LineWidth', 2);
        plot(arrayCenter(1), arrayCenter(2), 'o', 'MarkerSize', 10, 'Color', 'red', 'LineWidth', 3);
        set(gca,'YDir','normal')
        xlabel('X-coordinate [m]')
        ylabel('Y-coordinate [m]')

        % Create a new image with the current average and save it for the animation
        currentFrame = getframe(figureHandler);
        if frameCounter == 1 
            [imageToSave, colorsMap] = rgb2ind(currentFrame.cdata, 256, 'nodither');
        else
            [imageToSave(:, :, 1, frameCounter), colorsMap] = rgb2ind(currentFrame.cdata, colorsMap, 'nodither');
        end
        frameCounter = frameCounter + 1;
    end

    % Save animation to file
    fileNumber = 1;
    filenameToSave = sprintf('animation_%d.gif', fileNumber);
    while exist(filenameToSave, 'file')
        fileNumber = fileNumber + 1;
        filenameToSave = sprintf('animation_%d.gif', fileNumber);
    end
    imwrite(imageToSave, colorsMap, filenameToSave, 'DelayTime', 0, 'LoopCount', inf);
    close all

    x_coordinates = loadedMapData.x_coordinates;
    y_coordinates = loadedMapData.y_coordinates;
end

function outputData = limitDynamicRange(inputData, minimumThreshold, maximumThreshold)

    if min(min(inputData)) >= maximumThreshold
        maximumThreshold = 0;
    end
    for index1 = 1 : size(inputData, 1)
        for index2 = 1 : size(inputData, 2)
            if inputData(index1, index2) < minimumThreshold
                inputData(index1, index2) = minimumThreshold;
            else
                if inputData(index1, index2) > maximumThreshold
                    inputData(index1, index2) = maximumThreshold;
                end
            end
        end
    end
    outputData = inputData - max(max(inputData));

end

function saveCloudToFile(image_sources_list, tempVirtualSourcesCloudFile)

    save(tempVirtualSourcesCloudFile, 'image_sources_list');

end

function generateTempFileForRayTracing(room_boundaries, mic_position, source_position, boundary_collision_error_margin, include_real_source, ...
                                         max_reflection_order, mic_collision_error_margin, propagation_speed, sample_rate, scan_angle_max, ...
                                         scan_angle_min, scan_angle_resolution, boundaries_reflection_coefficients, boundaries_reflection_frequencies, ...
                                         tempFilename, show_animated_plot, animated_plot_pause_duration)

    % Save variables values to a file
    save(tempFilename, ...
                         'source_position', ...
                         'mic_position', ...
                         'max_reflection_order', ...
                         'mic_collision_error_margin', ...
                         'room_boundaries', ...
                         'boundaries_reflection_coefficients', ...
                         'boundaries_reflection_frequencies', ...
                         'scan_angle_min', ...
                         'scan_angle_resolution', ...
                         'scan_angle_max', ...
                         'show_animated_plot', ...
                         'animated_plot_pause_duration', ...
                         'include_real_source', ...
                         'sample_rate', ...
                         'propagation_speed', ...
                         'boundary_collision_error_margin');

end

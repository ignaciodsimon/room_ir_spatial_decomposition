function showWithAndWithoutWalls()

    PROPAGATION_SPEED = 343;
    POWER_MAP_FREQUENCIES = 250;%primes(2000);
    POWER_MAP_GRID_RESOLUTION = 0.1;
    x_limits = [-2 7];
    y_limits = [-2 11];

    noiseAmplitude = 1.0/2;

    rayTracingData = load('ray_tracing_1.mat');

    withWallsFilename = 'temp_cloud_1.mat';
    cleanCloud = load(withWallsFilename);

    % Create cloud with added noise
    withoutWallsFilename = 'noisy_cloud.mat';
    rng(1000);
%     cleanCloud.image_sources_list(:, 1 : 2) = cleanCloud.image_sources_list(:, 1 : 2) + randn(size(cleanCloud.image_sources_list(:, 1 : 2))) * noiseAmplitude;
    image_sources_list = cleanCloud.image_sources_list(8,:);
    save(withoutWallsFilename, 'image_sources_list');

    % Estimate power map for clean cloud
    cleanPowerMapFilename = 'generated_map_walls.mat';
    addpath('../../../room_geometry_estimator/');
    reconstructSoundField({withWallsFilename}, PROPAGATION_SPEED, POWER_MAP_FREQUENCIES, ...
                          cleanPowerMapFilename, x_limits, y_limits, POWER_MAP_GRID_RESOLUTION);

    % Estimate power map for noisy cloud
    noisyPowerMapFilename = 'generated_map_no_walls.mat';
    addpath('../../../room_geometry_estimator/');
    reconstructSoundField({withoutWallsFilename}, PROPAGATION_SPEED, POWER_MAP_FREQUENCIES, ...
                          noisyPowerMapFilename, x_limits, y_limits, POWER_MAP_GRID_RESOLUTION);

    % Load created maps
    cleanMapData = load(cleanPowerMapFilename);
    noisyMapData = load(noisyPowerMapFilename);

    cleanMapData.amplitudeMap = cleanMapData.amplitudeMap - max(max(cleanMapData.amplitudeMap));
    
    subplot(1,2,1)
    imagesc(cleanMapData.x_coordinates, cleanMapData.y_coordinates, cleanMapData.amplitudeMap)
    hold on
    set(gca,'YDir','normal')
    for index = 1 : size(rayTracingData.room_boundaries, 1)
        plot([rayTracingData.room_boundaries(index, 1) rayTracingData.room_boundaries(index, 3)], ...
             [rayTracingData.room_boundaries(index, 2) rayTracingData.room_boundaries(index, 4)], ...
             '--', 'LineWidth', 3.5, 'Color', [0.8 .1 .4 1]);
    end
    caxis([-45 -10])
    xlabel('X-coordinate [m]')
    ylabel('Y-coordinate [m]')
    title(sprintf('With boundaries, freq: %d [Hz]', POWER_MAP_FREQUENCIES))

    subplot(1,2,2)
    imagesc(noisyMapData.x_coordinates, noisyMapData.y_coordinates, noisyMapData.amplitudeMap)
    hold on
    set(gca,'YDir','normal')
    xlabel('X-coordinate [m]')
    ylabel('Y-coordinate [m]')
%     caxis([-45.5 -33])
    title(sprintf('Without boundaries, freq: %d [Hz]', POWER_MAP_FREQUENCIES))

    % Save plot to file
    set(gcf, 'PaperPosition', [-1.0 -0.35 13.4 8.05]);
    set(gcf, 'PaperSize', [11.5 7.5]);
    saveas(gcf, 'comparison_with_without_walls.pdf', 'pdf');
    close all
    
end

function executeDemonstration()

    PROPAGATION_SPEED = 343;
    POWER_MAP_FREQUENCIES = primes(2000);
    POWER_MAP_GRID_RESOLUTION = 0.1;
    x_limits = [-2 7];
    y_limits = [-2 11];

    noiseAmplitude = 1.0/2;
    
    rayTracingData = load('ray_tracing_1.mat');

    noNoiseCloudFilename = 'temp_cloud_1.mat';
    cleanCloud = load(noNoiseCloudFilename);

    % Create cloud with added noise
    noisyCloudFilename = 'noisy_cloud.mat';
    rng(1000);
    cleanCloud.image_sources_list(:, 1 : 2) = cleanCloud.image_sources_list(:, 1 : 2) + randn(size(cleanCloud.image_sources_list(:, 1 : 2))) * noiseAmplitude;
    image_sources_list = cleanCloud.image_sources_list;
    save(noisyCloudFilename, 'image_sources_list');

    % Estimate power map for clean cloud
    cleanPowerMapFilename = 'generated_map_no_noise.mat';
    addpath('../../../room_geometry_estimator/');
    reconstructSoundField({noNoiseCloudFilename}, PROPAGATION_SPEED, POWER_MAP_FREQUENCIES, ...
                          cleanPowerMapFilename, x_limits, y_limits, POWER_MAP_GRID_RESOLUTION);

    % Estimate power map for noisy cloud
    noisyPowerMapFilename = 'generated_map_noisy.mat';
    addpath('../../../room_geometry_estimator/');
    reconstructSoundField({noisyCloudFilename}, PROPAGATION_SPEED, POWER_MAP_FREQUENCIES, ...
                          noisyPowerMapFilename, x_limits, y_limits, POWER_MAP_GRID_RESOLUTION);

    % Load created maps
    cleanMapData = load(cleanPowerMapFilename);
    noisyMapData = load(noisyPowerMapFilename);

    subplot(1,2,1)
    imagesc(cleanMapData.x_coordinates, cleanMapData.y_coordinates, cleanMapData.amplitudeMap)
    hold on
    set(gca,'YDir','normal')
    xlabel('X-coordinate [m]')
    ylabel('Y-coordinate [m]')
    for index = 1 : size(rayTracingData.room_boundaries, 1)
        plot([rayTracingData.room_boundaries(index, 1) rayTracingData.room_boundaries(index, 3)], ...
             [rayTracingData.room_boundaries(index, 2) rayTracingData.room_boundaries(index, 4)], ...
             '--', 'LineWidth', 2.5, 'Color', [0.8 .1 .4 1]);
    end
    caxis([-44 -30])
    title('No virtual source location error')

    subplot(1,2,2)
    imagesc(noisyMapData.x_coordinates, noisyMapData.y_coordinates, noisyMapData.amplitudeMap)
    hold on
    set(gca,'YDir','normal')
    xlabel('X-coordinate [m]')
    ylabel('Y-coordinate [m]')
    for index = 1 : size(rayTracingData.room_boundaries, 1)
        plot([rayTracingData.room_boundaries(index, 1) rayTracingData.room_boundaries(index, 3)], ...
             [rayTracingData.room_boundaries(index, 2) rayTracingData.room_boundaries(index, 4)], ...
             '--', 'LineWidth', 2.5, 'Color', [0.8 .1 .4 1]);
    end
    caxis([-45.5 -33])
    title(sprintf('Location tolerance of ±%.3f [m]', noiseAmplitude))

    % Save plot to file
    set(gcf, 'PaperPosition', [-1.0 -0.35 13.4 8.05]);
    set(gcf, 'PaperSize', [11.5 7.5]);
    saveas(gcf, sprintf('geometry_comparison_error_%.3f.pdf', noiseAmplitude), 'pdf');
    close all
    
end

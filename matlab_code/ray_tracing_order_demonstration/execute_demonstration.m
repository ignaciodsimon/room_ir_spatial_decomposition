function execute_demonstration()

    roomGeometryFilename = 'room_geometry_definition.mat';

    % Execute the room geometry capture utility
    disp(' > Launching geometry capture utility ...')
%     import_room_geometry(roomGeometryFilename);

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
    MAX_REFLECTION_ORDER = 50;
    MIC_COLLISION_ERROR_MARGIN = 0.10;
    PROPAGATION_SPEED = 343;
    SAMPLE_RATE = 48000;
    SCAN_ANGLE_MAX = 180;
    SCAN_ANGLE_MIN = -180;
    SCAN_ANGLE_RES = 0.25;
    AMOUNT_OF_SENSORS = 1;
    ARRAY_RADIUS = 0.075;
    POWER_MAP_FREQUENCIES = primes(2000);
    x_limits = [-2 12];
    y_limits = [-2 12];
    NOISE_FLOOR = -55;

%     % Generate synthetic materials with the same absorption for all frequencies
%     disp(' > Generating synthetic materials ...')
%     ABSORPTION_FREQUENCIES = [125, 250, 500, 1000, 2000, 4000, 8000, 16000, 24000];
%     ABSORPTION_COEFFICIENTS = ones(length(loadedGeometry.boundaries_reflection_coefficients), length(ABSORPTION_FREQUENCIES)) * 0.90;
%     for index = 1 : size(ABSORPTION_COEFFICIENTS, 1)
%         ABSORPTION_COEFFICIENTS(index, :) = ABSORPTION_COEFFICIENTS(index, :) * loadedGeometry.boundaries_reflection_coefficients(index);
%     end
%     disp(sprintf('\b\b [DONE]'))
% 
%     disp(sprintf('\n[>] Main computation part started.\n'))
% 
%     % Update the new source position
%     currentSourcePosition = SOURCE_POSITIONS(1, :);
% 
%     % Save current simulation vars to file (for ray-tracing)
%     disp(' > Generating temp file for current simulation ...')
%     tempRaytracingFilename = sprintf('ray_tracing.mat');
%     generateTempFileForRayTracing(ROOM_BOUNDARIES, ARRAY_CENTER, currentSourcePosition, BOUNDARY_COLLISION_ERROR_MARGIN, ...
%                                   INCLUDE_REAL_SOURCE, MAX_REFLECTION_ORDER, MIC_COLLISION_ERROR_MARGIN, PROPAGATION_SPEED, ...
%                                   SAMPLE_RATE, SCAN_ANGLE_MAX, SCAN_ANGLE_MIN, SCAN_ANGLE_RES, ABSORPTION_COEFFICIENTS, ABSORPTION_FREQUENCIES, ...
%                                   tempRaytracingFilename);
% 
%     tic
%     % Perform ray-tracing
%     disp(' > Executing ray-tracing ...')
%     addpath('../image_sources/', '../image_sources/air_absorption/');
%     [image_sources_list, reflectogram] = mirror_images(tempRaytracingFilename);
%     save('reflectogram.mat', 'reflectogram');
% 
%     % Save clouds to file
%     tempVirtualSourcesCloudFile = sprintf('temp_cloud.mat');
%     saveCloudToFile(image_sources_list, tempVirtualSourcesCloudFile);
% 
%     % Generate the IRs for current cloud -> array
%     tempIRsFile = sprintf('temp_irs_%d.mat', MAX_REFLECTION_ORDER);
%     addpath('../array_simulation/');
%     execute_array_simulation(tempVirtualSourcesCloudFile, tempIRsFile, MAX_REFLECTION_ORDER, SAMPLE_RATE, PROPAGATION_SPEED, AMOUNT_OF_SENSORS, ...
%                              ARRAY_CENTER, ARRAY_RADIUS); 
%     toc

    % Load generated IRs and reflectogram
%     loadedIRsData = load(sprintf('temp_irs_%d.mat', MAX_REFLECTION_ORDER));
    loadedIRsData = load('temp_irs_50_reflection_0.9.mat');
    loadedIRs = loadedIRsData.microphonesData;
    loadedReflectogram = load('reflectogram.mat');
    reflectogram = loadedReflectogram.reflectogram;

    % Normalize and add noise
    loadedIRs = (loadedIRs / max(abs(loadedIRs))) + randn(size(loadedIRs))*10^(NOISE_FLOOR/20);
    
    % Generate plot
    plot([1 : 1 : size(loadedIRs, 2)] / SAMPLE_RATE * 1000, loadedIRs')
    grid on
    xlim([0 1000])
    ylim([-0.3 1.1]);
    xlabel('Time [ms]')
    ylabel('Amplitude [.]')
    set(gcf, 'PaperPosition', [-0.5 +0.1 9.7 3.0]);
    set(gcf, 'PaperSize', [8.5 3.0]);
    saveas(gcf, sprintf('ir_generated_max_order_%d.pdf', MAX_REFLECTION_ORDER), 'pdf');
    close all

end

function generateTempFileForRayTracing(room_boundaries, mic_position, source_position, boundary_collision_error_margin, include_real_source, ...
                                         max_reflection_order, mic_collision_error_margin, propagation_speed, sample_rate, scan_angle_max, ...
                                         scan_angle_min, scan_angle_resolution, boundaries_reflection_coefficients, boundaries_reflection_frequencies, ...
                                         tempFilename)


    show_animated_plot = 0;
    animated_plot_pause_duration = 0.2;

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

function saveCloudToFile(image_sources_list, tempVirtualSourcesCloudFile)

    save(tempVirtualSourcesCloudFile, 'image_sources_list');

end

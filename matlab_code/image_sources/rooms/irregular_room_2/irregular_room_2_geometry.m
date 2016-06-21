function irregular_room_geometry(display_geometry)
    % Generates the definitions file needed to execute the simulation. The
    % room created is a trapezoid with four different materials (wall,
    % wood, door, window).
    %
    % Joe.

    % Choose to include the direct ray as well
    include_real_source = 1;

    % Perform an interactive simulation or not
    show_animated_plot = 0;
    animated_plot_pause_duration = 0.01;

    % Simulation limit
    max_reflection_order = 10;

    % Collision precission
    boundary_collision_error_margin = 0.05;
    mic_collision_error_margin = 0.1;

    % General constants
    propagation_speed = 343;
    sample_rate = 48000;

    mic_position = [5, 6];
    source_position = [2, 3];

    % Scanning limits of the simulation
    scan_angle_max = 180;
    scan_angle_min = -180;
    scan_angle_resolution = 0.25;

    % Geometry definition
    room_boundaries = [0 0 0 7
                       0 7 2 8
                       2 8 6 8
                       6 8 6 3
                       6 3 0 0];

    % Absorption of boundaries
    boundaries_reflection_coefficients = [0.9900    0.9800    0.9800    0.9800    0.9700    0.9600  % Wall
                                          0.9900    0.9800    0.9800    0.9800    0.9700    0.9600  % Wall
                                          0.8400    0.8600    0.8900    0.9200    0.9200    0.9300  % Wood
                                          0.9968    0.9990    0.9996    0.9998    0.9997    0.9997  % Door
                                          0.6500    0.7500    0.8200    0.8800    0.9300    0.9600  % Window
                                          ];
    boundaries_reflection_frequencies = [125 250 500 1000 2000 4000];

    % Save data to .mat file
    save('irregular_room_geometry_2.mat', 'include_real_source', 'max_reflection_order', ...
         'boundary_collision_error_margin', 'mic_collision_error_margin', ...
         'propagation_speed', 'sample_rate', 'mic_position', 'source_position', ...
         'scan_angle_max', 'scan_angle_min', 'scan_angle_resolution', ...
         'room_boundaries', 'boundaries_reflection_coefficients', ...
         'boundaries_reflection_frequencies', 'show_animated_plot', ...
         'animated_plot_pause_duration');

    % Show generated room
    if display_geometry
        display_geometry_from_file('irregular_room_geometry_2.mat');
    end
end

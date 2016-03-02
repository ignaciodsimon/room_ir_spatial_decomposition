function synthetic_room_geometry(display_geometry)
    % Generates the definitions file needed to execute the simulation. The
    % room created is rectangular with two windows, a shelf and a door.
    % Different materials are considered: brick wall, wood, window and
    % door.
    %
    % Joe.

    % Choose to include the direct ray as well
    include_real_source = 1;

    % Perform an interactive simulation or not
    show_animated_plot = 0;
    animated_plot_pause_duration = 0.01;

    % Simulation limit
    max_reflection_order = 50;

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
    room_boundaries = [0 0 0 8          % Left wall (concrete)
                       0 8 7 8          % Top wall (wood)
                       7 8 7 7          % Right wall (concrete)
                       7 7 7 6          % Door
                       7 6 7 4          % Right wall (concrete)
                       7 4 8 4          % Closet (wood)
                       8 4 8 1
                       8 1 7 1
                       7 1 7 0          % Right wall (concrete)
                       7 0 6 0          % Bottom wall (concrete)
                       6 0 6 0.25       % Border of window (concrete)
                       6 0.25 4 0.25    % Window right
                       4 0.25 4 0       % Border of window (concrete)
                       4 0 3 0          % Bottom wall (concrete)
                       3 0 3 0.25       % Border of window (concrete)
                       3 0.25 1 0.25    % Window left
                       1 0.25 1 0       % Border of window (concrete)
                       1 0 0 0          % Bottom wall (concrete)
                       ];

    % Absorption of boundaries
    boundaries_reflection_coefficients = [0.9900    0.9800    0.9800    0.9800    0.9700    0.9600  % Wall
                                          0.8400    0.8600    0.8900    0.9200    0.9200    0.9300  % Wood
                                          0.9900    0.9800    0.9800    0.9800    0.9700    0.9600  % Wall
                                          0.9968    0.9990    0.9996    0.9998    0.9997    0.9997  % Door
                                          0.9900    0.9800    0.9800    0.9800    0.9700    0.9600  % Wall
                                          0.8400    0.8600    0.8900    0.9200    0.9200    0.9300  % Wood
                                          0.8400    0.8600    0.8900    0.9200    0.9200    0.9300  % Wood
                                          0.8400    0.8600    0.8900    0.9200    0.9200    0.9300  % Wood
                                          0.9900    0.9800    0.9800    0.9800    0.9700    0.9600  % Wall
                                          0.9900    0.9800    0.9800    0.9800    0.9700    0.9600  % Wall
                                          0.9900    0.9800    0.9800    0.9800    0.9700    0.9600  % Wall
                                          0.6500    0.7500    0.8200    0.8800    0.9300    0.9600  % Window
                                          0.9900    0.9800    0.9800    0.9800    0.9700    0.9600  % Wall
                                          0.9900    0.9800    0.9800    0.9800    0.9700    0.9600  % Wall
                                          0.9900    0.9800    0.9800    0.9800    0.9700    0.9600  % Wall
                                          0.6500    0.7500    0.8200    0.8800    0.9300    0.9600  % Window
                                          0.9900    0.9800    0.9800    0.9800    0.9700    0.9600  % Wall
                                          0.9900    0.9800    0.9800    0.9800    0.9700    0.9600  % Wall
                                          ];
    boundaries_reflection_frequencies = [125 250 500 1000 2000 4000];

    % Save data to .mat file
    save('synthetic_room_geometry.mat', 'include_real_source', 'max_reflection_order', ...
         'boundary_collision_error_margin', 'mic_collision_error_margin', ...
         'propagation_speed', 'sample_rate', 'mic_position', 'source_position', ...
         'scan_angle_max', 'scan_angle_min', 'scan_angle_resolution', ...
         'room_boundaries', 'boundaries_reflection_coefficients', ...
         'boundaries_reflection_frequencies', 'show_animated_plot', ...
         'animated_plot_pause_duration');

    % Show generated room
    if display_geometry
        display_geometry_from_file('synthetic_room_geometry.mat');
    end
end

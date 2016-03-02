function mirror_image_test()
    % Temporary script to test the old prototype of the mirror_image
    % function with a simple set of parameters. No longer used. Keep for
    % possible basic checks.
    %
    % Joe.

    source_position = [0, 1];                   % Position of source and receiver on the room
    mic_position = [0 0];                      %
                                                %
    max_reflection_order = 200;                   % Maximum amount of times a ray can bounce before we consider it extinguished
                                                %
    mic_collision_error_margin = 0.05;         % Collision margin for rays approaching the microphone
                                                %
    boundary_collision_error_margin = 0.02;    % Collision margin for rays approaching a boundary
                                                %
    irregular_room = [-2  -4  -2.5 5            % Definition of boundaries of the room, they dont need to be in
                      -2  -4   2  -2            % order, just be a closed contour.
                      -2.5 5   1.2   3.2        %
                       2  -2   1.2   3.2];      %
    square_room = [-2  -4  -2   2               %
                   -2  -4   2  -4               %
                   -2   2   2   2               %
                    2  -4   2   2];             %
    room_boundaries = square_room;              %
%     room_boundaries = circular_room(5, 20);
                                                %
    boundaries_reflection_coefficients = ones(size(room_boundaries, 1), 1);
%     boundaries_reflection_coefficients = [.7   % Reflection coefficients of each boundary
%                                           .8   %
%                                           .6   %
%                                           .5]; %
    scan_angle_min = -180;                      % Parameters for the angle scanning resolution
    scan_angle_resolution = 0.5;                %
    scan_angle_max = 180;                       %
                                                %
    show_animated_plot = 0;                     % Open a plot window that updates with all steps of each iteration
                                                %
    include_real_source = 1;                    % To include all rays that go directly from source to mic
                                                % There might be several, depending on the margin given by
                                                % "mic_collision_error_margin"
                                                %
%     boundary_intermediate_points_count = 4000;   % For the collision check
                                                %
%     ray_intermediate_points_count = 4000;        % For the ray launch
                                                %
    sample_rate = 48000;                        % For the reflectogram
    propagation_speed = 343;                    %

    [image_sources_list, reflectogram] = mirror_images(source_position,  mic_position, max_reflection_order, mic_collision_error_margin, ...
                                                                boundary_collision_error_margin, room_boundaries, boundaries_reflection_coefficients, ...
                                                                scan_angle_min, scan_angle_resolution, scan_angle_max, show_animated_plot, include_real_source, ...
                                                                sample_rate, propagation_speed);

    time_axis = [0 : 1/sample_rate : (length(reflectogram)-1)/sample_rate];

    figure()
    subplot(2,1,1)
    plot(mic_position(1), mic_position(2), 'o')
    hold on
    plot(source_position(1), source_position(2), 'x')
    for boundary = room_boundaries'
        plot([boundary(1) boundary(3)], [boundary(2) boundary(4)], 'Color', [.2 .8 .3])
    end
    for i = 1 : size(image_sources_list, 1)
        plot(image_sources_list(i,1), image_sources_list(i,2), 'o', 'MarkerFaceColor', [.2 .3 .8], 'Color', [1 .5 .3], 'MarkerSize', (1/40 + image_sources_list(i,3)) * 40);
    end
    xlabel('X-coordinate [m]')
    ylabel('Y-coordinate [m]')
    title('Image source location / strength');
    grid
    subplot(2,1,2)
    plot(time_axis, reflectogram);
    title('Obtained reflectogram (using fractional delay)')
    ylabel('Amplitude [.]')
    xlabel('Time [s]')
    grid

    save('computed_room_reflections.mat', 'image_sources_list', 'reflectogram', 'time_axis');
end

function boundaries = circular_room(radius, amount_segments)

    boundaries = zeros(amount_segments, 4);
    arc_angle = 360 / 180 * pi / amount_segments;
    for i = 1 : amount_segments

        boundaries(i,:) = [[cos(arc_angle * i) sin(arc_angle * i)*radius] [cos(arc_angle * (i+1)) sin(arc_angle * (i+1))*radius]];

    end
    
%     for boundary = boundaries'
%         plot([boundary(1) boundary(3)], [boundary(2) boundary(4)], 'Color', [.2 .8 .3])
%         hold on
%     end
end



%     source_position = [1, 1];                   % Position of source and receiver on the room
%     mic_position = [0 -2];                      %
%                                                 %
%     max_reflection_order = 10;                   % Maximum amount of times a ray can bounce before we consider it extinguished
%                                                 %
%     mic_collision_error_margin = 0.02;          % Collision margin for rays approaching the microphone
%                                                 %
%     boundary_collision_error_margin = 0.02;     % Collision margin for rays approaching a boundary
%                                                 %
%     irregular_room = [-2  -4  -2.5 3            % Definition of boundaries of the room, they dont need to be in
%                       -2  -4   2  -2            % order, just be a closed contour.
%                       -2.5 3   2   2            %
%                        2  -2   2   2];          %
%     square_room = [-2  -4  -2   2               %
%                    -2  -4   2  -4               %
%                    -2   2   2   2               %
%                     2  -4   2   2];             %
%     room_boundaries = square_room;              % 
%                                                 %
%     boundaries_reflection_coefficients = [0.7   % Reflection coefficients of each boundary
%                                           0.6   %
%                                           0.5   %
%                                           0.4]; %
%     scan_angle_min = -180;                      % Parameters for the angle scanning resolution
%     scan_angle_resolution = 0.5;                %
%     scan_angle_max = 180;                       %
%                                                 %
%     show_animated_plot = 0;                     % Open a plot window that updates with all steps of each iteration
%                                                 %
%     include_real_source = 1;                    % To include all rays that go directly from source to mic
%                                                 % There might be several, depending on the margin given by
%                                                 % "mic_collision_error_margin"
%                                                 %
%     boundary_intermediate_points_count = 400;   % For the collision check
%                                                 %
%     ray_intermediate_points_count = 400;        % For the ray launch
%                                                 %
%     sample_rate = 48000;                        % For the reflectogram
%     propagation_speed = 343;                    %

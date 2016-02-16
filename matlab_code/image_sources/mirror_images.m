function [image_sources_list, reflectogram] = mirror_images(room_definition_filename)
    % Script to find image sources using the ray tracing method.
    %
    % Notes: 
    %     Working for any 2D plane defined by arbitrary boundaries. The
    %     boundaries are not required to form a closed area. Instead, holes
    %     can exist, for example to allow for a (computationally efficient)
    %     completely absorbent material. The presence of inner obstacles is
    %     allowed, but the model does not contemplate transmission, only
    %     reflection, so it has to be used only as an approximation,
    %     specially in low frequency.
    %
    % The room definition is read from a .mat file that must contain the
    % following fields:
    %
    %     loaded_data
    %     source_position
    %     mic_position
    %     max_reflection_order
    %     mic_collision_error_margin
    %     boundary_collision_error_margin
    %     room_boundaries
    %     boundaries_reflection_coefficients
    %     scan_angle_min
    %     scan_angle_resolution
    %     scan_angle_max
    %     show_animated_plot
    %     include_real_source
    %     sample_rate
    %     propagation_speed
    %
    % Joe.

    % Load data from input file
    disp(sprintf(' > Loading data from file "%s"', room_definition_filename));
    loaded_data = load(room_definition_filename);
    source_position = loaded_data.source_position;
    mic_position = loaded_data.mic_position;
    max_reflection_order = loaded_data.max_reflection_order;
    mic_collision_error_margin = loaded_data.mic_collision_error_margin;
    boundary_collision_error_margin = loaded_data.boundary_collision_error_margin;
    room_boundaries = loaded_data.room_boundaries;
    boundaries_reflection_coefficients = loaded_data.boundaries_reflection_coefficients;
    scan_angle_min = loaded_data.scan_angle_min;
    scan_angle_resolution = loaded_data.scan_angle_resolution;
    scan_angle_max = loaded_data.scan_angle_max;
    show_animated_plot = loaded_data.show_animated_plot;
    include_real_source = loaded_data.include_real_source;
    sample_rate = loaded_data.sample_rate;
    propagation_speed = loaded_data.propagation_speed;

    % This stimation needs to be reviewed
    max_radius = find_max_radius(room_boundaries);

    % Range of scanning angles
    angular_range = [scan_angle_min : scan_angle_resolution: scan_angle_max];

    % Resulting data matrices
    image_sources_list = zeros(length(angular_range), 3);
    image_sources_count = 0;
    reflectogram = zeros(1, ceil(max_radius * max_reflection_order / propagation_speed * sample_rate));

    if show_animated_plot
        figure();
    end

    disp(sprintf(' > Throwing ray on direction: -000.0 [deg] ...'));
tic
    initial_ray_angle_index = 1;
    while (initial_ray_angle_index <= length(angular_range))
        initial_ray_angle = angular_range(initial_ray_angle_index);
        initial_ray_angle_index = initial_ray_angle_index + 1;

        % Show the boundaries, mic and source
        if show_animated_plot
            clf
            plot(mic_position(1), mic_position(2), 'o')
            hold on
            plot(source_position(1), source_position(2), 'x')
            for boundary = room_boundaries'
                plot([boundary(1) boundary(3)], [boundary(2) boundary(4)], 'Color', [.2 .8 .3])
            end
        end

        % New ray initial angle
        ray_origin = source_position;
        ray_direction = initial_ray_angle/180*pi;

        % Track to save all steps (bounces)
        current_track = zeros(max_reflection_order, 2);
        current_track_bounces_count = 0;

        % Continue bouncing current ray until limit conditions are reached
        disp(sprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%+06.1f [deg] ...', ray_direction/pi*180));
        amount_of_collisions = 0;
        wall_collision = 0;
        last_touched_boundary = [];
        while 1

            % Check if we reached the microphone with this ray
            [mic_reached, collision_point] = reached_position(ray_origin, ray_direction, mic_position, max_radius, mic_collision_error_margin);
            if mic_reached
                current_track_bounces_count = current_track_bounces_count + 1;
                current_track(current_track_bounces_count, 1) = sqrt((ray_origin(1) - collision_point(1))^2 + (ray_origin(2) - collision_point(2))^2);
                current_track(current_track_bounces_count, 2) = 1;
            else
                % Otherwise detect if it collides with a wall (it should)
                for boundary_index = [1 : size(room_boundaries, 1)]
                    boundary = room_boundaries(boundary_index, :);
                    boundary_limits = [boundary(1) boundary(2); boundary(3) boundary(4)];

                    % Check that the reached boundary is not the one
                    % it's actually on right now
                    if ~isequal(last_touched_boundary, boundary_limits)
                        [wall_collision, collision_point, reflection_direction, trazed_ray] = detect_collision(ray_origin, ray_direction, boundary_limits, max_radius, boundary_collision_error_margin);
                        if wall_collision
                            last_touched_boundary = boundary_limits;
                            amount_of_collisions = amount_of_collisions + 1;

                            current_track_bounces_count = current_track_bounces_count + 1;
                            current_track(current_track_bounces_count, 1) = sqrt((ray_origin(1) - collision_point(1))^2 + (ray_origin(2) - collision_point(2))^2);
                            current_track(current_track_bounces_count, 2) = boundaries_reflection_coefficients(boundary_index);
                            break
                        end
                    end
                end
            end

            % Plot the new segment: from ray_origin -> collision_point
            if (mic_reached || wall_collision) && show_animated_plot
                plot(trazed_ray(:, 1), trazed_ray(:,2), '.', 'Color', [.6 .2 .3], 'MarkerSize', 3);
                xlim([-5 5])
                ylim([-5 5])
%                 plot([ray_origin(1) collision_point(1)], [ray_origin(2) collision_point(2)], '.', 'Color', [.6 .2 .3], 'MarkerSize', 3)
                pause(0.1)
            end

            if wall_collision
                % Update the ray_origin, ray_direction
                ray_origin = collision_point;
                ray_direction = reflection_direction;
            end

            % Finish current ray if we reached the mic
            if mic_reached
                % Discard the real source if told so
                if include_real_source || (amount_of_collisions > 0)

                    % Add the new data to the reflectogram
                    total_distance = sum(current_track(:,1));
                    distance_attenuation = 1 / total_distance;
                    total_attenuation = prod(nonzeros(current_track(:,2))) * distance_attenuation;
                    total_delay = total_distance / propagation_speed;
%                     reflectogram(round(total_delay * sample_rate)) = reflectogram(round(total_delay * sample_rate)) + total_attenuation;
                    fractional_delta = pad_with_zeros( ...
                                            place_fractional_delta(total_delay * sample_rate), ...
                                            length(reflectogram));
                    reflectogram = reflectogram + (fractional_delta * total_attenuation);

                    % Save new image source position and amplitude
                    image_sources_count = image_sources_count +1;
                    image_souce_position = (total_distance * [cos(ray_direction + pi) sin(ray_direction + pi)]) + mic_position;
                    image_sources_list(image_sources_count, :) = [image_souce_position total_attenuation];

                    if show_animated_plot
                        plot(image_souce_position(1), image_souce_position(2), 'x')
                        plot([image_souce_position(1) mic_position(1)], [image_souce_position(2) mic_position(2)]);
                        pause
                    end
                end

                % Skip two degrees after a collision
                initial_ray_angle_index = initial_ray_angle_index + ceil(1 / scan_angle_resolution);
                break
            end

            % Finish current ray because max order was reached
            if amount_of_collisions >= max_reflection_order
                break
            end

            % Sanity check, boundaries should be closed ideally
            if ~mic_reached && ~wall_collision
                disp('** Unexpected. No collision and didnt reach the mic ... It could mean that there''s a hole in the boundaries ...');
                disp(sprintf('Throwing ray on direction: -000.0 [deg] ...'));
                break
            end
        end
    end
toc
    % Trim reflectogram to the maximum useful length
    for i = [0 : length(reflectogram)-1]
        if reflectogram(length(reflectogram) -i) ~= 0
            reflectogram = reflectogram(1 : length(reflectogram) -i);
            break;
        end
    end

    % Trim image sources list to the amount of sources found
    image_sources_list = image_sources_list(1 : image_sources_count, :);

    % Show plot with image sources and reflectogram
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
        plot(image_sources_list(i,1), image_sources_list(i,2), 'o', 'MarkerFaceColor', [.2 .3 .8], 'Color', [1 .5 .3], 'MarkerSize', image_sources_list(i,3) * 40);
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

    % Save results to a file
    save(sprintf('%s_processed.mat', room_definition_filename), 'image_sources_list', 'reflectogram', 'time_axis');

end

function [collision, collision_point, reflection_direction, trazed_ray] = detect_collision(ray_origin, ray_direction, boundary_limits, max_radius, error_margin)

    % Find the boundary intermediate points
    boundary_direction = atan2(boundary_limits(2,2) - boundary_limits(1,2), boundary_limits(2,1) - boundary_limits(1,1));
    boundary_length = sqrt((boundary_limits(2,2) - boundary_limits(1,2))^2 + (boundary_limits(2,1) - boundary_limits(1,1))^2);
    boundary_points_step_distance = 0.01;
    boundary_points = [boundary_limits(1,1) + cos(boundary_direction) * [0 : boundary_points_step_distance : boundary_length]' ...
                       boundary_limits(1,2) + sin(boundary_direction) * [0 : boundary_points_step_distance : boundary_length]' ];

    % Find the ray intermediate points
    ray_points_step_distance = 0.01;
    ray_points = [ray_origin(1) + cos(ray_direction) * [0 : ray_points_step_distance : max_radius]' ...
                  ray_origin(2) + sin(ray_direction) * [0 : ray_points_step_distance : max_radius]'];

    trazed_ray = ray_points;

    % Check if the ray gets close enough to any of the intermediate points of the boundary
    collision = 0;
    collision_point = [];
    reflection_direction = [];
    for i = [1 : length(ray_points)]
        if min(sqrt( (ray_points(i, 1) - boundary_points(:,1)).^2 + (ray_points(i, 2) - boundary_points(:,2)).^2 )) < error_margin
            collision = 1;
            collision_point = ray_points(i, :);
            reflection_direction = 2 * boundary_direction - ray_direction;
            return
        end        
    end
end

function [collision, collision_point] = reached_position(ray_origin, ray_direction, objective_point, max_radius, error_margin)

    % Creates intermediate points for the ray in the given direction
%     ray_points_count = round(max_radius / error_margin * 1.2);
    ray_points_step_distance = 0.01;
    ray_points = [ray_origin(1) + cos(ray_direction) * [0 : ray_points_step_distance : max_radius]' ...
                  ray_origin(2) + sin(ray_direction) * [0 : ray_points_step_distance : max_radius]'];

    % Check if the ray gets close enough to the objective point
    collision = 0;
    collision_point = [];
    for i = [1 : length(ray_points)]
        if min(sqrt( (ray_points(i, 1) - objective_point(:,1)).^2 + (ray_points(i, 2) - objective_point(:,2)).^2 )) < error_margin
            collision = 1;
            collision_point = ray_points(i, :);
        end
    end
end

function max_distance = find_max_radius(room_boundaries)

    % Get all vertices from the boundaries
    vertices = zeros(2*size(room_boundaries, 1), 2);
    for i = [1 : size(room_boundaries, 1)]
        vertices(2*i -1, :) = [room_boundaries(i,1) room_boundaries(i,2)];
        vertices(2*i, :) = [room_boundaries(i,3) room_boundaries(i,4)];
    end

    % Find the maximum distance between vertices
    max_distance = 0;
    for vertex1 = vertices
        for vertex2 = vertices
            if sqrt( (vertex2(1) - vertex1(1))^2 + (vertex2(2) - vertex1(2))^2 ) > max_distance
                max_distance = sqrt( (vertex2(1) - vertex1(1))^2 + (vertex2(2) - vertex1(2))^2 );
            end
        end
    end

    % Not sure why sometimes it requires this 200% / 400% length
    max_distance = max_distance * 2;
end

function padded_vector = pad_with_zeros(input_vector, desired_length)

    if length(input_vector) == desired_length
        padded_vector = input_vector;
        return
    end

    if length(input_vector) > desired_length
        padded_vector = input_vector(1 : desired_length);
    end
    
    if length(input_vector) < desired_length
        padded_vector = zeros(1, desired_length);
        padded_vector(1 : length(input_vector)) = input_vector;
    end

end

function output_ir = place_fractional_delta(fractionalDelay)

    output_ir = zeros(1, round(fractionalDelay * 2));
    for i = 0 : length(output_ir)-1

        arg = i - (fractionalDelay);
        if arg ~= 0
            output_ir(i+1) = sin(pi* arg) / (pi *arg);
        else
%             disp('The delay is actually an integer value.')
            output_ir(i+1) = 1;
        end
    end
end
function room_boundaries = import_room_geometry()
    % Script to generate room definitions. The following parameters are
    % obtained from the user using a graphical input:
    %
    %   - boundaries placement
    %   - boundaries reflection coefficient
    %   - source and microphone positions
    %
    % Other parameters are internally set to values that are known to work
    % well. They are:
    %
    %     max_reflection_order = 100;
    %     mic_collision_error_margin = 0.05;
    %     scan_angle_min = -180;
    %     scan_angle_resolution = 0.5;
    %     scan_angle_max = 180;
    %     show_animated_plot = 0;
    %     include_real_source = 1;
    %     sample_rate = 48000;
    %     propagation_speed = 343;
    %     boundary_collision_error_margin = 0.02;
    %
    % Joe.

    figureHandler = figure();
    vertices = [];
    while 1
        % Update the plot
        clf
        if ~isequal(vertices, [])
            plot(vertices(:,1), vertices(:,2), 'Color', [.8 .2 .3]);
        end
        xlabel('X-coordinate [m]');
        ylabel('Y-coordinate [m]');
        title(sprintf('Left-click -> add new vertex\nRight click -> Finish and close contour'));
        xlim([-5 5])
        ylim([-5 5])
        grid

        % Wait for the user mouse click
        [x, y, button] = ginput(1);
        if button == 3
            break;
        end

        % Add the new point to the vertex list
        vertices(size(vertices,1)+1, :) = [x, y];
    end

    clf(figureHandler);
    % Save it if there's information
    if ~isequal(vertices, [])

        % Create boundaries from captured vertices
        for i = [1 : size(vertices,1)-1]
            room_boundaries(i, :) = [vertices(i, :) vertices(i+1, :)];
        end
        room_boundaries(i+1, :) = [vertices(i+1, :) vertices(1, :)];

        boundaries_reflection_coefficients = zeros(size(room_boundaries, 1), 1);
        % Ask for boundaries coefficients
        for i = [1 : size(room_boundaries, 1)]
            clf
            for boundary = room_boundaries'
                plot([boundary(1) boundary(3)], [boundary(2) boundary(4)], 'Color', [.2 .8 .3])
                hold on
            end
            xlabel('X-coordinate [m]');
            ylabel('Y-coordinate [m]');
            xlim([-5 5])
            ylim([-5 5])
            grid

            % Show the current boundary on a marked color
            plot([room_boundaries(i, 1) room_boundaries(i, 3)], [room_boundaries(i, 2) room_boundaries(i, 4)], 'LineWidth', 2, 'Color', [.8 .2 .3])

            % Ask for user input
            boundaries_reflection_coefficients(i) = str2num(cell2mat(inputdlg(sprintf('Insert reflection coefficient for marked\nboundary, for example:  0.8'), 'Reflection coeffs.', 1, {'1.0'})));
        end

        clf
        for boundary = room_boundaries'
            plot([boundary(1) boundary(3)], [boundary(2) boundary(4)], 'Color', [.2 .8 .3])
            hold on
        end
        xlabel('X-coordinate [m]');
        ylabel('Y-coordinate [m]');
        xlim([-5 5])
        ylim([-5 5])
        grid

        % Ask for source and microphone positions
        title('Left click -> Place source');
        [x, y, button] = ginput(1);
        while button ~= 1
            [x, y, button] = ginput(1);
        end
        source_position = [x, y];
        plot(x, y, 'x', 'MarkerSize', 10, 'LineWidth', 3, 'Color', [.2 .3 .8]);

        title('Left click -> Place microphone');
        [x, y, button] = ginput(1);
        while button ~= 1
            [x, y, button] = ginput(1);
        end
        mic_position = [x, y];
        plot(x, y, 'o', 'MarkerSize', 10, 'LineWidth', 3, 'Color', [.8 .2 .4]);

        close(figureHandler);

        max_reflection_order = 100;
        mic_collision_error_margin = 0.05;
        scan_angle_min = -180;
        scan_angle_resolution = 0.5;
        scan_angle_max = 180;
        show_animated_plot = 0;
        include_real_source = 1;
        sample_rate = 48000;
        propagation_speed = 343;
        boundary_collision_error_margin = 0.02;

        % Save captured data to file
        output_filename = cell2mat(inputdlg(sprintf('Insert filename to save captured data to: '), 'Save data to file', 1, {'captured_data.mat'}));
        save(output_filename, ...
                             'source_position', ...
                             'mic_position', ...
                             'max_reflection_order', ...
                             'mic_collision_error_margin', ...
                             'room_boundaries', ...
                             'boundaries_reflection_coefficients', ...
                             'scan_angle_min', ...
                             'scan_angle_resolution', ...
                             'scan_angle_max', ...
                             'show_animated_plot', ...
                             'include_real_source', ...
                             'sample_rate', ...
                             'propagation_speed', ...
                             'boundary_collision_error_margin');
    end
end
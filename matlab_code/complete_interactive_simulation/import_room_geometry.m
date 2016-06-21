function room_boundaries = import_room_geometry(varargin)
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
    askForAbsorption = 0;

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
        xlim([-4 10])
        ylim([-4 10])
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
            if askForAbsorption
                clf
                for boundary = room_boundaries'
                    plot([boundary(1) boundary(3)], [boundary(2) boundary(4)], 'Color', [.2 .8 .3])
                    hold on
                end
                xlabel('X-coordinate [m]');
                ylabel('Y-coordinate [m]');
                xlim([-4 10])
                ylim([-4 10])
                grid

                % Show the current boundary on a marked color
                plot([room_boundaries(i, 1) room_boundaries(i, 3)], [room_boundaries(i, 2) room_boundaries(i, 4)], 'LineWidth', 2, 'Color', [.8 .2 .3])

                % Ask for user input
                boundaries_reflection_coefficients(i) = str2num(cell2mat(inputdlg(sprintf('Insert reflection coefficient for marked\nboundary, for example:  0.8'), 'Reflection coeffs.', 1, {'1.0'})));
                
            else
                boundaries_reflection_coefficients(i) = 1;
            end
        end

        clf
        for boundary = room_boundaries'
            plot([boundary(1) boundary(3)], [boundary(2) boundary(4)], 'Color', [.2 .8 .3])
            hold on
        end
        xlabel('X-coordinate [m]');
        ylabel('Y-coordinate [m]');
        xlim([-4 10])
        ylim([-4 10])
        grid

        % Ask for array position
        title('Left click -> Place microphones array');
        [x, y, button] = ginput(1);
        while button ~= 1
            [x, y, button] = ginput(1);
        end
        mic_position = [x, y];
        plot(x, y, 'o', 'MarkerSize', 10, 'LineWidth', 3, 'Color', [.8 .2 .4]);

        % Ask for source positions
        title(sprintf('Left click -> Place source positions\nRight click -> No more source positions'));
        source_positions = [];
        while 1
            [x, y, button] = ginput(1);
            if button == 1
                source_positions(size(source_positions, 1) + 1, :) = [x, y];
            else if button == 3
                    if size(source_positions, 1) > 0
                        break
                    end
                end
            end
            for index = 1 : size(source_positions, 1)
                plot(source_positions(index, 1), source_positions(index, 2), 'x', 'MarkerSize', 10, 'LineWidth', 3, 'Color', [.2 .3 .8]);
            end
        end

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
        if nargin == 0
            output_filename = cell2mat(inputdlg(sprintf('Insert filename to save captured data to: '), 'Save data to file', 1, {'captured_data.mat'}));
        else
            output_filename = char(varargin);
        end
        save(output_filename, ...
                             'source_positions', ...
                             'mic_position', ...
                             'room_boundaries', ...
                             'boundaries_reflection_coefficients');
    end
end
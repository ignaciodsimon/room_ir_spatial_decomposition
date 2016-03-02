function display_geometry_from_file(geometry_filename)

    % Load data from input file
    loaded_data = load(geometry_filename);
    mic_position = loaded_data.mic_position;
    source_position = loaded_data.source_position;
    room_boundaries = loaded_data.room_boundaries;

    figureHandler = figure();

    % Plot microphone and source
    plot(mic_position(1), mic_position(2), 'o', 'MarkerSize', 10, 'LineWidth', 2, 'Color', [.1 .7 .2])
    hold on
    plot(source_position(1), source_position(2), '^', 'MarkerSize', 10, 'LineWidth', 2, 'Color', [.7 .1 .2])

    % Plot boundaries
    for i = 1 : size(room_boundaries, 1);
        plot([room_boundaries(i, 1) room_boundaries(i, 3)], [room_boundaries(i, 2), room_boundaries(i, 4)], 'LineWidth', 2, 'Color', [.1 .2 .7])
    end

    % Format plot
    grid
    previous_x_lim = get(gca, 'xlim');
    xlim([previous_x_lim(1)-1 previous_x_lim(2)+1])
    previous_y_lim = get(gca, 'ylim');
    ylim([previous_y_lim(1)-1 previous_y_lim(2)+1])
    title_handler = title(sprintf('Geometry from file "%s"', geometry_filename));
    set(title_handler, 'interpreter', 'none')
    legend({'Microphone', 'Source'}, 'Location', 'best')
    xlabel('X-Coordinate [m]')
    ylabel('Y-Coordinate [m]')

    return

    % Save plot to PDF
    set(gcf, 'PaperPosition', [-0.42 +0.05 9.3 3.6]);
    set(gcf, 'PaperSize', [8.1 3.5]);
    saveas(gcf, 'displayed_geometry.pdf', 'pdf');
    close(figureHandler);

end

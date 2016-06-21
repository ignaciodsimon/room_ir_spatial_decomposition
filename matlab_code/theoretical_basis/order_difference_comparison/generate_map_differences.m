function generate_map_differences()
    % Uses three maps generated with different virtual image orders to
    % generate two difference-maps, showing that after a certain order,
    % considering additional image sources does not introduce much change
    % to the final reconstructred sound field.
    %
    % Joe.

    % Load input files
    map1_filename = 'power_map_maxorder_1_freq_200.mat';
    map2_filename = 'power_map_maxorder_5_freq_200.mat';
    map3_filename = 'power_map_maxorder_10_freq_200.mat';
    map1 = load(map1_filename);
    map2 = load(map2_filename);
    map3 = load(map3_filename);

    % Compute difference maps
    difference_1 = map2.amplitudeMap - map1.amplitudeMap;
    difference_2 = map3.amplitudeMap - map2.amplitudeMap;

    % -------------------- MAP LEFT --------------------

    subplot(1,2,1);
    % Power map
    imagesc(map1.x_coordinates, map1.y_coordinates, difference_1);
    set(gca,'YDir','normal')
    hold on
    % Crosses on top of sources
    for i = 1 : size(map1.inputSignals_coordinates, 1)

        x_pos = map1.inputSignals_coordinates(i, 1);
        y_pos = map1.inputSignals_coordinates(i, 2);

        plot(x_pos, y_pos, 'x', 'MarkerSize', 20, 'LineWidth', 6, 'Color', [1 1 1]);
        plot(x_pos, y_pos, 'x', 'MarkerSize', 10, 'LineWidth', 4, 'Color', [1 .2 .3]);
    end
    % Rest of the furniture of the plot
    xlabel('X axis [m]', 'FontSize', 14);
    ylabel('Y axis [m]', 'FontSize', 14);
    caxis([-30 0])
    colorBarHandler = colorbar('eastoutside', 'FontSize', 14);
    colorBarHandler.Label.String = 'Power [dB]';
    title('Comparison orders 1-5', 'FontSize', 14);

    % -------------------- MAP RIGHT --------------------

    subplot(1,2,2);
    % Power map
    imagesc(map1.x_coordinates, map1.y_coordinates, difference_2);
    set(gca,'YDir','normal')
    hold on
    % Crosses on top of sources
    for i = 1 : size(map1.inputSignals_coordinates, 1)

        x_pos = map1.inputSignals_coordinates(i, 1);
        y_pos = map1.inputSignals_coordinates(i, 2);

        plot(x_pos, y_pos, 'x', 'MarkerSize', 20, 'LineWidth', 6, 'Color', [1 1 1]);
        plot(x_pos, y_pos, 'x', 'MarkerSize', 10, 'LineWidth', 4, 'Color', [1 .2 .3]);
    end
    % Rest of the furniture of the plot
    xlabel('X axis [m]', 'FontSize', 14);
    ylabel('Y axis [m]', 'FontSize', 14);
    caxis([-30 0])
    colorBarHandler = colorbar('eastoutside', 'FontSize', 14);
    colorBarHandler.Label.String = 'Power [dB]';
    title('Comparison orders 5-10', 'FontSize', 14);

    % Save plot to file
    set(gcf, 'PaperPosition', [-0.5 0.10 9.4 3.25]);
    set(gcf, 'PaperSize', [8.5 3.5]);
    saveas(gcf, 'order_difference_comparison.pdf', 'pdf');
    close all

end

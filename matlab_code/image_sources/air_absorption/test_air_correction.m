function test_air_correction()
    % Script used to compute the air absorption filter for several
    % distances and plot them to a file, to show that with the distance not
    % only the overal level is reduced, but also the frequency response
    % changes.
    %
    % Joe.

    % Get the filters IRs by sending a delta as an input
    ir_500 = correct_air_absorption([1], 500);
    ir_1000 = correct_air_absorption([1], 1000);
    ir_2000 = correct_air_absorption([1], 2000);
    ir_4000 = correct_air_absorption([1], 4000);
    ir_8000 = correct_air_absorption([1], 8000);
    time_axis = [1 : length(ir_500)] / 48000;

    % Create plot with different filters
    figureHandler = figure();
    plot(time_axis, ir_500, 'LineWidth', 2)
    hold on
    plot(time_axis, ir_1000, 'LineWidth', 2)
    plot(time_axis, ir_2000, 'LineWidth', 2)
    plot(time_axis, ir_4000, 'LineWidth', 2)
    plot(time_axis, ir_8000, 'LineWidth', 2)
    grid
    legend({'d = 500 [m]', 'd = 1 [km]', 'd = 2 [km]', 'd = 4 [km]', 'd = 8 [km]'}, 'FontSize', 12)
    xlabel('Time [s]', 'FontSize', 12)
    ylabel('Amplitude [.]', 'FontSize', 12)
    set(gca, 'FontSize', 11);
    title(sprintf('Impulse response of atmospheric absorption filter\n(derivated from the values of ISO9613-2)'), 'FontSize', 12)

    % Saves plot to PDF
    set(gcf, 'PaperPosition', [-0.42 +0.05 9.1 3.4]);
    set(gcf, 'PaperSize', [8.1 3.5]);
    saveas(gcf, 'air_correction_filters.pdf', 'pdf');
%     close(figureHandler);

end

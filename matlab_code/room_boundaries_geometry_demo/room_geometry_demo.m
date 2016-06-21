function room_geometry_demo()

    boundaries = [0,0,0,8;0,8,7,8;7,8,7,7;7,7,7,6;7,6,7,4;7,4,8,4;8,4,8,1;8,1,7,1;7,1,7,0;7,0,6,0;6,0,6,0.250000000000000;6,0.250000000000000,4,0.250000000000000;4,0.250000000000000,4,0;4,0,3,0;3,0,3,0.250000000000000;3,0.250000000000000,1,0.250000000000000;1,0.250000000000000,1,0;1,0,0,0];

    for i = 1 : size(boundaries, 1)
        plot([boundaries(i, 1) boundaries(i, 3)], [boundaries(i, 2) boundaries(i, 4)], 'Color', 'blue', 'LineWidth', 1.5);
        hold on
    end
    xlim([-0.5 8.5])
    ylim([-0.5 8.5])
    grid
    xlabel('X-coordinate [m]')
    ylabel('Y-coordinate [m]')
    title('Example of a polygon-constructed enclosure boundaries')

    % Save plot to file
    set(gcf, 'PaperPosition', [-0.5 -0.4 9.4 8.25]);
    set(gcf, 'PaperSize', [8.5 7.5]);
    saveas(gcf, 'room_boundaries_demo.pdf', 'pdf');
    close all
end

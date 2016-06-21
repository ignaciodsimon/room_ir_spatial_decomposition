function comparePrecision()


    loadedMap = load('pentagonal_room_averaged_power_map.mat');
    loadedBoundaries = load('pentagonal_room_clicked_points.mat');

    % Find the maximum deviation of boundaries
    maxDeviation = max(max(abs(loadedMap.ROOM_BOUNDARIES - loadedBoundaries.estimatedBoundaries)));

    % Create plot
    imagesc(loadedMap.x_coordinates, loadedMap.y_coordinates, loadedMap.averagedPowerMap);
    hold on
    set(gca,'YDir','normal')
    xlabel('X-coordinate [m]')
    ylabel('Y-coordinate [m]')
    plot(loadedMap.ARRAY_CENTER(1), loadedMap.ARRAY_CENTER(2), 'o', 'MarkerSize', 10, 'Color', 'red', 'LineWidth', 3);
%     for index = 1 : size(loadedBoundaries.estimatedBoundaries)
%         plot([loadedBoundaries.estimatedBoundaries(index, 1) loadedBoundaries.estimatedBoundaries(index, 3)], ...
%              [loadedBoundaries.estimatedBoundaries(index, 2) loadedBoundaries.estimatedBoundaries(index, 4)], '--r', 'LineWidth', 4);
%         plot([loadedMap.ROOM_BOUNDARIES(index, 1) loadedMap.ROOM_BOUNDARIES(index, 3)], ...
%              [loadedMap.ROOM_BOUNDARIES(index, 2) loadedMap.ROOM_BOUNDARIES(index, 4)], ':g', 'LineWidth', 2);
%     end
    title(sprintf('Max dev.: %.3f [m]', maxDeviation));

    % Save plot to file
    set(gcf, 'PaperPosition', [-0.25 -0.25 5.7 7.5]);
    set(gcf, 'PaperSize', [5.0 7.0]);
    saveas(gcf, 'pentagonal_room_NO_boundaries_example.pdf', 'pdf');
    close all
% return

    
    loadedMap = load('square_room_averaged_power_map.mat');
    loadedBoundaries = load('square_room_clicked_points.mat');

    % Find the maximum deviation of boundaries
    maxDeviation = max(max(abs(loadedMap.ROOM_BOUNDARIES - loadedBoundaries.estimatedBoundaries)));

    % Create plot
    imagesc(loadedMap.x_coordinates, loadedMap.y_coordinates, loadedMap.averagedPowerMap);
    hold on
    set(gca,'YDir','normal')
    xlabel('X-coordinate [m]')
    ylabel('Y-coordinate [m]')
    plot(loadedMap.ARRAY_CENTER(1), loadedMap.ARRAY_CENTER(2), 'o', 'MarkerSize', 10, 'Color', 'red', 'LineWidth', 3);
%     for index = 1 : size(loadedBoundaries.estimatedBoundaries)
%         plot([loadedBoundaries.estimatedBoundaries(index, 1) loadedBoundaries.estimatedBoundaries(index, 3)], ...
%              [loadedBoundaries.estimatedBoundaries(index, 2) loadedBoundaries.estimatedBoundaries(index, 4)], '--r', 'LineWidth', 4);
%         plot([loadedMap.ROOM_BOUNDARIES(index, 1) loadedMap.ROOM_BOUNDARIES(index, 3)], ...
%              [loadedMap.ROOM_BOUNDARIES(index, 2) loadedMap.ROOM_BOUNDARIES(index, 4)], ':g', 'LineWidth', 2);
%     end
    title(sprintf('Max dev.: %.3f [m]', maxDeviation));

    % Save plot to file
    set(gcf, 'PaperPosition', [-0.25 -0.25 5.7 7.5]);
    set(gcf, 'PaperSize', [5.0 7.0]);
    saveas(gcf, 'square_room_NO_boundaries_example.pdf', 'pdf');
    close all




    % Show averaged image with an overlay of the real boundaries
    
    % Measure the maximum offset of all boundaries

    



end

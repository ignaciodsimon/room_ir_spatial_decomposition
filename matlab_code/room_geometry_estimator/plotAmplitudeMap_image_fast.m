function plotAmplitudeMap_image_fast(filenames, boundaries, signalFreqs, propagationSpeed, maxVirtualSourceOrder)
    % Displays the amplitude / power map generated using the function "generatePowerMapFast.m"
    % It includes the boundaries (if any) and array elements in the plot.
    %
    % Joe.

    figureHandler = figure('units', 'normalized', 'position', [.1 .1 .58 .85]);

    loadedData = load(filenames{1});

    averagedData = zeros(size(loadedData.amplitudeMap));
    for i = 1 : length(filenames)
        loadedData = load(filenames{i});
        averagedData = averagedData + loadedData.amplitudeMap;
    end
    averagedData = averagedData / length(filenames);

    averagedData = averagedData - max(max(averagedData));

    minimumThreshold = -70;
    maximumThreshold = -10;
    if min(min(averagedData)) >= maximumThreshold
        maximumThreshold = 0;
    end
    for index1 = 1 : size(averagedData, 1)
        for index2 = 1 : size(averagedData, 2)
            if averagedData(index1, index2) < minimumThreshold
                averagedData(index1, index2) = minimumThreshold;
            else
                if averagedData(index1, index2) > maximumThreshold
                    averagedData(index1, index2) = maximumThreshold;
                end
            end
        end
    end
    averagedData = averagedData - max(max(averagedData));

    % Power map
    imagesc(loadedData.x_coordinates, loadedData.y_coordinates, averagedData);
    set(gca,'YDir','normal')
    hold on

    % Crosses on top of sources
    for i = 1 : size(loadedData.inputSignals_coordinates, 1)

        x_pos = loadedData.inputSignals_coordinates(i, 1);
        y_pos = loadedData.inputSignals_coordinates(i, 2);

        plot(x_pos, y_pos, 'x', 'MarkerSize', 20, 'LineWidth', 6, 'Color', [1 1 1]);
        plot(x_pos, y_pos, 'x', 'MarkerSize', 10, 'LineWidth', 4, 'Color', [1 .2 .3]);
    end

    % Boundaries 
    for i = 1 : size(boundaries, 1)

%         plot([boundaries(i, 1) boundaries(i, 3)], ...
%               [boundaries(i, 2) boundaries(i, 4)], ...
%               'LineWidth', 6, 'Color', [1 1 1]);
% 
        plot([boundaries(i, 1) boundaries(i, 3)], ...
              [boundaries(i, 2) boundaries(i, 4)], ...
              'LineWidth', 4, 'Color', [1 .2 .3]);


%         plot([boundaries(i, 1) boundaries(i, 3)], ...
%               [boundaries(i, 2) boundaries(i, 4)], ...
%               '.', 'LineWidth', 1, 'Color', [1 .2 .3]);
    end

    % Rest of the furniture of the plot
    xlabel('X axis [m]', 'FontSize', 14);
    ylabel('Y axis [m]', 'FontSize', 14);
%     caxis([-40 25])
%     caxis([-60 0])

%     set(figureHandler, 'Position', [500 100 1200 800])

    colorBarHandler = colorbar('eastoutside', 'FontSize', 14);
    colorBarHandler.Label.String = 'Power [dB]';

    if maxVirtualSourceOrder >= 0
        maxVirtualSourceOrderString = sprintf(' - Max. order: %d', maxVirtualSourceOrder);
    else
        maxVirtualSourceOrderString = '';
    end
    
    if length(signalFreqs) > 1
        title(sprintf('Freq: (%.0f to %.0f) [Hz] - Propagation speed: %.1f [m/s]%s', ...
            signalFreqs(1), signalFreqs(length(signalFreqs)),...
            propagationSpeed, ...
            maxVirtualSourceOrderString), ...
            'FontSize', 14);
    else
        title(sprintf('Freq: %.0f [Hz] - Propagation speed: %.1f [m/s]%s', ...
            signalFreqs, ...
            propagationSpeed, ...
            maxVirtualSourceOrderString), ...
            'FontSize', 14);
    end

    % Save plot to file
    set(gcf, 'PaperPosition', [-0.5 -0.35 9.4 8.05]);
    set(gcf, 'PaperSize', [8.5 7.5]);
    saveas(gcf, sprintf('%s.pdf', filenames{1}), 'pdf');

%     close all
end

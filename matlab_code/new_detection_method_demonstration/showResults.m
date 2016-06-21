function showResults()

    disp(' >> Loading data and generating plot ...')
    loadedResults = load('results_data.mat');
    results = loadedResults.results;

    lastSNRValue = inf;
    currentCurve = [];
    snrValues = [];
    for index = 1 : size(results, 1)
        if results(index, 1) ~= lastSNRValue
            lastSNRValue = results(index, 1);
            snrValues(size(snrValues, 1) + 1, :) = lastSNRValue;

            % Plot the new curve
            if ~isempty(currentCurve)
                plot(currentCurve(:, 1) * 1000, currentCurve(:, 2), 'LineWidth', 2);
                hold on
                currentCurve = [];
            end
        else
            % Keep adding points to the current curve
            currentCurve(size(currentCurve, 1) + 1, :) = results(index, 2 : 3);
        end
    end
    % Plot the new curve
    if ~isempty(currentCurve)
        plot(currentCurve(:, 1) * 1000, currentCurve(:, 2), 'LineWidth', 2);
        hold on
        currentCurve = [];
    end

    % Generate legend texts
    legendEntries = cell(length(snrValues), 1);
    for index = 1 : length(snrValues)
        legendEntries{index} = sprintf('SNR: %.0f [dB]', snrValues(index));
    end
    legend(legendEntries);

    % General plot furniture
    grid on
    xlabel('Radius difference between the two sources [mm]')
    ylabel('Needed DOA difference to detect the two sources [degrees]')

    % Save plot to file
    set(gcf, 'PaperPosition', [-0.5 -0.2 9.4 6.05]);
    set(gcf, 'PaperSize', [8.5 5.5]);
    saveas(gcf, 'new_detection_resolution.pdf', 'pdf');
    close all
    disp('All done!')

end

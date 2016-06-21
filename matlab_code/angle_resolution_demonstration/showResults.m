function showResults()

    loadedData = load('results.mat');

    firstElementToDisplay = 1;
    amountToDisplay = 10;

%     for index = firstElementToDisplay : amountToDisplay
%         plot(loadedData.SNR_VALUES, loadedData.results(index, :), 'o', 'LineWidth', 1.5);
%         hold on
%     end
    interpolatedXaxis = [-2 : 1 : 25];
    for index = firstElementToDisplay : amountToDisplay
%         interpolatedYaxis = interp1(loadedData.SNR_VALUES, loadedData.results(index, :), interpolatedXaxis, 'pchip');
%         plot(interpolatedXaxis, interpolatedYaxis, 'x-', 'LineWidth', 1.5);
        plot(loadedData.SNR_VALUES, loadedData.results(index, :), 'x-', 'LineWidth', 1.5);
        hold on
    end
    grid on

    legendStrings = {};
    for index = firstElementToDisplay : amountToDisplay
        legendStrings{index - firstElementToDisplay + 1} = sprintf('%d sensors', loadedData.AMOUNT_OF_SENSORS(index));
    end

    ylim([0 105]);
    legend(legendStrings, 'FontSize', 12, 'Location', 'best');
    xlabel('SNR [dB]', 'FontSize', 12);
    ylabel('Successful DOA estimations [%]', 'FontSize', 12);

    % Save plot to file
    set(gcf, 'PaperPosition', [-0.5 0.05 9.7 3.05]);
    set(gcf, 'PaperSize', [8.5 3.0]);
    saveas(gcf, 'snr_vs_doa.pdf', 'pdf');
    close all

end

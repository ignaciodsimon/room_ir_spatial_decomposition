function combineAbsorptionPlots()

    data1 = load('estimated_absorption_window_0.0010.mat');
    data2 = load('estimated_absorption_window_0.0015.mat');
    data3 = load('estimated_absorption_window_0.0020.mat');

    semilogx([1 : 24000/length(data1.materialCharacteristics) : 24000], data1.materialCharacteristics, 'Color', [0.2 0.7 0.1], 'LineWidth', 4);
    hold on
    semilogx([1 : 24000/length(data1.estimatedAbsorption) : 24000], data1.estimatedAbsorption, '--', 'LineWidth', 2.5)
    semilogx([1 : 24000/length(data2.estimatedAbsorption) : 24000], data2.estimatedAbsorption, '--', 'LineWidth', 2.5)
    semilogx([1 : 24000/length(data3.estimatedAbsorption) : 24000], data3.estimatedAbsorption, '--', 'LineWidth', 2.5)
    grid on

    xlabel('Frequency [Hz]')
    ylabel('Absorption [.]')
    xlim([63 12000])
    ylim([0.5 1])

    legend({'Material characteristics', 'Estimated absorption (window 1.0 ms)', 'Estimated absorption (window 1.5 ms)', 'Estimated absorption (window 2.0 ms)'}, 'Location', 'best')
    % Save plot to file
    set(gcf, 'PaperPosition', [-0.5 0.05 9.7 4.05]);
    set(gcf, 'PaperSize', [8.5 4.0]);
    saveas(gcf, 'combined_absorption_estimations.pdf', 'pdf');
    close all
end

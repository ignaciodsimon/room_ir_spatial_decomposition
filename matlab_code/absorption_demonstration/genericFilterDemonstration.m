function genericFilterDemonstration()

    BOUNDARY_ABSORPTION = [ 125  250  500 1000 2000 4000 8000 16000 24000
                           0.20 0.22 0.30 0.40 0.49 0.56 0.600 0.63  0.65];

    % Plot of material characteristics
    semilogx(BOUNDARY_ABSORPTION(1, 1 : 6), BOUNDARY_ABSORPTION(2, 1 : 6), 'd-', 'LineWidth', 2, 'MarkerSize', 8)
    grid on
    xlim([20 20000])
    ylim([0 1])
    xlabel('Frequency [Hz]')
    ylabel('Absorption coefficient [.]')
    set(gcf, 'PaperPosition', [-0.5 0.125 9.7 3.0]);
    set(gcf, 'PaperSize', [8.5 3.0]);
    saveas(gcf, 'absorption_demonstration_material.pdf', 'pdf');
    close all

    inputSignal = [1 zeros(1, 47999)];
    [filteredSignal, filterGroupDelay] = generic_filter_design(inputSignal, ...
                                                               BOUNDARY_ABSORPTION(1,:), ...
                                                               ones(1, 9) - (BOUNDARY_ABSORPTION(2,:)));

    plot(filteredSignal, 'LineWidth', 2)
    grid on
    xlim([0 1200])
    ylim([-0.05 0.5])
    xlabel('Time [samples @ 48 kHz]')
    ylabel('Amplitude [.]')
    set(gcf, 'PaperPosition', [-0.5 0.125 9.7 3.0]);
    set(gcf, 'PaperSize', [8.5 3.0]);
    saveas(gcf, 'absorption_demonstration_ir.pdf', 'pdf');
    close all

    outputSpectrum = 20*log10(abs(fft(filteredSignal)));
    outputSpectrum = outputSpectrum(1 : round(length(outputSpectrum)/2));
    semilogx(outputSpectrum, 'LineWidth', 2);
    hold on
    semilogx(BOUNDARY_ABSORPTION(1, 1 : 6), 20*log10((ones(1, 6) - BOUNDARY_ABSORPTION(2, 1 : 6))) -0.65, 'd', 'MarkerSize', 8, 'LineWidth', 2)
    grid on
    ylim([-12 0])
    xlim([20 20000])
    legend({'Obtained filter frequency response', 'Material absorption characteristics'})
    xlabel('Frequency [Hz]')
    ylabel('Frequency response module [dB]')
    set(gcf, 'PaperPosition', [-0.5 0.125 9.7 3.0]);
    set(gcf, 'PaperSize', [8.5 3.0]);
    saveas(gcf, 'absorption_demonstration_filter_response.pdf', 'pdf');
    close all
    

end

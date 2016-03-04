function display_results_freq_response()

    % Load IR
    loaded_data = load('simulation_results.mat');

    % Spectrum of loaded IR
    spectrum = 20*log10(abs(fft(loaded_data.convolved_ir, 48000)));

    % Produce an average that uses a bigger window as frequency grows, for
    % easier display on the plot
    evolving_average = zeros(size(spectrum));
    for i = 1 : length(spectrum)
        if (i+i) > length(spectrum)
            break
        end
        evolving_average(i) = mean(spectrum(i - round(i*0.1): i + round(i*0.1)));
    end
    evolving_average = evolving_average(1 : i-1);

    % Show results
    semilogx(spectrum(1:length(spectrum)/2) -3);
    hold on
    semilogx(evolving_average -3, 'LineWidth', 2)
    xlim([0 24000])
    ylim([-45 15])
    grid on
    xlabel('Frequency [Hz]')
    ylabel('Amplitude [dB rel. 1 kHz]')
    legend({'Spectrum of IR (module)', 'Averaged display'}, 'Location', 'best', 'FontSize', 12)

end

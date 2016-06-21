function evolving_average = averageSpectrum(inputSpectrum)

    % Produce an average that uses a bigger window as frequency grows, for
    % easier display on the plot
    evolving_average = zeros(size(inputSpectrum));
    for i = 1 : length(inputSpectrum)
        if (i+i) > length(inputSpectrum)
            break
        end
        evolving_average(i) = mean(inputSpectrum(i - round(i * 0.1): i + round(i * 0.1)));
    end
    evolving_average = evolving_average(1 : i-1);

    % Show results
%     semilogx(inputSpectrum(1:length(inputSpectrum)/2) -3);
%     hold on
%     semilogx(evolving_average -3, 'LineWidth', 2)
%     xlim([0 24000])
%     ylim([-45 15])
%     grid on
%     xlabel('Frequency [Hz]')
%     ylabel('Amplitude [dB rel. 1 kHz]')
%     legend({'Spectrum of IR (module)', 'Averaged display'}, 'Location', 'best', 'FontSize', 12)


end

function testFilter()

    BOUNDARY_ABSORPTION = [ 125  250  500 1000 2000 4000 8000 16000 24000
%                            0.99 0.99 0.99 0.99 0.99 0.99];
                           0.20 0.30 0.40 0.50 0.60 0.70 0.8 0.9 0.99];

    inputSignal = [zeros(1,501) 1 zeros(1, 1000)];

    [filteredSignal, filterGroupDelay] = generic_filter_design(inputSignal, ...
                                                               BOUNDARY_ABSORPTION(1,:), ...
                                                               ones(1, size(BOUNDARY_ABSORPTION, 2)) - BOUNDARY_ABSORPTION(2,:));
    filteredSignal = filteredSignal(filterGroupDelay : length(filteredSignal));
    filteredSignal = pad_with_zeros(filteredSignal, length(inputSignal));

    subplot(2,2,1)
    plot(inputSignal)
    grid on
    title('Input signal')
    
    subplot(2,2,3)
    plot(filteredSignal)
    grid on
    title('Filtered signal')

    subplot(2,2,2)
    semilogx(20*log10(abs(fft(inputSignal))));
    ylim([-60 20])
    grid on
    title('Spectrum input signal')
    
    subplot(2,2,4)
    semilogx(20*log10(abs(fft(filteredSignal,48000))));
    hold on
    semilogx(BOUNDARY_ABSORPTION(1,:), 20*log10(ones(1, size(BOUNDARY_ABSORPTION, 2)) - BOUNDARY_ABSORPTION(2, :)), 'x');
    ylim([-60 20])
    grid on
    title('Spectrum filtered signal vs desired filter shape')
    legend({'Spectrum filtered signal', 'Desired filter points'})
    
end

function padded_vector = pad_with_zeros(input_vector, desired_length)

    if length(input_vector) == desired_length
        padded_vector = input_vector;
        return
    end

    if length(input_vector) > desired_length
        padded_vector = input_vector(1 : desired_length);
    end
    
    if length(input_vector) < desired_length
        padded_vector = zeros(1, desired_length);
        padded_vector(1 : length(input_vector)) = input_vector;
    end

end

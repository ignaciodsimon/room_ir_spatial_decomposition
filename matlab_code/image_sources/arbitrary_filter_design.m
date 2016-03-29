function filtered_signal = arbitrary_filter_design(input_signal, bands_frequencies, bands_levels)
    % Creates a linear-phase filter based on the bands level given and uses
    % it to filter input signal. It tries to produce a causal filter and
    % compensates for the group delay.
    %
    % Note: The input signal is assumed to be sampled at 48 kHz.
    %
    % Joe.

    SAMPLE_FREQ = 48000;

    % Interpolate the input bands levels to obtain a smooth response over
    % the complete spectrum
    freqBins = [1 : SAMPLE_FREQ/2];
    interpolated_bands_levels = zeros(length(freqBins), 2);
    interpolated_bands_levels(:, 1) = freqBins;
    interpolated_bands_levels(:, 2) = spline(bands_frequencies, bands_levels, freqBins);

    % Obtain the module of the frequency response from the given band levels
    filter_module = [interpolated_bands_levels(1,2)    % Consider the first entry also the value for DC
                     interpolated_bands_levels(:,2)    % First part of the spectrum
                                                       % Mirrored part of the spectrum
                     flipud(interpolated_bands_levels(1 : length(interpolated_bands_levels)-1, 2))];

    % Insert by default a big group delay that will assure the IR to be
    % causal. It will be removed afterwards.
    group_delay_samples = 10000;

    % Design the phase of the frequency response to be linear
    filter_phase = [1 : length(filter_module)/2]' / (length(filter_module)/2) * group_delay_samples * -pi;
    filter_phase = [filter_phase
                    filter_phase];

    % Obtain the filter IR using the inverse FFT
    filter_response = filter_module .* (cos(filter_phase) + 1i*sin(filter_phase));
    filter_ir = real(ifft(filter_response))';
    filter_ir = filter_ir(1 : group_delay_samples * 2);    

    % Filter input signal with generated IR
    filtered_signal = conv(input_signal, filter_ir);

    % Remove the introduced group delay
    filtered_signal = filtered_signal(group_delay_samples + 1: length(input_signal) + group_delay_samples -1 +1);

%     close all
%     subplot(3,1,1)
%     plot(input_signal)
%     grid
%     subplot(3,1,2)
%     plot(filtered_signal)
%     grid
%     subplot(3,1,3)
%     plot(filter_ir)
%     grid

end
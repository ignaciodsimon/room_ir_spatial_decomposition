function filtering_tests()

%     sinc_function = sinc([-100 : 0.01 : 100]);
% 
% %     semilogx(20*log10(abs(fft(sinc_function))))
%     plot(unwrap(angle(fft(sinc_function))))
%     grid
%     return






























%     data_x = [125 250 500 1000 2000 4000];
%     data_y = [  1   1   1  0.05  0.05  0.005];
% 
%     x_values = [1 : 100: 24000];
%     y_values = interp1(data_x, data_y, x_values, 'pchip', 'extrap');
% 
% 
% %     x_values_2 = [1 : 1: 24000];
% %     y_values_2 = interp1(x_values, y_values, x_values_2, 'spline');
%     
%     
%     
% %     subplot(2,1,1)
%     semilogx(data_x, data_y, 'x-')
%     ylim([0 1.5])
%     grid
% %     subplot(2,1,2)
%     hold on
%     semilogx(x_values, y_values, 'x-')
% %     semilogx(x_values_2, y_values_2 - 0.25, '-')
%     ylim([0 1.5])
% %     grid
% 
% 
% %       'linear'   - (default) linear interpolation
% %       'nearest'  - nearest neighbor interpolation
% %       'next'     - next neighbor interpolation
% %       'previous' - previous neighbor interpolation
% %       'spline'   - piecewise cubic spline interpolation (SPLINE)
% %       'pchip'    - shape-preserving piecewise cubic interpolation
% %       'cubic'    - same as 'pchip'
% %       'v5cubic'  - the cubic interpolation from MATLAB 5, which does not
% %                    extrapolate and uses 'spline' if X is not equally
% %                    spaced.
% 
% 
% 
% 
% 
% 
% 
% 
%     return

    input_signal = [1 zeros(1, 48000)];

    freq_bands = [125 250 500 1000 2000 4000];
    bands_levels = [1 1 1 0.1 0.1 0.1];

    filtered_signal = arbitrary_filter_design(input_signal, freq_bands, bands_levels);
%     filtered_signal(3001) = filtered_signal(3000);

    spectrum_input = 20*log10(abs(fft(input_signal)));
    spectrum_output = 20*log10(abs(fft(filtered_signal)));

    subplot(2,3,1)
    plot(input_signal)
    grid
    subplot(2,3,4)
    plot(filtered_signal)
    grid

    subplot(2,3,3)
    plot(angle(fft(input_signal)))
    grid
    subplot(2,3,6)
    plot(unwrap(angle(fft(filtered_signal))))
    grid


    subplot(2,3,2)
    semilogx(spectrum_input)
    grid
    ylim([-20 20])
    subplot(2,3,5)
    semilogx(spectrum_output)
    grid
%     ylim([-20 20])
    

end

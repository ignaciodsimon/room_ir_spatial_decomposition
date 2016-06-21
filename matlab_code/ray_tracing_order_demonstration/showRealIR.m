function showRealIR()

    [real_ir, SAMPLE_RATE] = audioread('meas_Aw_1_5s_02.wav');
    real_ir = real_ir(1300 : length(real_ir));
    plot([1 : 1 : length(real_ir)] / SAMPLE_RATE * 1000, real_ir); grid on
    xlim([0 1000])
    xlabel('Time [ms]')
    ylabel('Amplitude [.]')
    set(gcf, 'PaperPosition', [-0.5 +0.1 9.7 3.0]);
    set(gcf, 'PaperSize', [8.5 3.0]);
    saveas(gcf, 'real_ir.pdf', 'pdf');
    close all

end

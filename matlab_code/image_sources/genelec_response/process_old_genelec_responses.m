function process_old_genelec_responses()
    % Uses some old measured impulse responses from a Genelec loudspeaker
    % on the anechoic chamber to produce an averaged and trimmed one. This
    % should not be taken as a reference of the loudspeaker, but as a guide
    % to the expected behavior. It's only used for rough estimation
    % purposes.
    %
    % Joe.

    % Load the old measured IRs from the Genelec on the anechoic chamber
    disp('> Loading Genelec impulse responses ...');
    loaded_data = load('old_genelec_IRs.mat');
    ir1 = loaded_data.ir1;
    ir2 = loaded_data.ir2;
    ir3 = loaded_data.ir3;
    ir4 = loaded_data.ir4;
    ir5 = loaded_data.ir5;

    % Trim IRs to a shorter length and remove most of the group delay
    disp('> Averaging and trimming ...');
    initial_trim = 450 + 30;
    final_trim = 2000;
    ir1 = ir1(1 + initial_trim : 1 + initial_trim + final_trim);
    ir2 = ir2(4 + initial_trim : 4 + initial_trim + final_trim);
    ir3 = ir3(9 + initial_trim : 9 + initial_trim + final_trim);
    ir4 = ir4(1 + initial_trim - 14: 1 + initial_trim -14 + final_trim);
    ir5 = ir5(1 + initial_trim - 10: 1 + initial_trim -10 + final_trim);

    % Average *all* IRs (number 3 is left out cause it seems not well
    % matched on time)
    averaged_ir = (ir1 + ir2 + ir4 + ir5) / 4;
    averaged_ir = averaged_ir / max(abs(averaged_ir));

    subplot(2,1,1)
    plot(averaged_ir)
    xlabel('Time [samples @ 48 kHz]')
    ylabel('Amplitude [.]')
    xlim([-10 1000])
    grid on
    title('IR from Genelec 1031A loudspeaker')
    subplot(2,1,2)
    module = abs(fft(averaged_ir, 48000));
    semilogx(20*log10(module) -7.44)
    xlim([0 24000])
    xlabel('Frequency [Hz]')
    ylabel('Amplitude [dB rel. 1 kHz]')
    grid on
    title('Associated frequency response (module)')

    disp('> Saving to file ...');
    ir = averaged_ir;
    save('genelec_IR_processed.mat', 'ir');
end

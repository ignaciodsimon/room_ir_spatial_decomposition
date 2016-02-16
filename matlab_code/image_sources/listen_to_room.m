function listen_to_room(processed_data_filename)

    disp('Loading processed data ...');
    loaded_data = load(processed_data_filename);
    disp('Loading input audio ...');
    [input_audio, sample_rate] = audioread('bakunin.wav');
    if sample_rate ~= 48000
        disp('Input audio is not sampled at 48 kHz.');
        return
    end

    disp('Processing audio with generated room ...');
    processed_audio = conv(loaded_data.reflectogram, input_audio);
    processed_audio = processed_audio / max(abs(processed_audio));

    disp('Saving processed audio to file ...');
    audiowrite('processed_audio.wav', processed_audio, sample_rate);

    disp('[Original audio]');
    p = audioplayer(input_audio, sample_rate);
    play(p)
    disp('  Press a key to continue ...');
    pause
    stop(p)

    disp('[Original audio]');
    p = audioplayer(processed_audio, sample_rate);
    play(p)
    disp('  Press a key to continue ...');
    pause
    stop(p)

end

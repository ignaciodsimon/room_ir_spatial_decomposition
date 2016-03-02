function display_results(results_filename)
    % Generates a plot showing the results of the simulation. It includes
    % the boundaries of the room, the source and microphone positions and
    % the found image sources. The reflectogram is also plotted with a time
    % axis.
    %
    % Joe.

    % Load data from simulation result file
    loaded_data = load(results_filename);
    reflectogram = loaded_data.reflectogram;
    sample_rate = loaded_data.sample_rate;
    mic_position = loaded_data.mic_position;
    source_position = loaded_data.source_position;
    room_boundaries = loaded_data.room_boundaries;
    image_sources_list = loaded_data.image_sources_list;

    % Show plot with image sources and reflectogram
    figure();
    time_axis = [0 : 1/sample_rate : (length(reflectogram)-1)/sample_rate];
    figure()
    subplot(2,1,1)
    plot(mic_position(1), mic_position(2), 'o')
    hold on
    plot(source_position(1), source_position(2), 'x')
    for curr_boundary = room_boundaries'
        plot([curr_boundary(1) curr_boundary(3)], [curr_boundary(2) curr_boundary(4)], 'Color', [.2 .8 .3])
    end
    for i = 1 : size(image_sources_list, 1)
        plot(image_sources_list(i,1), image_sources_list(i,2), 'o', 'MarkerFaceColor', [.3 .6 1], 'Color', [0 0 0], 'MarkerSize', image_sources_list(i,3) * 100 + 2);
    end
    xlabel('X-coordinate [m]')
    ylabel('Y-coordinate [m]')
    title('Image source location / strength');
    grid
    subplot(2,1,2)
    plot(time_axis, reflectogram);
    title('Obtained reflectogram (using fractional delay)')
    ylabel('Amplitude [.]')
    xlabel('Time [s]')
    grid

end

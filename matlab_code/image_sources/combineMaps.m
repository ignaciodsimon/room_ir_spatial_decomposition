function combineMaps()

    mapsFilenames = {'simulation_results_no_diffusion_rectangular_room_1.mat', ...
                     'simulation_results_no_diffusion_rectangular_room_2.mat', ...
                     'simulation_results_no_diffusion_rectangular_room_4.mat', ...
                     'simulation_results_no_diffusion_rectangular_room_5.mat', ...
                     };


    combinedMap = [];
    for index = 1 : length(mapsFilenames)
        newMapData = load(char(mapsFilenames(index)));

        if isempty(combinedMap)
            combinedMap = newMapData;
        else
            combinedMap.image_sources_list = [combinedMap.image_sources_list
                                              newMapData.image_sources_list];
            combinedMap.source_position =    [combinedMap.source_position
                                              newMapData.source_position];
        end
    end

    convolved_ir = combinedMap.convolved_ir;
    image_sources_list = combinedMap.image_sources_list;
    mic_position = combinedMap.mic_position;
    reflectogram = combinedMap.reflectogram;
    room_boundaries = combinedMap.room_boundaries;
    sample_rate = combinedMap.sample_rate;
    source_position = combinedMap.source_position;
    time_axis = combinedMap.time_axis;
    save('simulation_results_no_diffusion_rectangular_room_combinedMap_several.mat', 'convolved_ir', 'image_sources_list', 'mic_position', ...
         'reflectogram', 'room_boundaries', 'sample_rate', 'source_position', ...
         'time_axis');
end

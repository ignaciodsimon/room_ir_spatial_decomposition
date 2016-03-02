% Run simulation using the pre-defined rectangular room with the materials
% definitions taken from the book "Environmental and Architectural Acoustics".
%
% Joe.

% Update simulation parameters
synthetic_room_geometry(0);

% Perform simulation and obtain results
[image_sources_list, reflectogram] = mirror_images('synthetic_room_geometry.mat');

% Save results to a file
loaded_data = load('synthetic_room_geometry.mat');
time_axis = [0 : 1/loaded_data.sample_rate : (length(reflectogram)-1)/loaded_data.sample_rate];
sample_rate = loaded_data.sample_rate;
mic_position = loaded_data.mic_position;
source_position = loaded_data.source_position;
room_boundaries = loaded_data.room_boundaries;
save('simulation_results.mat', 'image_sources_list', 'reflectogram', 'time_axis', ...
     'sample_rate', 'mic_position', 'source_position', 'room_boundaries', 'image_sources_list');

% Displays results
display_results('simulation_results.mat')

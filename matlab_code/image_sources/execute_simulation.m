% Run simulation using the pre-defined rectangular room with the materials
% definitions taken from the book "Environmental and Architectural Acoustics".
%
% Joe.

% Update simulation parameters
disp('>> Regenerating simulation parameters ...')
% addpath('rooms/synthetic_room/');
% synthetic_room_geometry(0);
% addpath('rooms/irregular_room_2/');
% irregular_room_2_geometry(0);
addpath('rooms/rectangular_room/');
rectangular_room_geometry(0);

% Perform simulation and obtain results
disp(sprintf('>> Performing simulation ...\n'))
% [image_sources_list, reflectogram] = mirror_images('synthetic_room_geometry.mat');
% [image_sources_list, reflectogram] = mirror_images('irregular_room_geometry_2.mat');
[image_sources_list, reflectogram] = mirror_images('rectangular_room_geometry.mat');

% Obtain response with loudspeaker IR
disp(sprintf('\n>> Convolving reflectogram to generate IR ...'));
genelec_data = load('genelec_response/genelec_IR_processed.mat');
convolved_ir = conv(reflectogram, genelec_data.ir);
convolved_ir = convolved_ir(1:length(reflectogram));

% Save results to a file
disp('>> Saving all results to a file');
% loaded_data = load('synthetic_room_geometry.mat');
% loaded_data = load('irregular_room_geometry_2.mat');
loaded_data = load('rectangular_room_geometry.mat');
time_axis = [0 : 1/loaded_data.sample_rate : (length(reflectogram)-1)/loaded_data.sample_rate];
sample_rate = loaded_data.sample_rate;
mic_position = loaded_data.mic_position;
source_position = loaded_data.source_position;
room_boundaries = loaded_data.room_boundaries;
% save('simulation_results_no_diffusion.mat', 'image_sources_list', 'reflectogram', ...
%      'time_axis', 'sample_rate', 'mic_position', 'source_position', ...
%      'room_boundaries', 'image_sources_list', 'convolved_ir');
save('simulation_results_no_diffusion_rectangular_room.mat', 'image_sources_list', 'reflectogram', ...
     'time_axis', 'sample_rate', 'mic_position', 'source_position', ...
     'room_boundaries', 'image_sources_list', 'convolved_ir');

% Displays results
disp('>> Plotting results.')
% display_results('simulation_results_no_diffusion.mat')
display_results('simulation_results_no_diffusion_rectangular_room.mat')

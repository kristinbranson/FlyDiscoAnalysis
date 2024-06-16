experiment_folder_path = '/groups/branson/bransonlab/taylora/flydisco/pflySingleMale_pfly_SlowRamp_71G01_UASChrimsonVenusX0070_20240320T102203' ;
dot_movie_path = fullfile(experiment_folder_path, 'dotmovie.ufmf') ;
%playfmf([], dotmovie_path) ;
movie_path = fullfile(experiment_folder_path, 'movie.ufmf') ;
%playfmf([], movie_path) ;

% ellipseTrajectory.csv contains a look-up table from indices to <x,y>
% actualTraj.csv contains a LUT index for each frame
ellipse_trajectory_csv_path = fullfile(experiment_folder_path, 'ellipseTrajectory.csv') ;
%ellipse_trajectory_struct = importdata(ellipse_trajectory_csv_path, ',', 1) ;
ellipse_trajectory = readmatrix(ellipse_trajectory_csv_path) ;
dimension_count = 2 ;
raw_xy_from_lut_index = [1 -1]' .* ellipse_trajectory(:,1:dimension_count)' ;  % 2 x 471, y sign convention is opposite of what we want
pfly_radius_in_mm = mean(vecnorm(raw_xy_from_lut_index)) ;  % mm
theta_hat_from_lut_index = raw_xy_from_lut_index / pfly_radius_in_mm ; % want a unit vector at each index

figure('color', 'w') ;
plot(theta_hat_from_lut_index(1,:), theta_hat_from_lut_index(2,:), 'k.') ;
axis equal
title('LUT (ellipseTrajectory.csv) x-y values') ;
xlabel('x') ;
ylabel('y') ;
% looks like a circle, early points are on the x-axis and points make a CCW
% circle

theta_from_lut_index = unwrap(atan2(theta_hat_from_lut_index(2,:), theta_hat_from_lut_index(1,:))) ;
% theta starts and zero and increases steadily until it is equal to
% 2*pi-epsilon

heading_angle_from_lut_index = -ellipse_trajectory(:,3)' ;

% How many frames in the video
movie_header = ufmf_read_header(movie_path)
frame_count = movie_header.nframes   % -> 61237, same as line count of stamp_log_cam0.txt
dot_movie_header = ufmf_read_header(dot_movie_path)
dot_movie_frame_count = dot_movie_header.nframes  % 61232, same as line count of stamp_log_cam1.txt
min_frame_count = min(frame_count, dot_movie_frame_count) ;

%% 
actual_trajectory_csv_path = fullfile(experiment_folder_path, 'actualTraj.csv') ;
actual_trajectory = load_actual_traj_csv(actual_trajectory_csv_path) ;
% 1st col a timestamp, 2nd col the lut index, 3rd col 'on' or 'off' depending
% on whether the fake fly is on or off.
timestamp_from_tick_index = actual_trajectory(:,1)' - actual_trajectory(1,1) ;  % s
lut_index_from_tick_index = actual_trajectory(:,2)' ;
is_pfly_on_from_tick_index = actual_trajectory(:,3)' ;
tick_count = numel(timestamp_from_tick_index)  % -> 108914.  NB: Lots more than the frame count
  % we use "tick" here instead of "frame" just because there are more of them
  % than the frame_count

%%
figure('color', 'w') ;
plot(timestamp_from_tick_index,'k') ;
xlabel('Tick index') ;
ylabel('Timestamp (s)') ;
title('pfly timestamp vs tick index (from actualTraj.csv)')

dtimestamp_from_tick_index = diff(timestamp_from_tick_index) ;  % s
figure('color', 'w') ;
plot(timestamp_from_tick_index(1:end-1),dtimestamp_from_tick_index,'k') ;
xlabel('Timestamp (s)') ;
ylabel('Timestamp delta (s)') ;
title('pfly Timestamp delta vs timestamp')

mean_tick_interval = mean(dtimestamp_from_tick_index)   % s
mean_tick_rate = 1/mean_tick_interval  % hz
%nominal_timestamp_from_tick_index = mean_tick_interval * (0:(tick_count-1))' ;  % s

%%
theta_hat_from_tick_index = theta_hat_from_lut_index(:,lut_index_from_tick_index) ;

%tick_index_from_frame_index = (1:frame_count)' ;

%theta_hat_from_frame_index = theta_hat_from_tick_index(tick_index_from_frame_index,:) ;
figure('color', 'w') ;
plot(theta_hat_from_tick_index(1,:), theta_hat_from_tick_index(2,:), 'k.') ;
title('x-y values, tick index explicit (from actualTraj.csv)') ;
xlabel('x') ;
ylabel('y') ;
axis equal

% In dotmovie.ufmf, the first frame with a visible pfly is 5, the last is
% 61205 (using 1-based indexing).

% What do the timestamps in the timestamp files look like?
movie_timestamp_file_path = fullfile(experiment_folder_path, 'stamp_log_cam0.txt') ;
timestamp_from_movie_frame_index = importdata(movie_timestamp_file_path)' ;
dot_movie_timestamp_file_path = fullfile(experiment_folder_path, 'stamp_log_cam1.txt') ;
timestamp_from_dot_movie_frame_index = importdata(dot_movie_timestamp_file_path)' ;
figure('color', 'w') ;
plot(timestamp_from_movie_frame_index(1:min_frame_count), timestamp_from_dot_movie_frame_index(1:min_frame_count), 'k') ;
xlabel('Movie timestamp (s)') ;
ylabel('Dot movie timestamp (s)') ;
title('Camera timestamp comparison') ;

%%
camera_settings_frame_rate = 59.99542999267578125   % Hz, as specified in cameraSettings0.json, cameraSettings1.json, approx == 7863721/2^17
camera_settings_frame_interval = 1/camera_settings_frame_rate 
camera_settings_timestamp_from_movie_frame_index = camera_settings_frame_interval * (0:(frame_count-1)) ;
camera_settings_timestamp_from_dot_movie_frame_index = camera_settings_frame_interval * (0:(dot_movie_frame_count-1)) ;

movie_timestamp_diff = timestamp_from_movie_frame_index - camera_settings_timestamp_from_movie_frame_index ;
dot_movie_timestamp_diff = timestamp_from_dot_movie_frame_index - camera_settings_timestamp_from_dot_movie_frame_index ;

figure('color', 'w') ;
line_handles = ...
  plot(camera_settings_timestamp_from_movie_frame_index, 1000*movie_timestamp_diff, 'r', ...
       camera_settings_timestamp_from_dot_movie_frame_index, 1000*dot_movie_timestamp_diff, 'b') ;
xlabel('Nominal timestamp (s)') ;
ylabel('Timestamp diff from nominal (ms)') ;
legend(line_handles, {'Movie', 'Dot movie'}, 'location', 'northwest') ;
title('Timestamps vs nominal time') ;

% Let's check the pfly position in actualTraj.csv vs the position in
% dotmovie.ufmf

% First interpolate to get pfly x,y at each movie frame timestamp
theta_hat_from_dot_movie_frame_index = interp1(timestamp_from_tick_index', theta_hat_from_tick_index', timestamp_from_dot_movie_frame_index')' ;
  % dot_movie_frame_count x 2 
fractional_is_pfly_on_from_dot_movie_frame_index = interp1(timestamp_from_tick_index', is_pfly_on_from_tick_index', timestamp_from_dot_movie_frame_index')' ;
is_pfly_on_from_dot_movie_frame_index = (fractional_is_pfly_on_from_dot_movie_frame_index>0.5) ;

% plot those  
figure('color', 'w') ;
line_handles = ...
  plot(timestamp_from_dot_movie_frame_index, theta_hat_from_dot_movie_frame_index(1,:), 'r', ...
       timestamp_from_dot_movie_frame_index, theta_hat_from_dot_movie_frame_index(2,:), 'b') ;
xlabel('Movie timestamp (s)') ;
ylabel('Coordinate') ;
legend(line_handles, {'x', 'y'}, 'location', 'northeast') ;
title('Interpolated position vs movie timestamp') ;
xlim([800 850]) ;  % zoom in on a random part---too much data otherwise
% Looks reasonable

%%
% Load the flytracker info about the arena locations and size
ft_calibration_path = fullfile(experiment_folder_path, 'flytracker-calibration.mat') ;
calib = load_anonymous(ft_calibration_path) ;
raw_center_from_arena_index = flipud((calib.centroids)')   % 2x9, pels
arena_count = size(raw_center_from_arena_index, 2) 
arena_radius = calib.r  % scalar, pels
% One issue: turns out these are correct for the main movie, but the dot movie
% is shifted somewhat.
center_from_arena_index = raw_center_from_arena_index + [0 16]' ;  % shift for dot-movie
pixels_per_mm = calib.PPM ;  % pixels per mm
% pfly_radius = 0.85*125  % pels
pfly_radius = pixels_per_mm * pfly_radius_in_mm  % pels

% Generate the pfly position for each arena, for each frame
big_center_from_arena_index = reshape(center_from_arena_index, [dimension_count 1 arena_count]) ;
position_from_dot_movie_frame_index_from_arena_index = big_center_from_arena_index + pfly_radius * theta_hat_from_dot_movie_frame_index ;

% Read a frame in the middle of the dot movie
dot_movie_frame_index = 60000 ;
dot_movie_frame = ufmf_read_frame(dot_movie_header, dot_movie_frame_index) ;
[f, a, ih] = imglance(dot_movie_frame) ;  % show the frame
set_figure_size_in_pixels(f, [1100 1100]) ;
a.Position = [ 0 0 1 1 ] ;

% For the current frame, get the pfly positions
position_from_arena_index = reshape(position_from_dot_movie_frame_index_from_arena_index(:,dot_movie_frame_index,:), ...
                                    [dimension_count arena_count]) ;

% Plot the arena centers on the image
center_dots_handle = ...
  line('parent', a, ...
       'xdata', center_from_arena_index(1,:), ...
       'ydata', center_from_arena_index(2,:), ...
       'linestyle', 'none', ...
       'marker', 'o', ...
       'markersize', 9, ...
       'markerfacecolor', 'b', ...
       'markeredgecolor', 'none') ;

% Plot the arena boundaries
theta_from_index = linspace(0,2*pi,361) ;
for arena_index = 1 : arena_count ,
  x = center_from_arena_index(1,arena_index) ;
  y = center_from_arena_index(2,arena_index) ;
  line(a, x + arena_radius*cos(theta_from_index), y + arena_radius*sin(theta_from_index)', 'color', [0 0.7 0], 'linewidth', 2) ;
end

% Plot those on the image
pfly_dots_handle = ...
  line('parent', a, ...
       'xdata', position_from_arena_index(1,:), ...
       'ydata', position_from_arena_index(2,:), ...
       'linestyle', 'none', ...
       'marker', 'o', ...
       'markersize', 9, ...
       'markerfacecolor', 'r', ...
       'markeredgecolor', 'none') ;

%%
% Plot that every n frames
stride_count = 1000 ;
dot_movie_frame_index_from_output_movie_frame_index = (1:stride_count:dot_movie_frame_count) ;
output_movie_frame_count = numel(dot_movie_frame_index_from_output_movie_frame_index) ;

% Get the output movie started
output_movie_avi_file_name = 'pfly-test-result.avi' ;
profile = 'Motion JPEG AVI';
vw = VideoWriter(output_movie_avi_file_name, profile) ;
vw.FrameRate = 2 ;
vw.Quality = 100 ;
vw.open() ;

for output_movie_frame_index = 1 : output_movie_frame_count ,
  dot_movie_frame_index = dot_movie_frame_index_from_output_movie_frame_index(output_movie_frame_index) ;
  dot_movie_frame = ufmf_read_frame(dot_movie_header, dot_movie_frame_index) ;
  % For the current frame, get the pfly positions
  position_from_arena_index = reshape(position_from_dot_movie_frame_index_from_arena_index(:,dot_movie_frame_index,:), ...
                                      [dimension_count arena_count]) ;
  
  % plot
  if output_movie_frame_index == 1 ,
    [f, a, ih] = imglance(dot_movie_frame) ;  % show the frame
    set_figure_size_in_pixels(f, [1100 1100]) ;
    a.Position = [ 0 0 1 1 ] ;
    % Plot the arena centers on the image
    center_dots_handle = ...
      line('parent', a, ...
           'xdata', center_from_arena_index(1,:), ...
           'ydata', center_from_arena_index(2,:), ...
           'linestyle', 'none', ...
           'marker', 'o', ...
           'markersize', 9, ...
           'markerfacecolor', 'b', ...
           'markeredgecolor', 'none') ;
    % Plot the arena boundaries
    theta_from_index = linspace(0,2*pi,361) ;
    for arena_index = 1 : arena_count ,
      x = center_from_arena_index(1,arena_index) ;
      y = center_from_arena_index(2,arena_index) ;
      line(a, x + arena_radius*cos(theta_from_index), y + arena_radius*sin(theta_from_index)', 'color', [0 0.7 0], 'linewidth', 2) ;
    end
    % Plot those on the image
    pfly_dots_handle = ...
      line('parent', a, ...
           'xdata', position_from_arena_index(1,:), ...
           'ydata', position_from_arena_index(2,:), ...
           'linestyle', 'none', ...
           'marker', 'o', ...
           'markersize', 9, ...
           'markerfacecolor', 'r', ...
           'markeredgecolor', 'none') ;
  else
    % if not 1st output frame
    ih.CData = dot_movie_frame ;
    pfly_dots_handle.XData = position_from_arena_index(1,:) ;
    pfly_dots_handle.YData = position_from_arena_index(2,:) ;
  end

  % Write the frame to the video
  output_movie_frame = getframe(a);
  if output_movie_frame_index == 1 ,
    height = size(output_movie_frame.cdata,1);
    width = size(output_movie_frame.cdata,2);
  else
    output_movie_frame.cdata = tweak_rgb_image_size(output_movie_frame.cdata, height, width) ;  % sometimes subsequent frames aren't quite the right size.
  end
  vw.writeVideo(output_movie_frame) ;
end
vw.close() ;


%%

%
% Compress the .avi file to h.264
%

% Tweak the frame size to be multiples of 4
newheight = 4*ceil(height/4);
newwidth = 4*ceil(width/4);

% Generate the two command-line commands we need
mp4_file_path = replace_extension(output_movie_avi_file_name, '.mp4') ;
ffmpeg_command = 'env -u LD_LIBRARY_PATH /usr/bin/ffmpeg' ;  
cmd = ...
  sprintf('%s -i %s -y -c:v h264 -pix_fmt yuv420p -s %dx%d -b:v 1600k -f mp4 %s',...
          ffmpeg_command, output_movie_avi_file_name, newwidth, newheight, mp4_file_path) ;

% Run the second command, deal with any error
system_with_error_handling(cmd) ;

% If both commands succeeded, clean up the intermediate files that are no longer needed
delete(output_movie_avi_file_name);

%%

% Read in the trx file produced by FlyTracker
flytracker_trx_file_path = fullfile(experiment_folder_path, 'movie_JAABA/trx.mat') ;
s = load(flytracker_trx_file_path);
trx_timestamp_from_movie_frame_index = s.timestamps - s.timestamps(1) ;
ft_trx = s.trx ;  % ft_ is for FlyTracker

% How do these timestamps compare to the ones from the camera timestamp file?
trx_timestamp_diff_from_movie_frame_index = trx_timestamp_from_movie_frame_index - timestamp_from_movie_frame_index ;
figure('color', 'w') ; 
plot(timestamp_from_movie_frame_index, 1000*trx_timestamp_diff_from_movie_frame_index, 'k') ;
xlabel('Timestamp from stamp_log file (s)', 'interpreter', 'none') ;
ylabel('Diff with trx timestamps (ms)') ;
title('Timestamps in trx file vs stamp_log file', 'interpreter', 'none') ;
% Goes from zero up to -80 ms at end.  That's several frame interval's worth.

% What is the frame interval from the .trx file?
mean_trx_timestamp_interval = mean(diff(trx_timestamp_from_movie_frame_index))
sd_trx_timestamp_interval = std(diff(trx_timestamp_from_movie_frame_index))
% Mean is 0.0166666666666667, SD is 3.0481100383374e-14
% So that's consistent with a frame rate of 60 Hz, which is not *exactly*
% right.  Does that add up to 80 ms difference by the end?

timestamp_diff_at_end_in_theory_in_ms = 1000*(frame_count*mean_trx_timestamp_interval - frame_count*camera_settings_frame_interval)  % ms
% About -77.7 ms.  So that's most of the difference.  Sigh.

% Think we should use the timestamps from the stamp-log file.  They're likely
% to be more accurate.

% Read in the fake-fly config, to get fake fly ellipse dimensions
fake_fly_params_file_path = fullfile(experiment_folder_path, 'ellipseTrajectory_config.csv') ;
fake_fly_params =  readtable(fake_fly_params_file_path); 

%%
% Create the trx struct array for the pflies
maximum_real_fly_id = max([ft_trx.id]) ;
arena_center_shift = [ 0 16 ]' ;  % pels
do_debug = false ; 
fake_trx = ...
  fake_fly_trx_from_various_inputs(dot_movie_path, timestamp_from_dot_movie_frame_index, ellipse_trajectory, actual_trajectory, ...
                                   calib, fake_fly_params, maximum_real_fly_id, arena_center_shift, do_debug)

trx = [ ft_trx fake_trx ] ;
timestamps = (1/calib.FPS) * ( 0 : frame_count )' ;
save('trx-with-pfly.mat', 'trx', 'timestamps') ;

make_ctrax_result_movie('moviename', 'pflySingleMale_pfly_SlowRamp_71G01_UASChrimsonVenusX0070_20240320T102203/movie.ufmf', ...
                        'trxname', 'trx-with-pfly.mat', ...
                        'aviname', 'trx-with-pfly.avi', ...
                        'doshowsex', false, ...
                        'doesYAxisPointUp', false) ;

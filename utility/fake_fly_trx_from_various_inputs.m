function fake_trx = fake_fly_trx_from_various_inputs(movie_path, timestamp_from_frame_index, ellipse_trajectory, actual_trajectory, ...
                                                     flytracker_calibration, fake_fly_params, maximum_real_fly_id, arena_center_shift, do_debug)

% Computes the trx struct array for the fake flies.  This is a pure function
% unless do_debug is true, in which case the movie file is read from and a
% debugging movie is output.
%
% movie_path: the path to the movie file, usually named movie.ufmf
% timestamp_from_frame_index: The time of each movie frame, in seconds.
%   Normally obtained from reading the stamp_log_cam<n>.txt file for the movie.
%   This can be somewhat different than the timestamps you would get from
%   reading the fps and frame count from flytracker_calibration and
%   constructing a timeline.  Mainly because sometimes people tell FlyTracker
%   the fps is 60, when it's actually e.g. 59.99542999267578125.
% ellipse_trajectory: 1st col the x coord, 2nd col the y coord, 3rd col the
%   heading angle in radians.  Number of rows is the number of entries in the
%   look-up table (LUT).  x,y are in typical cartesian coordinates (y axis
%   points up), heading angle uses usual convention (angle relative to x-axis,
%   increasing counterclockwise)
% actual_trajectory: 1st col a timestamp (in seconds); 2nd col the LUT index;
%   3rd col 0 or 1, depending on whether the pfly ("projected fly") is visible or not
% flytracker_calibration: The "calib" structure from a FlyTracker calibration
%   file.
% fake_fly_params: table or struct containing 'maj_axis' and 'min_axis'
%   columns/fields that give the fake fly major/minor axis length.  This is
%   the full length, the semi- axis length or quarter- axis length.  In mm.
% maximum_real_fly_id: What the tin says.  Used to make the fake fly ids
%   distinct.
% arena_center_shift: 2x1, in pixels.  How much the arena centers are shifted
%   relative to values in flytracker_calibration.  Useful for movies that show
%   the pfly dots.  Defaults to [0 0]' if empty.
% do_debug: If true, debugging plots are made and a debug movie is written
%   into the current directory.  Defaults to false if missing or empty.

% Deal with args
if ~exist('arena_center_shift', 'var') || isempty(arena_center_shift) ,
  arena_center_shift = [0 0]' ;
end
if ~exist('do_debug', 'var') || isempty(do_debug) ,
  do_debug = false ;
end
do_debug = false ;  % debug code is not working right now

%experiment_folder_path = '/groups/branson/bransonlab/taylora/flydisco/pflySingleMale_pfly_SlowRamp_71G01_UASChrimsonVenusX0070_20240320T102203' ;
%dotmovie_path = fullfile(experiment_folder_path, 'dotmovie.ufmf') ;
%playfmf([], dotmovie_path) ;
%movie_path = fullfile(experiment_folder_path, 'movie.ufmf') ;
%playfmf([], movie_path) ;

% ellipseTrajectory.csv contains a look-up table from indices to <x,y>
% actualTraj.csv contains a LUT index for each frame
%ellipse_trajectory_csv_path = fullfile(experiment_folder_path, 'ellipseTrajectory.csv') ;
%ellipse_trajectory_struct = importdata(ellipse_trajectory_csv_path, ',', 1) ;
%ellipse_trajectory = readmatrix(ellipse_trajectory_csv_path) ;
dimension_count = 2 ;
raw_pfly_xy_from_lut_index = [1 -1]' .* ellipse_trajectory(:,1:dimension_count)' ;  % 2 x 471, y sign convention is opposite of what we want
pfly_radius_in_mm = mean(vecnorm(raw_pfly_xy_from_lut_index)) ;  % mm
pfly_theta_hat_from_lut_index = raw_pfly_xy_from_lut_index / pfly_radius_in_mm ; % want a unit vector at each index

if do_debug ,
  figure('color', 'w') ;
  plot(pfly_theta_hat_from_lut_index(1,:), pfly_theta_hat_from_lut_index(2,:), 'k.') ;
  axis equal
  title('LUT (ellipseTrajectory.csv) x-y values') ;
  xlabel('x') ;
  ylabel('y') ;
  % looks like a circle, early points are on the x-axis and points make a CCW
  % circle
end

pfly_theta_from_lut_index = unwrap(atan2(pfly_theta_hat_from_lut_index(2,:), pfly_theta_hat_from_lut_index(1,:))) ;
% theta starts and zero and increases steadily until it is equal to
% 2*pi-epsilon
pfly_heading_angle_from_lut_index = -ellipse_trajectory(:,3)' ;
assert(all(abs(pfly_heading_angle_from_lut_index-pfly_theta_from_lut_index+pi/2)<1e-9))
% the pfly starts at 3 o'clock, goes around CCW, so the heading should be
% theta-pi/2 in the "CRTesian" (y-axis points down) coordinate system.
pfly_heading_hat_from_lut_index = [ cos(pfly_heading_angle_from_lut_index) ; sin(pfly_heading_angle_from_lut_index) ] ;

% How many frames in the video
frame_count = numel(timestamp_from_frame_index) ;  %#ok<NASGU> 

%% 
%actual_trajectory_csv_path = fullfile(experiment_folder_path, 'actualTraj.csv') ;
%actual_trajectory = load_actual_traj_csv(actual_trajectory_csv_path) ;
% 1st col a timestamp, 2nd col the lut index, 3rd col 'on' or 'off' depending
% on whether the fake fly is on or off.
timestamp_from_tick_index = actual_trajectory(:,1)' - actual_trajectory(1,1) ;  % s
lut_index_from_tick_index = actual_trajectory(:,2)' ;
is_pfly_on_from_tick_index = actual_trajectory(:,3)' ;
% tick_count = numel(timestamp_from_tick_index) ;  % Usually many more ticks than frames
%   % we use "tick" here instead of "frame" just because there are more of them
%   % than the frame_count

%%
if do_debug ,
  dtimestamp_from_tick_index = diff(timestamp_from_tick_index) ;  % s
  figure('color', 'w') ;
  plot(timestamp_from_tick_index,'k') ;
  xlabel('Tick index') ;
  ylabel('Timestamp (s)') ;
  title('pfly timestamp vs tick index (from actualTraj.csv)')

  figure('color', 'w') ;
  plot(timestamp_from_tick_index(1:end-1),dtimestamp_from_tick_index,'k') ;
  xlabel('Timestamp (s)') ;
  ylabel('Timestamp delta (s)') ;
  title('pfly Timestamp delta vs timestamp')
end

%mean_tick_interval = mean(dtimestamp_from_tick_index) ;  % s
%mean_tick_rate = 1/mean_tick_interval  % hz

%%
pfly_theta_hat_from_tick_index = pfly_theta_hat_from_lut_index(:,lut_index_from_tick_index) ;
pfly_heading_hat_from_tick_index = pfly_heading_hat_from_lut_index(:,lut_index_from_tick_index) ;

if do_debug ,
  figure('color', 'w') ;
  plot(pfly_theta_hat_from_tick_index(1,:), pfly_theta_hat_from_tick_index(2,:), 'k.') ;
  title('Pfly x-y values, tick index explicit') ;
  xlabel('x') ;
  ylabel('y') ;
  axis equal
end

%%

% First interpolate to get pfly x,y at each movie frame timestamp
pfly_theta_hat_from_frame_index = interp1(timestamp_from_tick_index', pfly_theta_hat_from_tick_index', timestamp_from_frame_index')' ;
  % frame_count x 2 
fractional_is_pfly_on_from_frame_index = interp1(timestamp_from_tick_index', is_pfly_on_from_tick_index', timestamp_from_frame_index')' ;
is_pfly_on_from_frame_index = (fractional_is_pfly_on_from_frame_index>0.5) ;

% Do the same for the heading
pfly_heading_hat_from_frame_index = interp1(timestamp_from_tick_index', pfly_heading_hat_from_tick_index', timestamp_from_frame_index')' ;
  % frame_count x 2 

% plot those  
if do_debug ,
  figure('color', 'w') ;
  line_handles = ...
    plot(timestamp_from_frame_index, pfly_theta_hat_from_frame_index(1,:), 'r', ...
         timestamp_from_frame_index, pfly_theta_hat_from_frame_index(2,:), 'b') ;
  xlabel('Movie timestamp (s)') ;
  ylabel('Coordinate') ;
  legend(line_handles, {'x', 'y'}, 'location', 'northeast') ;
  title('Interpolated position vs movie timestamp') ;
  xlim([800 850]) ;  % zoom in on a random part---too much data otherwise
end

%%
% Load the flytracker info about the arena locations and size
%ft_calibration_path = fullfile(experiment_folder_path, 'flytracker-calibration.mat') ;
%flytracker_calibration = load_anonymous(ft_calibration_path) ;
raw_center_from_arena_index = flipud((flytracker_calibration.centroids)') ;   % 2x9, pels
arena_count = size(raw_center_from_arena_index, 2) ; 
%arena_radius = flytracker_calibration.r  % scalar, pels
% One issue: turns out these are correct for the main movie, but the dot movie
% is shifted somewhat.
%center_from_arena_index = raw_center_from_arena_index + [0 16]' ;  % shift for dot-movie
center_from_arena_index = raw_center_from_arena_index + arena_center_shift ;
pixels_per_mm = flytracker_calibration.PPM ;  % pixels per mm
% pfly_radius = 0.85*125  % pels
pfly_radius = pixels_per_mm * pfly_radius_in_mm ;  % pels

% Generate the pfly position for each arena, for each frame
big_center_from_arena_index = reshape(center_from_arena_index, [dimension_count 1 arena_count]) ;
pfly_position_from_frame_index_from_arena_index = big_center_from_arena_index + pfly_radius * pfly_theta_hat_from_frame_index ;

if do_debug, 
  % Read a frame in the middle of the dot movie
  frame_index = round((frame_count+1)/2) ;
  dot_movie_frame = ufmf_read_frame(dot_movie_header, frame_index) ;
  [f, a, ih] = imglance(dot_movie_frame) ;  % show the frame
  set_figure_size_in_pixels(f, [1100 1100]) ;
  a.Position = [ 0 0 1 1 ] ;
  
  % For the current frame, get the pfly positions
  position_from_arena_index = reshape(pfly_position_from_frame_index_from_arena_index(:,frame_index,:), ...
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
         'markeredgecolor', 'none') ;  %#ok<NASGU> 
  
  % Plot the arena boundaries
  theta_from_index = linspace(0,2*pi,361) ;
  for arena_index = 1 : arena_count ,  %#ok<FXUP> 
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
  frame_index_from_debug_movie_frame_index = (1:stride_count:frame_count) ;
  debug_movie_frame_count = numel(frame_index_from_debug_movie_frame_index) ;
  
  % Get the output movie started
  output_movie_avi_file_name = sprintf('%s.avi', mfilename()) ;
  profile = 'Motion JPEG AVI';
  vw = VideoWriter(output_movie_avi_file_name, profile) ;
  vw.FrameRate = 2 ;
  vw.Quality = 100 ;
  vw.open() ;
  
  for debug_movie_frame_index = 1 : debug_movie_frame_count ,
    frame_index = frame_index_from_debug_movie_frame_index(debug_movie_frame_index) ;
    dot_movie_frame = ufmf_read_frame(dot_movie_header, frame_index) ;
    % For the current frame, get the pfly positions
    position_from_arena_index = reshape(pfly_position_from_frame_index_from_arena_index(:,frame_index,:), ...
                                        [dimension_count arena_count]) ;
    
    % plot
    if debug_movie_frame_index == 1 ,
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
             'markeredgecolor', 'none') ;  %#ok<NASGU> 
      % Plot the arena boundaries
      theta_from_index = linspace(0,2*pi,361) ;
      for arena_index = 1 : arena_count ,  %#ok<FXUP> 
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
    if debug_movie_frame_index == 1 ,
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
end

%%

% % Read in the trx file produced by FlyTracker
% flytracker_trx_file_path = fullfile(experiment_folder_path, 'movie_JAABA/trx.mat') ;
% s = load(flytracker_trx_file_path);
% trx_timestamp_from_frame_index = s.timestamps - s.timestamps(1) ;
% ft_trx = s.trx ;  % ft_ is for FlyTracker
% 
% % How do these timestamps compare to the ones from the camera timestamp file?
% trx_timestamp_diff_from_frame_index = trx_timestamp_from_frame_index - timestamp_from_frame_index ;
% figure('color', 'w') ; 
% plot(timestamp_from_frame_index, 1000*trx_timestamp_diff_from_frame_index, 'k') ;
% xlabel('Timestamp from stamp_log file (s)', 'interpreter', 'none') ;
% ylabel('Diff with trx timestamps (ms)') ;
% title('Timestamps in trx file vs stamp_log file', 'interpreter', 'none') ;
% % Goes from zero up to -80 ms at end.  That's several frame interval's worth.
% 
% % What is the frame interval from the .trx file?
% mean_trx_timestamp_interval = mean(diff(trx_timestamp_from_frame_index))
% sd_trx_timestamp_interval = std(diff(trx_timestamp_from_frame_index))
% % Mean is 0.0166666666666667, SD is 3.0481100383374e-14
% % So that's consistent with a frame rate of 60 Hz, which is not *exactly*
% % right.  Does that add up to 80 ms difference by the end?
% 
% timestamp_diff_at_end_in_theory_in_ms = 1000*(frame_count*mean_trx_timestamp_interval - frame_count*camera_settings_frame_interval)  % ms
% % About -77.7 ms.  So that's most of the difference.  Sigh.

% Think we should use the timestamps from the stamp-log file.  They're likely
% to be more accurate.

% Read in the fake-fly config
% fake_fly_params_file_path = fullfile(experiment_folder_path, 'ellipseTrajectory_config.csv') ;
% fake_fly_params =  readtable(fake_fly_params_file_path); 
fake_fly_a_in_mm = fake_fly_params.maj_axis/4 ;
  % mm, length of longest line segment that spans ellipse (equal to diameter for a circle)
fake_fly_b_in_mm = fake_fly_params.min_axis/4 ; 
  % mm, length of longest line segment that spans ellipse and is perpendicular to major axis (equal to diameter for a circle)
fake_fly_a = pixels_per_mm * fake_fly_a_in_mm ;
fake_fly_b = pixels_per_mm * fake_fly_b_in_mm ;

% Figure out all the other trx fields we'll need
[~, movie_leaf_name] = fileparts2(movie_path) ;
moviename = movie_path ;
experiment_folder_path = fileparts(movie_path) ;
moviefile = fullfile(experiment_folder_path, 'movie_JAABA', movie_leaf_name) ;
firstframe = find(is_pfly_on_from_frame_index, 1) ;
off = 1 - firstframe ;  % ??
endframe = find(is_pfly_on_from_frame_index, 1, 'last') ;
nframes = endframe - firstframe + 1 ;
fps = flytracker_calibration.FPS ;  % This is 60 Hz, slightly different than camera settings frame rate
pxpermm = pixels_per_mm ;
% id will need to be per-pfly
sex = 'b' ;  % Is this right?
frame_count = numel(timestamp_from_frame_index) ;
ft_timestamp_from_frame_index = (1/fps) * (0:(frame_count-1)) ;  
  % The timestamps implied by the FT calibration can be slightly different than
  % the actual one from the stamp_log_cam<n>.txt file.
timestamps = ft_timestamp_from_frame_index(is_pfly_on_from_frame_index) ;
dt = diff(timestamps) ;

%pfly_theta_from_frame_index = atan2(pfly_theta_hat_from_frame_index(2,:), pfly_theta_hat_from_frame_index(1,:)) ;
pfly_heading_angle_from_frame_index = atan2(pfly_heading_hat_from_frame_index(2,:), pfly_heading_hat_from_frame_index(1,:)) ;

% Define a function to map from pfly indices to scalar trx structs
function trx = trx_from_pfly_index(pfly_index)
  arena_index = pfly_index ;

  nframes = sum(is_pfly_on_from_frame_index) ;
  x = pfly_position_from_frame_index_from_arena_index(1, is_pfly_on_from_frame_index, arena_index) ;
  y = pfly_position_from_frame_index_from_arena_index(2, is_pfly_on_from_frame_index, arena_index) ;
  theta = pfly_heading_angle_from_frame_index(is_pfly_on_from_frame_index) ;
  a = repmat(fake_fly_a, [1 nframes]) ;
  b = repmat(fake_fly_b, [1 nframes]) ;
  xwingl = nan(1, nframes) ;
  ywingl = nan(1, nframes) ;
  xwingr = nan(1, nframes) ;
  ywingr = nan(1, nframes) ;
  mm_per_pixel = 1/pixels_per_mm ;
  x_mm = mm_per_pixel * x ;
  y_mm = mm_per_pixel * y ;
  a_mm = repmat(mm_per_pixel * fake_fly_a, [1 nframes]) ;
  b_mm = repmat(mm_per_pixel * fake_fly_b, [1 nframes]) ;
  theta_mm = theta ;
  id = maximum_real_fly_id + pfly_index ;

  trx = ...
    struct('moviename', {moviename}, ...
           'moviefile', {moviefile}, ...
           'firstframe', {firstframe}, ...
           'off', {off}, ...
           'endframe', {endframe}, ...
           'nframes', {nframes}, ...
           'fps', {fps}, ...
           'pxpermm', {pxpermm}, ...
           'id', {id}, ...
           'sex', {sex}, ...
           'timestamps', {timestamps}, ...
           'dt', {dt}, ...
           'x', {x}, ...
           'y', {y}, ...
           'theta', {theta}, ...
           'a', {a}, ...
           'b', {b}, ...
           'xwingl', {xwingl}, ...
           'ywingl', {ywingl}, ...
           'xwingr', {xwingr}, ...
           'ywingr', {ywingr}, ...
           'x_mm', {x_mm}, ...
           'y_mm', {y_mm}, ...
           'theta_mm', {theta_mm}, ...
           'a_mm', {a_mm}, ...
           'b_mm', {b_mm}, ...
           'is_pfly', {true}) ;
end  % function

% Use the just-defined function to create the fake_trx struct array
pfly_index_from_arena_index = (1:arena_count) ;
fake_trx = arrayfun(@trx_from_pfly_index, pfly_index_from_arena_index) ;

end

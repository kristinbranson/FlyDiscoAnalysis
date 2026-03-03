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
%   look-up table (LUT).  x,y are in mm, in typical cartesian coordinates (origin at
%   center of arena, x axis points towards right of frame, y axis points
%   towards top of frame), heading angle uses usual math-class convention
%   (angle in radians, relative to x-axis, increasing counterclockwise).
% actual_trajectory: 1st col a timestamp (in seconds); 2nd col the LUT index;
%   3rd col 0 or 1, depending on whether the pfly ("projected fly") is visible or not
% flytracker_calibration: The "calib" structure from a FlyTracker calibration
%   file.
% fake_fly_params: table or struct containing 'maj_axis' and 'min_axis'
%   columns/fields that give the fake fly major/minor axis length.  This is
%   the full length, *not* the semi- axis length or quarter- axis length.  In
%   mm.
% maximum_real_fly_id: A fly "id" is a nonnegative (possibly strictly
%   positive?) integer that uniquely identifies that fly.  This gives the largest 
%   fly id of all the real flies.  The smallest fly id of the projected flies will be 
%   this number plus one.
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

% ellipseTrajectory.csv contains a look-up table from indices to <x,y>
% actualTraj.csv contains a LUT index for each frame
dimension_count = 2 ;
lut_count = size(ellipse_trajectory,1) ;  %#ok<NASGU>  % Number of entries in the LUT
pfly_position_in_mm_from_lut_index = [1 -1]' .* ellipse_trajectory(:,1:dimension_count)' ;  % mm, 2 x lut_count
  % pfly_position_from_lut_index is 2 x number of LUT entries, uses convention that
  % origin is at center of arena, x axis points to right of frame, y axis points
  % to *bottom* of frame.  (We're going to call this "CRTesian coordinates".)

% Plot the x-y coordinates of the pfly at each LUT row  
if do_debug ,
  figure('color', 'w') ;
  plot(pfly_position_in_mm_from_lut_index(1,:), pfly_position_in_mm_from_lut_index(2,:), 'k.') ;
  axis equal ;
  set(gca, 'YDir', 'reverse') ;
  title('LUT (from ellipseTrajectory.csv) x-y values') ;
  xlabel('x (mm)') ;
  ylabel('y (mm)') ;
  % For the original pfly protocol, looks like a circle, early points are on the
  % x-axis and points make a CCW circle.  But later experiments do different
  % things.
end

% % Show a movie of the pfly positions, as a function of LUT index
% if do_debug
%   fig = figure('color', 'w') ;
%   ax = axes('Parent', fig) ;
%   plot(ax, pfly_position_in_mm_from_lut_index(1,:), pfly_position_in_mm_from_lut_index(2,:), '.', 'Color', 0.9*[1 1 1]) ;
%   axis(ax, 'equal') ;
%   set(ax, 'YDir', 'reverse') ;
%   title(ax, 'LUT (from ellipseTrajectory.csv) x-y values') ;
%   xlabel(ax, 'x (mm)') ;
%   ylabel(ax, 'y (mm)') ;
%   drawnow() ;
% 
%   set(ax, 'XLimMode', 'manual', 'YLimMode', 'manual') ;
%   dot = line('Parent', ax, 'XData', [], 'YData', [], 'LineStyle', 'none', 'Marker', '.', 'MarkerFaceColor', 'k') ;
%   for i = 1 : lut_count
%     set(dot, 'XData', pfly_position_in_mm_from_lut_index(1,i), 'YData', pfly_position_in_mm_from_lut_index(2,i)) ;
%     drawnow();
%     pause(0.01) ;
%   end
% end

% Get the angle of the pfly tail-to-head vector for each LUT index.
% Note that ellipse_trajectory(lut_index,3) is an angle, in radians, relative to the x
% axis, increasing counterclockwise.
pfly_heading_angle_from_lut_index = -ellipse_trajectory(:,3)' ;
  % pfly_heading_angle_from_lut_index is 1 x lut_count, angle in radians,
  % relative to x axis, increasing *clockwise*. (to be consistent with CRTesian
  % coordinates).

% Compute the unit vector describing the fly tail-to-head vector at each LUT
% index.  This is in CRTesian coordinates (y-axis points to bottom of frame).
pfly_heading_hat_from_lut_index = [ cos(pfly_heading_angle_from_lut_index) ; sin(pfly_heading_angle_from_lut_index) ] ;

% How many frames in the video
frame_count = numel(timestamp_from_frame_index) ;  %#ok<NASGU> 

% Extract individual signals frm the actual_trajectory arrays
timestamp_from_tick_index = actual_trajectory(:,1)' - actual_trajectory(1,1) ;  % s
lut_index_from_tick_index = actual_trajectory(:,2)' ;
is_pfly_on_from_tick_index = actual_trajectory(:,3)' ;
tick_count = numel(timestamp_from_tick_index) ;  % Usually many more ticks than frames
  % we use "tick" here instead of "frame" just because there are more of them
  % than the frame_count

% Plot some stuff to sanity-check the timestamps
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

% Extract the centroid position and heading angle (as a unit vector) for each tick.
pfly_position_in_mm_from_tick_index = pfly_position_in_mm_from_lut_index(:,lut_index_from_tick_index) ;  % 2 x tick_count
pfly_heading_hat_from_tick_index = pfly_heading_hat_from_lut_index(:,lut_index_from_tick_index) ;  % 2 x tick_count

% Show a movie of the pfly positions, as a function of tick index
if do_debug
  fig = figure('color', 'w') ;
  ax = axes('Parent', fig) ;
  plot(ax, pfly_position_in_mm_from_tick_index(1,:), pfly_position_in_mm_from_tick_index(2,:), '.', 'Color', 0.9*[1 1 1]) ;
  axis(ax, 'equal') ;
  set(ax, 'YDir', 'reverse') ;
  title(ax, 'pFly positions') ;
  xlabel(ax, 'x (mm)') ;
  ylabel(ax, 'y (mm)') ;
  drawnow() ;

  set(ax, 'XLimMode', 'manual', 'YLimMode', 'manual') ;
  dot = line('Parent', ax, 'XData', [], 'YData', [], 'LineStyle', 'none', 'Marker', '.', 'MarkerFaceColor', 'k') ;
  for i = 1 : round(tick_count/10)
    set(dot, 'XData', pfly_position_in_mm_from_tick_index(1,i), 'YData', pfly_position_in_mm_from_tick_index(2,i)) ;
    drawnow();
  end
end

% Plot the x-y coordinates of the pfly at each tick, with time coded by color
if do_debug ,
  desired_sample_count = 200 ;
  stride = round(tick_count / desired_sample_count) ;
  tick_index_from_sample_index = 1:stride:tick_count ;
  % tick_index_from_sample_index = 1:10:450 ;  
  sample_count = numel(tick_index_from_sample_index) ;
  color_from_sample_index = turbo(sample_count) ;
  fig = figure('color', 'w') ;
  ax = axes('Parent', fig, 'DataAspectRatio', [1 1 1], 'YDir', 'reverse') ;
  line('Parent', ax, 'XData', pfly_position_in_mm_from_tick_index(1,:), 'YData', pfly_position_in_mm_from_tick_index(2,:), ...
       'LineStyle', 'none', 'Marker', '.', 'Color', 0.9*[1 1 1]) ;

  for sample_index = 1 : sample_count 
    tick_index = tick_index_from_sample_index(sample_index) ;
    pfly_position_in_mm = pfly_position_in_mm_from_tick_index(:,tick_index) ;
    pfly_heading_hat = pfly_heading_hat_from_tick_index(:,tick_index) ;
    heading_scale_factor = 1 ;    
    clr = color_from_sample_index(sample_index,:) ;
    line('Parent', ax, 'Color', clr, ...
         'XData', pfly_position_in_mm(1)+[0 heading_scale_factor*pfly_heading_hat(1)], ...
         'YData', pfly_position_in_mm(2)+[0 heading_scale_factor*pfly_heading_hat(2)]) ;    
    line('Parent', ax, 'LineStyle', 'none', 'Marker', 'o', 'MarkerSize', 4, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', clr, ...
         'XData', pfly_position_in_mm(1), 'YData', pfly_position_in_mm(2)) ;
  end
  title(ax, 'pFly positions (blue->early, red->late)') ;
  xlabel(ax, 'x (mm)') ;
  ylabel(ax, 'y (mm)') ;
end

% For frames where the fly is moving at a decent speed, plot the heading
% vector angle vs the angle of the velocity vector.  These should agree most of
% the time.
if do_debug
  pfly_dxy_in_mm_from_intertick_index = diff(pfly_position_in_mm_from_tick_index, [], 2) ;  % 2 x (tick_count-1)
  pfly_ds_in_mm_from_intertick_index = sqrt(sum(pfly_dxy_in_mm_from_intertick_index.^2,1)) ;  % mm, step size from tick to tick
  ds_in_mm_threshold = 0.03 ;  % mm/tick, determined empirically
  is_fast_enough_from_intertick_index = (pfly_ds_in_mm_from_intertick_index >= ds_in_mm_threshold) ;
  pfly_heading_hat_from_intertick_index = (pfly_heading_hat_from_tick_index(:,1:end-1) + pfly_heading_hat_from_tick_index(:,2:end))/2 ;
  pfly_heading_angle_from_intertick_index = atan2(pfly_heading_hat_from_intertick_index(2,:), pfly_heading_hat_from_intertick_index(1,:)) ;
  pfly_dxy_angle_from_intertick_index = atan2(pfly_dxy_in_mm_from_intertick_index(2,:), pfly_dxy_in_mm_from_intertick_index(1,:)) ;
  
  figure('color', 'w');
  plot(pfly_heading_angle_from_intertick_index(is_fast_enough_from_intertick_index), ...
       pfly_dxy_angle_from_intertick_index(is_fast_enough_from_intertick_index), ...
       'k.') ;
  xlabel('Heading angle (rad)') ;
  ylabel('Velocity angle (rad)') ;
end

% Interpolate to get pfly position and heading at each movie frame timestamp.
% Also transpose so that rows correspond to frames.
pfly_position_in_mm_from_frame_index = interp1(timestamp_from_tick_index', pfly_position_in_mm_from_tick_index', timestamp_from_frame_index')' ;
  % 2 x frame_count
pfly_heading_hat_from_frame_index = interp1(timestamp_from_tick_index', pfly_heading_hat_from_tick_index', timestamp_from_frame_index')' ;
  % 2 x frame_count 
fractional_is_pfly_on_from_frame_index = interp1(timestamp_from_tick_index', is_pfly_on_from_tick_index', timestamp_from_frame_index')' ;
is_pfly_on_from_frame_index = (fractional_is_pfly_on_from_frame_index>0.5) ;  % 1 x frame_count

% Plot the the per-frame quantities to check them
if do_debug ,
  figure('color', 'w') ;
  line_handles = ...
    plot(timestamp_from_frame_index, pfly_position_in_mm_from_frame_index(1,:), 'r', ...
         timestamp_from_frame_index, pfly_position_in_mm_from_frame_index(2,:), 'b') ;
  xlabel('Movie timestamp (s)') ;
  ylabel('Coordinate (mm)') ;
  legend(line_handles, {'x', 'y'}, 'location', 'northeast') ;
  title('Interpolated position vs movie timestamp') ;
  xlim([800 850]) ;  % zoom in on a random part---too much data otherwise
end

% Extract the flytracker info about the arena locations and size
raw_center_from_arena_index = flipud((flytracker_calibration.centroids)') ;   % 2x9, pels
arena_count = size(raw_center_from_arena_index, 2) ; 
center_from_arena_index = raw_center_from_arena_index + arena_center_shift ;
pixels_per_mm = flytracker_calibration.PPM ;  % pixels per mm
pfly_position_in_pels_from_frame_index = pixels_per_mm * pfly_position_in_mm_from_frame_index ;  % pels, 2 x frame_count

% Generate the pfly position for each arena, for each frame
big_center_from_arena_index = reshape(center_from_arena_index, [dimension_count 1 arena_count]) ;
pfly_position_from_frame_index_from_arena_index = big_center_from_arena_index + pfly_position_in_pels_from_frame_index ;  % 2 x frame_count x arena_count

% Extract things we need from the fake-fly config
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
timestamps = ft_timestamp_from_frame_index(firstframe:endframe) ;
dt = diff(timestamps) ;

% Convert the per-frame headings from unit vectors to angles in radians
% (angles increasing clockwise from the x axis).
pfly_heading_angle_from_frame_index = atan2(pfly_heading_hat_from_frame_index(2,:), pfly_heading_hat_from_frame_index(1,:)) ;

% Define a function to map from pfly indices to scalar trx structs
function trx = trx_from_pfly_index(pfly_index)
  arena_index = pfly_index ;

  nframes = endframe - firstframe + 1 ;
  x = pfly_position_from_frame_index_from_arena_index(1, firstframe:endframe, arena_index) ;
  y = pfly_position_from_frame_index_from_arena_index(2, firstframe:endframe, arena_index) ;
  theta = pfly_heading_angle_from_frame_index(firstframe:endframe) ;
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
  is_pfly_on = is_pfly_on_from_frame_index(firstframe:endframe) ;

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
           'is_pfly', {true}, ...
           'is_pfly_on', is_pfly_on) ;
end  % function

% Use the just-defined function to create the fake_trx struct array
pfly_index_from_arena_index = (1:arena_count) ;
fake_trx = arrayfun(@trx_from_pfly_index, pfly_index_from_arena_index) ;

end  % function

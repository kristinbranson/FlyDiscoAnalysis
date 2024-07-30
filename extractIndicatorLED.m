function indicatordata = extractIndicatorLED(expdir, dataloc_params, indicator_params, ledIndicatorPoints, has_explicit_masks, mask_from_mask_index, debug)

% Deal with optional args
if ~exist('debug', 'var') || isempty(debug) ,
  debug = false ;
end

% Read the video metadata, make a frame-reading function 
moviefile = fullfile(expdir,dataloc_params.moviefilestr);
[readfcn,~,~,headerinfo] = get_readframe_fcn(moviefile);
nr = headerinfo.nr ;
nc = headerinfo.nc ;
nframes = headerinfo.nframes ;
timestamps = headerinfo.timestamps(:)' ;  % want row vector

% Compute some things we'll need for extracting the LED signals from each
% movie frame.
% If we're using masks, that means to compute the overlap of the mask and the
% movie frame.
% If not using masks, things are simpler.
% (Values in ledIndicatorPoints are guaranteed to be integral.)
if has_explicit_masks ,
  mask_width = size(mask_from_mask_index,2) ;
  mask_height = size(mask_from_mask_index,1) ;
  mask_count = size(mask_from_mask_index,3) ;
else
  ledwindowr = indicator_params.indicatorwindowr ;
  mask_height = 2*ledwindowr + 1 ;
  mask_width = 2*ledwindowr + 1 ;
  mask_count = 1 ;
end
mask_shape = [ mask_height mask_width ] ;
frame_shape = [ nr nc ] ;  % yx order
% There are three rectangles to deal with: the frame, the mask, and the
% trimmed mask.  The center of the mask in frame coordinates is given by
% ledIndicatorPoints, which are all integers.  The center of the mask in mask
% coordinates is the middle of the mask.  The trimmed mask is the part of the
% mask that overlaps with the frame.  In what follows, the bounds of an image
% in x/y, in a particular reference frame, are given by a "bound", a
% two-element row vector in yx order.  The bounds endpoints are inclusive.  An
% image has a lower bound and an upper bound.  We also talk about "anchors",
% which are two-element row vectors in yx order that determine the location of
% something.  The center of an image is one kind of anchor, but an anchor is
% not always a center.  They are often the center of that thing, but not
% always.  We also deal with "extents".  The rule is that anchor + extent =
% bound.  And image has a lower extent and an upper extent.  Bounds, anchors,
% and centers are all relative to a coordinate system, but extents are not.
% Lower extents are typically all-negative, upper extents are typically all-postive.
% All elements of anchors, bounds, and extents are integers.
% The lower extent of an image with odd dimensions will be the negative of the
% upper extent, when the image center is the anchor.
mask_center_in_mask_coords = [ round((mask_height+1)/2) round((mask_width+1)/2) ] ;  % yx order
mask_hi_extent = mask_shape - mask_center_in_mask_coords ;  % yx order
mask_lo_extent = [1 1] - mask_center_in_mask_coords ;  % yx order
mask_center_in_frame_coords = [ ledIndicatorPoints(2) ledIndicatorPoints(1) ] ;  % yx order
mask_lo_bound_in_frame_coords = mask_center_in_frame_coords + mask_lo_extent ;
mask_hi_bound_in_frame_coords = mask_center_in_frame_coords + mask_hi_extent ;
trimmed_mask_lo_bound_in_frame_coords = max([1 1], mask_lo_bound_in_frame_coords) ;
trimmed_mask_hi_bound_in_frame_coords = min(mask_hi_bound_in_frame_coords, frame_shape) ;
trimmed_mask_lo_extent = trimmed_mask_lo_bound_in_frame_coords - mask_center_in_frame_coords ;
trimmed_mask_hi_extent = trimmed_mask_hi_bound_in_frame_coords - mask_center_in_frame_coords ;
trimmed_mask_shape = trimmed_mask_hi_extent - trimmed_mask_lo_extent + [1 1] ;
trimmed_mask_anchor_in_frame_coords = mask_center_in_frame_coords ;
trimmed_mask_anchor_in_mask_coords = mask_center_in_mask_coords ;
trimmed_mask_lo_bound_in_mask_coords = trimmed_mask_anchor_in_mask_coords + trimmed_mask_lo_extent ;
trimmed_mask_hi_bound_in_mask_coords = trimmed_mask_anchor_in_mask_coords + trimmed_mask_hi_extent ;
trimmed_mask_lo_bound_in_frame_coords = trimmed_mask_anchor_in_frame_coords + trimmed_mask_lo_extent ;
trimmed_mask_hi_bound_in_frame_coords = trimmed_mask_anchor_in_frame_coords + trimmed_mask_hi_extent ;
if has_explicit_masks ,
  trimmed_mask_from_mask_index = ...
    mask_from_mask_index(trimmed_mask_lo_bound_in_mask_coords(1):trimmed_mask_hi_bound_in_mask_coords(1), ...
                         trimmed_mask_lo_bound_in_mask_coords(2):trimmed_mask_hi_bound_in_mask_coords(2), ...
                         :) ;
else
  trimmed_mask_from_mask_index = true(trimmed_mask_shape) ;
end

% Plot the middle frame along with the area(s) to be used for signal
% extraction
if debug ,
  frame = readfcn(min(max(1,round(nframes/2)),nframes)) ;
  f = figure('color', 'w') ;
  a = axes(f) ;
  imshow(frame, 'Parent', a) ;
  a.YDir = 'reverse' ;
  a.Colormap = parula(256) ;
  corners = [ trimmed_mask_lo_bound_in_frame_coords(2)-0.5 trimmed_mask_lo_bound_in_frame_coords(1)-0.5 ; ...
              trimmed_mask_hi_bound_in_frame_coords(2)+0.5 trimmed_mask_lo_bound_in_frame_coords(1)-0.5 ; ...
              trimmed_mask_hi_bound_in_frame_coords(2)+0.5 trimmed_mask_hi_bound_in_frame_coords(1)+0.5 ; ...
              trimmed_mask_lo_bound_in_frame_coords(2)-0.5 trimmed_mask_hi_bound_in_frame_coords(1)+0.5 ; ...
              trimmed_mask_lo_bound_in_frame_coords(2)-0.5 trimmed_mask_lo_bound_in_frame_coords(1)-0.5 ] ; 
  x_data = corners(:,1)' ;
  y_data = corners(:,2)' ;
  line('Parent', a, 'XData', x_data, 'YData', y_data, 'Color', 'w', 'LineWidth', 3) ;
  line('Parent', a, 'XData', x_data, 'YData', y_data, 'Color', 'r', 'LineWidth', 1) ;  % 0.5 is the default line width
  if has_explicit_masks ,
    for mask_index = 1 : mask_count ,
      trimmed_mask = trimmed_mask_from_mask_index(:,:,mask_index) ;
      boundary_from_island_index = bwboundaries(trimmed_mask) ;  % boundary_count x 2, in yx order
      boundary = boundary_from_island_index{1} ;  % hopefully there's only one
      boundary_mean = mean(boundary, 1) ; % 1x2, yx order
      line('Parent', a, ...
           'XData', trimmed_mask_lo_bound_in_frame_coords(2)+1+boundary(:,2), ...
           'YData', trimmed_mask_lo_bound_in_frame_coords(1)-1+boundary(:,1), ...
           'Color', 'w', ...
           'LineWidth', 3) ;
      line('Parent', a, ...
           'XData', trimmed_mask_lo_bound_in_frame_coords(2)+1+boundary(:,2), ...
           'YData', trimmed_mask_lo_bound_in_frame_coords(1)-1+boundary(:,1), ...
           'Color', 'r', ...
           'LineWidth', 1) ;
      text(a, ...
           trimmed_mask_lo_bound_in_frame_coords(2)+1+boundary_mean(2), ...
           trimmed_mask_lo_bound_in_frame_coords(1)+1+boundary_mean(1), ...
           sprintf('%d', mask_index), ...
           'Color', 'w', ...
           'FontWeight', 'bold', ...
           'HorizontalAlignment', 'center') ;      
      text(a, ...
           trimmed_mask_lo_bound_in_frame_coords(2)+1+boundary_mean(2), ...
           trimmed_mask_lo_bound_in_frame_coords(1)+1+boundary_mean(1), ...
           sprintf('%d', mask_index), ...
           'Color', 'r', ...
           'HorizontalAlignment', 'center') ;      
    end
  end
  a.Box = 'on' ;
  a.Layer = 'top' ;
  title(a, 'Example frame and LED detection window', 'Fontsize', 11) ;
  drawnow() ;
end

% Figure out how to combine all the pixels in each region/mask, in each frame
% Will use either max or mean
if isfield(indicator_params, 'reducefun') ,
  reducefun = indicator_params.reducefun ;
  if strcmp(reducefun, 'mean') ,
    do_reduce_with_max = false ;
  elseif strcmp(reducefun, 'max') ,
    do_reduce_with_max = true ;
  else
    if ischar(reducefun)
      error('indicator_params.reducefun is ''%s''.  It should be ''mean'' or ''max''.', reducefun) ;
    else
      error('indicator_params.reducefun is of class %s.  It should be a char array with value ''mean'' or ''max''.', class(reducefun)) ;
    end
  end
else
  % Use mean, for backwards-compatibility
  do_reduce_with_max = false ;
end

% Read in video frames, extract the LED indicator region, compute the mean+max
% of that region for each frame
pixel_max_signal = zeros(mask_count,nframes) ;
pixel_mean_signal = zeros(mask_count,nframes) ;
fprintf('Extracting LED indicator signal from video frames...\n') ;
pbo = progress_bar_object(nframes) ;
for frame_index = 1:nframes
  frame = readfcn(frame_index);
  detail = double(frame(trimmed_mask_lo_bound_in_frame_coords(1):trimmed_mask_hi_bound_in_frame_coords(1), ...
                        trimmed_mask_lo_bound_in_frame_coords(2):trimmed_mask_hi_bound_in_frame_coords(2))) ;
  masked_detail = detail .* trimmed_mask_from_mask_index ;  % detail_height x detail_width x mark_count
  pixel_max_signal(:,frame_index) = reshape(max(masked_detail, [], [1 2]), ...
                                            [mask_count 1]) ;
  pixel_mean_signal(:,frame_index) = reshape(mean(masked_detail, [1 2]), ...
                                             [mask_count 1]) ;
  if (mod(frame_index,1000) == 0) ,    
    pbo.update(1000) ;
  end
end
pbo.update(inf) ;
fprintf('Done extracting LED indicator signal from video frames.\n') ;

% Make some debug plots, if called for
if debug ,
  % Plot the raw extracted signals
  f = figure('color', 'w') ;
  set_figure_size([10 6]) ;
  a = axes(f) ;
  plot(a, 1:nframes, pixel_mean_signal) ;
  xlabel(a, 'Frame index') ;
  ylabel(a, 'LED mean signal (counts)') ;
  title(a, 'LED mean signal', 'Fontsize', 11) ;
  drawnow() ;  

  f = figure('color', 'w') ;
  set_figure_size([10 6]) ;
  a = axes(f) ;
  plot(a, 1:nframes, pixel_max_signal) ;
  xlabel(a, 'Frame index') ;
  ylabel(a, 'LED max signal (counts)') ;
  title(a, 'LED max signal', 'Fontsize', 11) ;
  drawnow() ;    
end

% Select which signal to go forward with
if do_reduce_with_max ,
  analog_led_signal = pixel_max_signal ;
else
  analog_led_signal = pixel_mean_signal ;
end

% Report the range, standard deviation of the LED signal.
analog_led_signal_max = max(analog_led_signal, [], 2) ;
analog_led_signal_min = min(analog_led_signal, [], 2) ;
analog_led_signal_range = analog_led_signal_max - analog_led_signal_min ;
analog_led_signal_sd = std(analog_led_signal, [], 2) ;
for mask_index = 1 : mask_count ,
  fprintf('LED signal %d range is %g counts, standard deviation is %g counts\n', ...
          mask_index, analog_led_signal_range(mask_index), analog_led_signal_sd(mask_index)) ;
end

% Error if any of the LED signals have very low range
minimum_signal_range = 5 ;
for mask_index = 1 : mask_count ,
  if analog_led_signal_range(mask_index) < minimum_signal_range ,
    fprintf('Insufficient variation in LED signal %d to reliably detect edges\n', mask_index) ;
  end
end
if any(analog_led_signal_range < minimum_signal_range) ,
    error('Insufficient variation in at least one LED signal to reliably detect edges') ;
end

% % Find start and ends of light stimulus.
%
% For pulsed stimuli, group pulses together by mean(pad length) is either 1/0 before/after start/stop.
% Pad needs to be more than 1/2 duty cycle of pulses / sampling frequency.
% Offtime can't be less than pad+1 or it will fail.
analog_led_threshold = 0.6*analog_led_signal_max + 0.4*analog_led_signal_min ;  % mask_count x 1
  % Changed to 60/40 to help with AO pfly exps in summer 2024.  Only needed b/c
  % we're not properly reading all three flydisco LEDs properly, but taking a
  % signal that gets one corner of all three.  They're using a 2nd LED to
  % indicate when the pfly is shown, which is on for most frames except possible
  % a few at the beginning and end.  using a 60/40 threshold usually allows us
  % to extract the optogenetic stimulus LED even if this other LED is active.
raw_binary_led_signal = (analog_led_signal>analog_led_threshold) ;  % mask_count x nframes
mask_count = size(raw_binary_led_signal, 1) ;
if mask_count==1 ,
  [binary_led_signal, startframe, endframe] = compute_pulse_envelope(raw_binary_led_signal, indicator_params.pad) ;
else
  [binary_led_signal, startframe, endframe] = compute_multiple_pulse_envelopes(raw_binary_led_signal, indicator_params.pad) ;
end

% For debugging, compare indicatordigital to thresholded meanimage
if debug ,
  for mask_index = 1 : mask_count ,
    f = figure('color', 'w') ;
    set_figure_size([10 6]) ;
    a1 = axes(f, 'Position', [0.1 0.6 0.8 0.3]) ;
    a2 = axes(f, 'Position', [0.1 0.1 0.8 0.3]) ;
    l1 = line('Parent', a1, 'XData', 1:nframes, 'YData', analog_led_signal(mask_index,:), 'Color', [0 0.8 0]) ;
    l2 = line('Parent', a1, 'XData', [1 nframes], 'YData', analog_led_threshold(mask_index)*[1 1], 'Color', [242 140 40]/255, 'LineStyle', '--') ;  
    l3 = line('Parent', a2, 'XData', 1:nframes, 'YData', double(raw_binary_led_signal(mask_index,:)), 'Color', 'b') ;
    l4 = line('Parent', a2, 'XData', 1:nframes, 'YData', double(binary_led_signal(mask_index,:)), 'Color', 'r') ;
    xlabel(a2, 'Frame index') ;
    ylabel(a1, 'LED signal (counts)') ;
    ylabel(a2, 'LED binary') ;
    ylim(a2, [-0.05 1.05]) ;
    legend(a1, [l1 l2], {'signal', 'threshold'}) ;
    legend(a2, [l3 l4], {'signal', 'envelope'}) ;
    title(a1, sprintf('LED %d signal, binary, and envelope', mask_index), 'Fontsize', 11) ;
    drawnow() ;
  end
end

% Put stuff into indicatorLED
indicatorLED = [] ;
indicatorLED.startframe = startframe ;
indicatorLED.endframe = endframe ;
indicatorLED.StartEndStatus = ( [analog_led_signal(1) analog_led_signal(end)] > analog_led_threshold ) ;
indicatorLED.indicatordigital = binary_led_signal ;
if iscell(startframe) ,
  indicatorLED.starttimes = cellfun(@(indices)(timestamps(indices)), startframe, 'UniformOutput', false) ;
  indicatorLED.endtimes = cellfun(@(indices)(timestamps(indices)), endframe, 'UniformOutput', false) ;
else
  indicatorLED.starttimes = timestamps(startframe) ;  % row vector
  indicatorLED.endtimes = timestamps(endframe) ;  % row vector
end

% Put stuff, including indicatorLED, into indicatordata
indicatordata = struct() ;
indicatordata.rect = [trimmed_mask_lo_bound_in_frame_coords(2) trimmed_mask_lo_bound_in_frame_coords(1) trimmed_mask_shape(2) trimmed_mask_shape(1)] ;
indicatordata.has_explicit_masks = has_explicit_masks ;
if has_explicit_masks ,
  indicatordata.trimmed_mask_lo_bound_in_frame_coords = trimmed_mask_lo_bound_in_frame_coords ;
  indicatordata.trimmed_mask_hi_bound_in_frame_coords = trimmed_mask_hi_bound_in_frame_coords ;
  indicatordata.trimmed_mask_from_mask_index = trimmed_mask_from_mask_index ;
end
indicatordata.maximage = pixel_max_signal ;
indicatordata.meanimage = pixel_mean_signal ;
indicatordata.analog_led_signal = analog_led_signal ;
indicatordata.indicatorLED = indicatorLED ;

end

function indicatordata = extractIndicatorLED(expdir, dataloc_params, indicator_params, ledIndicatorPoints, debug)

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

% Compute the rectangle that will be used for the LED signal
ledwindowr = indicator_params.indicatorwindowr;
x = max(1,ledIndicatorPoints(1)-ledwindowr);
y = max(1,ledIndicatorPoints(2)-ledwindowr);
height = min(ledwindowr*2,nr-y);
width = min(ledwindowr*2,nc-x);

% Read in video frames, extract the LED indicator region, compute the mean+max
% of that region for each frame
maximage = zeros(1,nframes) ;
meanimage = zeros(1,nframes);
fprintf('Extracting LED indicator signal from video frames...\n') ;
pbo = progress_bar_object(nframes) ;
for i = 1:nframes
  frame = readfcn(i);
  detail = frame(y:y+height,x:x+width);
  maximage(i) = max(detail(:));
  meanimage(i) = mean(detail(:));
  if (mod(i,1000) == 0) ,
    pbo.update(1000) ;
    %disp(round((i/nframes)*100))
  end
end
pbo.update(inf) ;
fprintf('Done extracting LED indicator signal from video frames.\n') ;

% plot the middle frame along with the area used for LED detection
if debug ,
  frame = readfcn(min(max(1,round(nframes/2)),nframes)) ;
  f = figure('color', 'w') ;
  a = axes(f) ;
  imshow(frame, 'Parent', a) ;
  a.YDir = 'reverse' ;
  a.Colormap = parula(256) ;
  corners = [x       y        ; ...
    x+width y        ; ...
    x+width y+height ; ...
    x       y+height ; ...
    x       y        ] ;
  x_data = corners(:,1)' ;
  y_data = corners(:,2)' ;
  line('Parent', a, 'XData', x_data, 'YData', y_data, 'Color', 'r') ;
  a.Box = 'on' ;
  a.Layer = 'top' ;
  title(a, 'Example frame and LED detection window', 'Fontsize', 11) ;

  % For debugging, plot the derived LED signals
  f = figure('color', 'w') ;
  set_figure_size([10 6]) ;
  a1 = axes(f, 'Position', [0.1 0.6 0.8 0.3]) ;
  a2 = axes(f, 'Position', [0.1 0.1 0.8 0.3]) ;
  line('Parent', a1, 'XData', 1:nframes, 'YData', meanimage, 'Color', 'r') ;
  line('Parent', a2, 'XData', 1:nframes, 'YData', maximage, 'Color', 'b') ;
  xlabel(a2, 'Frame index') ;
  ylabel(a1, 'Mean') ;
  ylabel(a2, 'Max') ;
  title(a1, 'Mean and max of LED window for each frame', 'Fontsize', 11) ;
  drawnow() ;
end

% Report the standard deviation of the LED signal.
% Error if it's crazy low.
%meanimage_range = max(meanimage)-min(meanimage) ;
meanimage_max = quantile(meanimage,0.999) ;  % quantiles more robust than max/min
meanimage_min = quantile(meanimage,0.001) ;
meanimage_range = meanimage_max-meanimage_min ;  % more robust than max/min
meanimage_sd = std(meanimage) ;
fprintf('LED signal range is %g counts, standard deviation is %g counts\n', meanimage_range, meanimage_sd) ;
if meanimage_range < 5 ,
  error('Insufficient variation in LED signal to reliably detect edges') ;
end

% % Find start and ends of light stimulus.
%
% For pulsed stimuli, group pulses together by mean(pad length) is either 1/0 before/after start/stop.
% Pad needs to be more than 1/2 duty cycle of pulses / sampling frequency.
% Offtime can't be less than pad+1 or it will fail.
IRthreshold = 0.6*meanimage_max + 0.4*meanimage_min ;  
  % Changed to 60/40 to help with AO pfly exps in summer 2024.  Only needed b/c
  % we're not properly reading all three flydisco LEDs properly, but taking a
  % signal that gets one corner of all three.  They're using a 2nd LED to
  % indicate when the pfly is shown, which is on for most frames except possible
  % a few at the beginning and end.  using a 60/40 threshold usually allows us
  % to extract the optogenetic stimulus LED even if this other LED is active.
meanimagethresh = (meanimage>IRthreshold) ;
[indicatordigital, startframe, endframe] = compute_pulse_envelope(meanimagethresh, indicator_params.pad) ;

% For debugging, compare indicatordigital to thresholded meanimage
if debug ,
  f = figure('color', 'w') ;
  set_figure_size([10 6]) ;
  a1 = axes(f, 'Position', [0.1 0.6 0.8 0.3]) ;
  a2 = axes(f, 'Position', [0.1 0.1 0.8 0.3]) ;
  l1 = line('Parent', a1, 'XData', 1:nframes, 'YData', meanimage, 'Color', [0 0.8 0]) ;
  l2 = line('Parent', a1, 'XData', [1 nframes], 'YData', IRthreshold*[1 1], 'Color', [242 140 40]/255, 'LineStyle', '--') ;  
  l3 = line('Parent', a2, 'XData', 1:nframes, 'YData', double(meanimagethresh), 'Color', 'b') ;
  l4 = line('Parent', a2, 'XData', 1:nframes, 'YData', double(indicatordigital), 'Color', 'r') ;
  xlabel(a2, 'Frame index') ;
  ylabel(a1, 'LED mean (counts)') ;
  ylabel(a2, 'LED binary') ;
  ylim(a2, [-0.05 1.05]) ;
  legend(a1, [l1 l2], {'signal', 'threshold'}) ;
  legend(a2, [l3 l4], {'signal', 'envelope'}) ;
  title(a1, 'LED mean, binary, and envelope', 'Fontsize', 11) ;
  drawnow() ;
end

% Put stuff into indicatorLED
indicatorLED = [] ;
indicatorLED.startframe = startframe ;
indicatorLED.endframe = endframe ;
indicatorLED.StartEndStatus = ( [meanimage(1) meanimage(end)] > IRthreshold ) ;
indicatorLED.indicatordigital = indicatordigital ;
indicatorLED.starttimes = headerinfo.timestamps(startframe)' ;
indicatorLED.endtimes = headerinfo.timestamps(endframe)' ;

% Put stuff, including indicatorLED, into indicatordata
indicatordata = struct() ;
indicatordata.rect = [x y width height] ;
indicatordata.maximage = maximage ;
indicatordata.meanimage = meanimage ;
indicatordata.indicatorLED = indicatorLED ;

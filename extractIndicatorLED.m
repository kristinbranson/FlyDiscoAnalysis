function indicatordata = extractIndicatorLED(expdir, dataloc_params, indicator_params, ledIndicatorPoints)

% % Load in the LED indicator points
% registrationdatafile = fullfile(expdir,dataloc_params.registrationmatfilestr);
% load(registrationdatafile,'ledIndicatorPoints');

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

% % plot the middle frame along with the area used for LED detection
% frame = readfcn(min(max(1,round(nframes/2)),nframes)) ;
% f = figure('color', 'w') ;
% a = axes(f) ;
% imshow(frame, 'Parent', a) ;
% a.YDir = 'normal' ;  % this is how FlyBowl/FlyBubble/FlyDisco videos are typically viewed, seemingly
% a.Colormap = parula(256) ;
% corners = [x       y        ; ...
%            x+width y        ; ...
%            x+width y+height ; ...
%            x       y+height ; ...
%            x       y        ] ;
% x_data = corners(:,1)' ;
% y_data = corners(:,2)' ;
% line('Parent', a, 'XData', x_data, 'YData', y_data, 'Color', 'r') ;
% a.Box = 'on' ;
% a.Layer = 'top' ;
% 
% % For debugging, plot the derived LED signals
% f = figure('color', 'w') ;
% set_figure_size([10 6]) ;
% a1 = axes(f, 'Position', [0.1 0.6 0.8 0.3]) ;
% a2 = axes(f, 'Position', [0.1 0.1 0.8 0.3]) ;
% line('Parent', a1, 'XData', 1:nframes, 'YData', meanimage, 'Color', 'r') ;
% line('Parent', a2, 'XData', 1:nframes, 'YData', maximage, 'Color', 'b') ;
% xlabel(a2, 'Frame index') ;
% ylabel(a1, 'Mean') ;
% ylabel(a2, 'Max') ;

% Report the standard deviation of the LED signal.
% Error if it's crazy low.
meanimage_range = max(meanimage)-min(meanimage) ;
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
IRthreshold = (max(meanimage)+min(meanimage))/2 ;
meanimagethresh = (meanimage>IRthreshold) ;
[indicatordigital, startframe, endframe] = compute_pulse_envelope(meanimagethresh, indicator_params.pad) ;

% % For debugging, compare indicatordigital to thresholded meanimage
% f = figure('color', 'w') ;
% set_figure_size([10 6]) ;
% a = axes(f) ;
% line('Parent', a, 'XData', 1:nframes, 'YData', double(meanimagethresh), 'Color', 'b') ;
% line('Parent', a, 'XData', 1:nframes, 'YData', double(indicatordigital), 'Color', 'r') ;
% xlabel(a, 'Frame index') ;
% ylabel(a, 'LED') ;
% ylim(a, [-0.05 1.05]) ;

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

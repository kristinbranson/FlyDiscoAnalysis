function [success,msgs] = FlyDiscoDectectIndicatorLedOnOff(expdir,varargin)


success = true;
msgs = [];

version = '0.3';
timestamp = datestr(now,'yyyymmddTHHMMSS');

[analysis_protocol,settingsdir,datalocparamsfilestr] = ...
  myparse(varargin,...
  'analysis_protocol','current_bubble',...
  'settingsdir', internal_settings_folder_path(), ...
  'datalocparamsfilestr','dataloc_params.txt');

%% read in the data locations
datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
dataloc_params = ReadParams(datalocparamsfile);

% Ignore the indicator log file parameter, just log to stdout
logfid = 1 ;
% if isfield(dataloc_params,'indicator_logfilestr'),
%   logfile = fullfile(expdir,dataloc_params.indicator_logfilestr);
%   logfid = fopen(logfile,'a');
%   if logfid < 1,
%     warning('Could not open log file %s\n',logfile);
%     logfid = 1;
%   end
% else
%   logfid = 1;
% end

fprintf(logfid,'\n\n***\nRunning %s version %s analysis_protocol %s at %s\n',mfilename(), version, analysis_protocol, timestamp) ;
%% read in indicator params

indicatorparamsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.indicatorparamsfilestr);
if ~exist(indicatorparamsfile,'file')
  error('Indictor params file %s does not exist',indicatorparamsfile);
end

indicator_params = ReadParams(indicatorparamsfile);
%% read in registration params and check for optogenetic flag

registrationparamsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.registrationparamsfilestr);
if ~exist(registrationparamsfile,'file')
  error('Registration params file %s does not exist',registrationparamsfile);
end

registration_params = ReadParams(registrationparamsfile);

if isfield(registration_params,'OptogeneticExp')
  if registration_params.OptogeneticExp
    DoLEDdetection = true;
    fprintf(logfid,'\n Optogenetic flag set to: %1.0f, running indicator detection',registration_params.OptogeneticExp);
  else
    DoLEDdetection = false;
    fprintf(logfid,'\n Optogenetic flag set to: %1.0f, skipping indicator detection',registration_params.OptogeneticExp)
  end
else
  error('Registration params optogenetic flag %s does not exist',registrationparamsfile);
end


%% make LED window
indicatordata = struct;
%only run for optogentic experiments

if DoLEDdetection
  
  registrationdatafile = fullfile(expdir,dataloc_params.registrationmatfilestr);
  %
  load(registrationdatafile,'ledIndicatorPoints');
  % name of movie file
  moviefile = fullfile(expdir,dataloc_params.moviefilestr);
  ledwindowr = indicator_params.indicatorwindowr;
  [readfcn,~,fid,headerinfo] = get_readframe_fcn(moviefile);
  x = max(1,ledIndicatorPoints(1)-ledwindowr);
  y = max(1,ledIndicatorPoints(2)-ledwindowr);
  height = min(ledwindowr*2,headerinfo.nr-y);
  width = min(ledwindowr*2,headerinfo.nc-x);
  
  
  maximage = zeros(1,headerinfo.nframes);
  meanimage = zeros(1,headerinfo.nframes);
  for i = 1:headerinfo.nframes
    %im = uint8(zeros(width+1,height+1));
    tmp = readfcn(i);
    im = tmp(y:y+height,x:x+width);
    maximage(i) = max(im(:));
    meanimage(i) = mean(im(:));
    if (mod(i,1000) == 0) ,
      disp(round((i/headerinfo.nframes)*100))
    end
    
  end
  
  indicatordata.rect = [x y width height];
  %indicatordata.im = im;
  indicatordata.maximage = maximage;
  indicatordata.meanimage = meanimage;
  if fid > 0
    fclose(fid);
  end
  %% find start and ends of light stimulus.
  % For pulsed stimuli, group pulses together by mean(pad length) is either 1/0 before/after start/stop.
  % Pad needs to be more than 1/2 duty cycle of pulses / sampling frequency.
  % Offtime can't be less than pad+1 or it will fail.
  
  try
    IRthreshold = (max(meanimage)+min(meanimage))/2 ;
    % pad need to be more than 1/2 duty cycle / sampling frequency
    pad = indicator_params.pad;
    indicatorLED = [];
    thresholdimageorig = meanimage > IRthreshold;
    indicatorLED.StartEndStatus = logical([thresholdimageorig(1),thresholdimageorig(end)]);
    
    % pad threshold image
    thresholdimage = zeros(1,numel(meanimage)+2*pad);
    thresholdimage(pad+1:end-pad) = thresholdimageorig;
    %find the start of light pulses
    onpulsecount = 0;
    for i = pad+1:(numel(thresholdimage)-pad)
      if (thresholdimage(i) == 1) && (mean(thresholdimage(i-pad:i-1)) == 0) ,
        onpulsecount = onpulsecount+1;
        indicatorLED.startframe(onpulsecount) = i-pad;
      end
    end
    %find the end of light pulses
    offpulsecount = 0;
    for i = pad+1:numel(thresholdimage)-pad
      if (thresholdimage(i) == 0) && thresholdimage(i-1) == 1 && (mean(thresholdimage(i+1:i+pad)) == 0) ,
        offpulsecount = offpulsecount+1;
        indicatorLED.endframe(offpulsecount) = i-pad;
      end
    end
    
    % make digital verion from start and end - pulses grouped together
    indicatordigital = zeros(1,headerinfo.nframes);
    if ~isempty(indicatorLED.startframe) && ~indicatorLED.StartEndStatus(1) && ~indicatorLED.StartEndStatus(2) % movie starts and ends with indicator off
      for i = 1:numel(indicatorLED.startframe)
        indicatordigital(indicatorLED.startframe(i):indicatorLED.endframe(i)-1) = 1;
      end
      % not tested, I don't have these kind of movies
    elseif ~isempty(indicatorLED.startframe) && indicatorLED.StartEndStatus(1) && indicatorLED.StartEndStatus(2) %movie starts and ends with indicator on
      indicatordigital(1:indicatorLED.endframe(1)-1) = 1;
      for i = 1:numel(indicatorLED.startframe) -1
        j = i+1;
        indicatordigital(indicatorLED.startframe(i):indicatorLED.endframe(j)-1) = 1;
      end
      indicatordigital(indicatorLED.startframe(end):headerinfo.nframes) = 1;
    elseif ~isempty(indicatorLED.startframe) && ~indicatorLED.StartEndStatus(1) && indicatorLED.StartEndStatus(2) %movie starts with indicator off and ends with indicator on
      for i = 1:numel(indicatorLED.startframe)
        indicatordigital(indicatorLED.startframe(i):indicatorLED.endframe(i)-1) = 1;
      end
      indicatordigital(indicatorLED.startframe(end):headerinfo.nframes) = 1;
      
    elseif ~isempty(indicatorLED.startframe) && indicatorLED.StartEndStatus(1) && ~indicatorLED.StartEndStatus(2) %movie starts with indicator on and ends with indicator off
      indicatordigital(1:indicatorLED.endframe(1)-1) = 1;
      for i = 1:numel(indicatorLED.startframe) -1
        j = i+1;
        indicatordigital(indicatorLED.startframe(i):indicatorLED.endframe(j)-1) = 1;
      end
    else
      fprintf(logfid,'\n Error creating indicatordigital');
      msgs{end+1} = 'Error creating indicatordigital' ;
    end
    
    indicatorLED.indicatordigital = indicatordigital;
    indicatorLED.starttimes = headerinfo.timestamps(indicatorLED.startframe)';
    indicatorLED.endtimes = headerinfo.timestamps(indicatorLED.endframe)';
    indicatordata.indicatorLED = indicatorLED;
  catch exception
    imagediff = max(meanimage)-min(meanimage);
    fprintf(logfid, '\n Error calculating indicator on-off: Error identifier was "%s", error message was "%s"', exception.identifier, exception.message) ;
    fprintf(logfid, '\n Possibly relevant: meanimage value is small: %g', imagediff) ;
    success = false ;
    msgs{end+1} = 'Error calculating indicator LED on-off' ;
    return
  end
  
end

% %% plot figure
% figure, plot(meanimage)
% hold on, plot(indicatorLED.startframe,meanimage(indicatorLED.startframe),'rx')
% hold on, plot(indicatorLED.endframe,meanimage(indicatorLED.endframe),'gx')

%% compare to timestamps from file

%% reconstruct protocol

%% save indicator data to mat file
indicatordatamatfile = fullfile(expdir,dataloc_params.indicatordatafilestr);

didsave = false;
try
  if exist(indicatordatamatfile,'file'),
    delete(indicatordatamatfile);
  end
  save(indicatordatamatfile,'-struct','indicatordata');
  didsave = true;
catch ME
  warning('Could not save indicator data to mat file: %s',getReport(ME));
  success = false;
  msgs{end+1} = ['Could not save indicator data to mat file: %s',getReport(ME)];
end

if didsave,
  fprintf(logfid,'\n Saved indicator data to file %s\n',indicatordatamatfile);
else
  fprintf(logfid,'\n Could not save indicator data to mat file:\n%s\n',getReport(ME));
end
%% close log

fprintf(logfid,'\n Finished running %s at %s.\n',mfilename(), datestr(now,'yyyymmddTHHMMSS'));

if logfid > 1,
  fclose(logfid);
end

%% close figures

if isdeployed,
  delete(findall(0,'type','figure'));
end

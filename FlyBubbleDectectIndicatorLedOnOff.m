function [success,msgs] = FlyBubbleDectectIndicatorLedOnOff(expdir,varargin)


success = true;
msgs = [];

version = '0.2';
timestamp = datestr(now,'yyyymmddTHHMMSS');

[analysis_protocol,settingsdir,datalocparamsfilestr,registrationmatfilestr] = ...
  myparse(varargin,...
  'analysis_protocol','20150428_flybubble_centralcomplex',...
  'settingsdir','/groups/branson/home/robiea/Code_versioned/FlyBubbleAnalysis/settings',...
  'datalocparamsfilestr','dataloc_params.txt',...
  'registrationmatfilestr','registrationdata.mat');

%% read in the data locations
datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
dataloc_params = ReadParams(datalocparamsfile);
   
if isfield(dataloc_params,'indicator_logfilestr'),
  logfile = fullfile(expdir,dataloc_params.indicator_logfilestr);
  logfid = fopen(logfile,'a');
  if logfid < 1,
    warning('Could not open log file %s\n',logfile);
    logfid = 1;
  end
else
  logfid = 1;
end

fprintf(logfid,'\n\n***\nRunning FlyBubbleDectectIndicatorLedOnOff version %s analysis_protocol %s at %s\n',version,analysis_protocol,timestamp);
%% read in indicator params

indicatorparamsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.indicatorparamsfilestr);
if ~exist(indicatorparamsfile,'file')
  error('Indictor params file %s does not exist',indicatorparamsfile);
end

indicator_params = ReadParams(indicatorparamsfile);

%% make LED window

registrationdatafile = fullfile(expdir,dataloc_params.registrationmatfilestr);
% 
load(registrationdatafile,'ledIndicatorPoints');
% name of movie file
moviefile = fullfile(expdir,dataloc_params.moviefilestr);
ledwindowr = indicator_params.indicatorwindowr;
[readfcn,nframes,fid,headerinfo] = get_readframe_fcn(moviefile);
x = max(1,ledIndicatorPoints(1)-ledwindowr);
y = max(1,ledIndicatorPoints(2)-ledwindowr);
height = min(ledwindowr*2,headerinfo.nr-y);
width = min(ledwindowr*2,headerinfo.nc-x);


maximage = zeros(1,headerinfo.nframes);
meanimage = zeros(1,headerinfo.nframes);
for i = 1:headerinfo.nframes
    %im = uint8(zeros(width+1,height+1));
    tmp = readfcn(i);
    im = tmp(x:x+width,y:y+height);
    maximage(i) = max(im(:));
    meanimage(i) = mean(im(:));
%     if (mod(i,1000) == 0);
%       disp(round((i/headerinfo.nframes)*100))
%    end
    
end
indicatordata = struct;
indicatordata.rect = [x y width height];
%indicatordata.im = im;
indicatordata.maximage = maximage;
indicatordata.meanimage = meanimage;
if fid > 0
fclose(fid);
end 

%%
% imagediff = max(meanimage)-min(meanimage);
% if imagediff <= 10
%     error('Change is indicator meanimage value is small: %d, no indicator light',imagediff);
% else
try
        
    IRthreshold = max(meanimage)-min(meanimage)/2;
    % pad need to be more than 1/2 duty cycle / sampling frequency
    pad = indicator_params.pad;
    indicatorLED = [];
%     thresholdimageorig = zeros(1,numel(meanimage));
%     idx = [];
%     idx = find(meanimage > IRthreshold );
%     thresholdimageorig(idx) = 1;
    thresholdimageorig = meanimage > IRthreshold;
    indicatorLED.StartEndStatus = logical([thresholdimageorig(1),thresholdimageorig(end)]);
    indicatorLED.indicatordigital = thresholdimageorig;
    % pad threshold image
    thresholdimage = zeros(1,numel(meanimage)+2*pad);
    thresholdimage(pad+1:end-pad) = thresholdimageorig;
    %find the start of light pulses
    onpulsecount = 0;
    for i = pad+1:(numel(thresholdimage)-pad)
        if (thresholdimage(i) == 1) && (mean(thresholdimage(i-pad:i-1)) == 0);
            onpulsecount = onpulsecount+1;
            indicatorLED.startframe(onpulsecount) = i-pad;
        end
        
    end
    %find the end of light pulses
    offpulsecount = 0;
    for i = pad+1:numel(thresholdimage)-pad
        if (thresholdimage(i) == 0) && thresholdimage(i-1) == 1 && (mean(thresholdimage(i+1:i+pad)) == 0);
            offpulsecount = offpulsecount+1;
            indicatorLED.endframe(offpulsecount) = i-pad;
        end
    end
    indicatorLED.starttimes = headerinfo.timestamps(indicatorLED.startframe)';
    indicatorLED.endtimes = headerinfo.timestamps(indicatorLED.endframe)';
    
    indicatordata.indicatorLED = indicatorLED;
    catch
        imagediff = max(meanimage)-min(meanimage);
        fprintf(logfid,'Error calculating indicator on-off: meanimage value is small: %d, no indicator light',imagediff)
        success = false;
        msgs = {'Error calculating indicator LED on-off: no indicator detected'};
        return;
end

%% plot figure
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
  fprintf(logfid,'Saved indicator data to file %s\n',indicatordatamatfile);
else
  fprintf(logfid,'Could not save indicator data to mat file:\n%s\n',getReport(ME));
end
%% close log

fprintf(logfid,'Finished running FlyBowlRegisterTrx at %s.\n',datestr(now,'yyyymmddTHHMMSS'));

if logfid > 1,
  fclose(logfid);
end

%% close figures

if isdeployed,
  delete(findall(0,'type','figure'));
end

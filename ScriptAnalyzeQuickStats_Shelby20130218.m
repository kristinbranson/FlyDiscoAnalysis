%% set up path

addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
addpath /groups/branson/bransonlab/projects/olympiad/FBDC;
rootdatadir = '/groups/branson/bransonlab/projects/CourtshipBowls/data/Shelby';
outdir = '/groups/branson/bransonlab/projects/CourtshipBowls/DataCapture';
params_file = '../FlyBowlDataCapture/FlyBowlDataCaptureParams_EP00014.txt';

%% parameters

tmp = dir(fullfile(rootdatadir,'shelbymbscreen*'));
expdirs = cellfun(@(x) fullfile(rootdatadir,x),{tmp([tmp.isdir]).name},'UniformOutput',false);

% name of ufmf diagnostics file
UFMFDiagnosticsFileStr = 'ufmf_diagnostics.txt';
% minimum fraction of experiment directories required to store stream
%ufmfstream_minfracexps = .9;
% mat file to save stats to 
savename = fullfile(outdir,['QuickStats_Shelby_Stats_',datestr(now,30),'.mat']);

%% read the parameter file

GUIi = 1;

% comment character in params file
comment_char = '#';

% open the parameter file
fid = fopen(params_file,'r');
if fid < 0,
  s = sprintf('Could not read in parameters file %s',params_file);
  uiwait(errordlg(s,'Error reading parameters'));
  error(s); %#ok<SPERR>
end

params = struct;
try
  % read each line
  while true,
    s = fgetl(fid);
    if ~ischar(s), break; end
    
    % remove extra white space
    s = strtrim(s);
    
    % skip comments
    if isempty(s) || s(1) == comment_char,
      continue;
    end
    
    % split at ,
    v = regexp(s,',','split');
    
    % first value is the parameter name, rest are parameter values
    params.(v{1}) = v(2:end);
  end
  fclose(fid);
  
  % some values are numbers
  numeric_params = {'PreAssayHandling_CrossDate_Range',...
    'PreAssayHandling_SortingDate_Range','PreAssayHandling_StarvationDate_Range',...
    'PreAssayHandling_SortingHour_Range','PreAssayHandling_StarvationHour_Range',...
    'PreAssayHandling_SortingHour_Interval','PreAssayHandling_StarvationHour_Interval',...
    'Imaq_ROIPosition','NFlies','RecordTime','PreviewUpdatePeriod',...
    'MetaData_RoomTemperatureSetPoint','MetaData_RoomHumiditySetPoint',...
    'FrameRatePlotYLim','TempPlotYLim','DoQuerySage','Imaq_FrameRate',...
    'Imaq_Shutter','Imaq_Gain','Imaq_Brightness','TempProbePeriod','TempProbeChannels',...
    'TempProbeReject60Hz','DoRecordTemp','NPreconSamples','gdcamPreviewFrameInterval',...
    'Imaq_MaxFrameRate','UFMFPrintStats','UFMFStatStreamPrintFreq',...
    'UFMFStatComputeFrameErrorFreq','UFMFStatPrintTimings',...
    'UFMFMaxFracFgCompress','UFMFMaxBGNFrames','UFMFBGUpdatePeriod',...
    'UFMFBGKeyFramePeriod','UFMFMaxBoxLength','UFMFBackSubThresh',...
    'UFMFNFramesInit','UFMFBGKeyFramePeriodInit','ColormapPreview',...
    'ScanLineYLim','MinFliesLoadedTime','MaxFliesLoadedTime',...
    'PreAssayHandling_FlipUsed','WishListRange',...
    'DoSyncBarcode','flip_days','CheckBarcode','CoupleCameraTempProbeStart'};
  for i = 1:length(numeric_params),
    if isfield(params,numeric_params{i}),
      params.(numeric_params{i}) = str2double(params.(numeric_params{i}));
    else
      fprintf('Parameter %s not set in parameter file.\n',numeric_params{i});
    end
  end
  
  % some values are not lists
  notlist_params = {'Imaq_Adaptor','Imaq_DeviceName','Imaq_VideoFormat',...
    'FileType','MetaData_AssayName',...
    'MetaData_Effector','MetaData_Gender','MetaDataFileName','MovieFilePrefix','LogFileName',...
    'UFMFLogFileName','UFMFStatFileName','PreconSensorSerialPort',...
    'DoRotatePreviewImage',...
    'QuickStatsStatsFileName',...
    'ScreenType','ScreenReason',...
    'Assay_Room'};
  for i = 1:length(notlist_params),
    fn = notlist_params{i};
    if ischar(params.(fn){1}),
      tmp = params.(fn){1};
      for j = 2:length(params.(fn)),
        tmp = [tmp,',',params.(fn){j}]; %#ok<AGROW>
      end
      params.(fn) = tmp;
    else
      params.(fn) = cat(2,params.(fn){:});
    end
  end
  
  % parameters that are selected by GUI instance
  GUIInstance_params = {'OutputDirectory','TmpOutputDirectory','HardDriveName'};
  for i = 1:length(GUIInstance_params),
    fn = GUIInstance_params{i};
    j = mod(GUIi-1,length(params.(fn)))+1;
    if iscell(params.(fn)),
      params.(fn) = params.(fn){j};
    else
      params.(fn) = params.(fn)(j);
    end
  end

catch ME
  uiwait(errordlg({'Error parsing parameter file:',getReport(ME)},'Error reading parameters'));
  rethrow(ME);
end

params.MovieFileStr = sprintf('%s.%s',params.MovieFilePrefix,params.FileType);


%% create quickstats params

ComputeQuickStatsParams = {...
  'UFMFDiagnosticsFileStr',params.UFMFStatFileName,...
  'MovieFileStr',params.MovieFileStr,...
  'MetaDataFileStr',params.MetaDataFileName,...
  'FigHandle',100,...
  'GUIInstance',GUIi,...
  'SaveFileStr','QuickStats.png',...
  'SaveDataStr','QuickStats.txt'...
  'ScanLineYLim',params.ScanLineYLim...
  };

%% compute quick stats for these experiments

figpos = [100,100,1000,800];

for i = 1:numel(expdirs),
  expdir = expdirs{i};
  [QuickStats,success,errmsg,warnings] = computeQuickStats(expdir,...
    ComputeQuickStatsParams{:},...
    'figpos',figpos);
  if success,
    close(QuickStats.fig);
    close(QuickStats.showufmf_handle);
  else
    fprintf('Error: %s\n',errmsg);
  end
end

%% do it
quickstats_stats = AnalyzeQuickStats(...
  'UFMFDiagnosticsFileStr',UFMFDiagnosticsFileStr,...
  'savename',savename,...
  'rootdir',rootdatadir,...
  'expdirs',expdirs);
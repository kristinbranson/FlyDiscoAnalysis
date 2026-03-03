% TestCtraxStatsBase

%% set up path

if ispc,
  addpath ../JCtrax/filehandling/
  addpath ../JCtrax/misc/
else
  addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling/
  addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc/
end
addpath ../FlyBowlDataCapture

%% data

% scratched polycarbonate
if ispc,
  rootdatadir = 'E:\Data\FlyBowl\polycarbonate_scratched';
else
  rootdatadir = '/groups/branson/bransonlab/tracking_data/olympiad/FlyBowl/polycarbonate_scratched';
end
expdir = 'DL-wildtype_TrpA_Rig1Plate01BowlC_20101112T152624';

%% parameters

detectregistrationparams = struct;
detectregistrationparams.method = 'circle_manual';
detectregistrationparams.circleImageType = 'canny';
detectregistrationparams.circleRLim = [.4,.55];
detectregistrationparams.circleXLim = [.3,.6];
detectregistrationparams.circleYLim = [.3,.6];
detectregistrationparams.circleImageThresh = [];
detectregistrationparams.circleCannyThresh = [];
detectregistrationparams.circleCannySigma =[];
detectregistrationparams.circleNXTry = 50;
detectregistrationparams.circleNYTry = 50;
detectregistrationparams.circleNRTry = 50;
detectregistrationparams.maxDistCornerFrac_BowlLabel = .17;
detectregistrationparams.featureRadius = 25;
detectregistrationparams.circleRadius_mm = 63.5;
detectregistrationparams.bowlMarkerPairTheta_true = -3*pi/4;
detectregistrationparams.debug = true;
detectregistrationparams.isBowlMarker = false;

%% initialize data

obj = CtraxStatsBase('rootdatadir',rootdatadir,'ctraxfilestr','fixerrors_results.mat',...
  'NOWRITEACCESS',false,'detectregistrationparams',detectregistrationparams);
obj.AddExpDir(expdir);

%% derive some more measurements
obj.ComputeSpeedMeasurements();
obj.ComputeLandmarkMeasurements();

%% histogramming parameters
obj.histogrammeasurements_ploterrorbars = 'std';
obj.histogrammeasurements_plotperfly = true;
obj.histogrammeasurements_plotperexp = true;
obj.histogrammeasurements_nbins = 25;

warning('off','MATLAB:Axes:NegativeDataInLogAxis');

%% histogramming two measurements parameters

%% 

obj.HistogramTwoMeasurements('x_mm','y_mm');

%% histogram stuff

if true,

% fields to histogram
fns_hist = {'velmag','du_cor','absdv_cor','absdtheta','absyaw','dist2wall','ddist2wall','angle2wall','dangle2wall'};

% special bin mode
%binmode_log = {'velmag','absdv_cor','absdtheta','absyaw','dist2wall','angle2wall'};
binmode_log = {};
binmode_logabs = {};
               
%binmode_logabs = {'du_cor','ddist2wall','dangle2wall'};

% special output functions
outputfun = struct;
outputfun.angle2wall = @abs;

% left-most edge should be 0
zerolim = [binmode_log,{'angle2wall'}];

% special lim_prctile
lim_prctile = struct;
lim_prctile.absdv_cor = [0,95];
lim_prctile.absdtheta = [0,95];

% linear scale for y
ylinearscale = {};

hfig = 10;
for i = 1:length(fns_hist),

  fn = fns_hist{i};
  
  fprintf('\n** %s **\n',fn);
  
  % specify parameters per fn
  if ismember(fn,binmode_log),
    binmode_curr = 'log';
  elseif ismember(fn,binmode_logabs),
    binmode_curr = 'logabs';
  else
    binmode_curr = 'linear';
  end
  fprintf('  bin mode: %s\n',binmode_curr);
  if isfield(outputfun,fn),
    outputfun_curr = outputfun.(fn);
    fprintf('  output function: %s\n',func2str(outputfun_curr));
  else
    outputfun_curr = [];
  end
  lim_curr = nan(1,2);
  if ismember(fn,zerolim),
    fprintf('  forcing binning to begin at 0\n');
    lim_curr(1) = 0;
  end
  lim_prctile_curr = obj.histogrammeasurements_lim_prctile;
  if isfield(lim_prctile,fn),
    lim_prctile_curr = lim_prctile.(fn);
  end
  fprintf('  percentiles of data for choosing binning limits: [%.1f,%.1f]\n',lim_prctile_curr);

  ylogscale_curr = ~ismember(fn,ylinearscale);
  if ylogscale_curr,
    fprintf('  setting y-axis scale to log\n');
  else
    fprintf('  setting y-axis scale to linear\n');
  end
  obj.HistogramMeasurement(fn,'hfig',hfig,'binmode',binmode_curr,...
    'outputfun',outputfun_curr,'lim',lim_curr,'lim_prctile',lim_prctile_curr,...
    'ylogscale',ylogscale_curr,'ploterrorbars','std');

  hfig = hfig + 1;
  
  pause(5);
end

end
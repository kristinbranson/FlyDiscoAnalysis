%% set up path
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling/
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc/
addpath ../FlyBowlDataCapture

%% data

% scratched polycarbonate
rootdatadir = '/groups/branson/bransonlab/tracking_data/olympiad/FlyBowl/polycarbonate_scratched';
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

obj = CtraxDiagnostics('rootdatadir',rootdatadir,...
  'ctraxfilestr','fixerrors_results.mat',...
  'NOWRITEACCESS',false,...
  'detectregistrationparams',detectregistrationparams);
obj.AddExpDir(expdir);

%% derive some more measurements
obj.ComputeSpeedMeasurements();
obj.ComputeLandmarkMeasurements();

%% where in the arena do the flies spend time?

hfig = 1;
obj.CenterPositionHeatmap('hfig',hfig,'doplot',true,'nbins',50,...
  'minspeed',1.3,'fraclogscale',true);


%obj.CenterPositionHeatmapPolar('hfig',hfig,'doplot',true,'nbins_r',5,'nbins_theta',20);
%hfig = hfig+1;
%obj.CenterPositionBinEntriesPolar('hfig',hfig,'doplot',true,'nbins_r',5,'nbins_theta',20);

%% is there a bias toward one side of the arena

hfig = 5;
obj.BowlAngleBias('doplot',true,'hfig',hfig);
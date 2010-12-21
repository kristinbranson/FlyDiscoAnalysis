%% set up path
if ispc,
  addpath ../JCtrax/filehandling/
  addpath ../JCtrax/misc/
  addpath ../FlyBowlDataCapture
  addpath(genpath('../FlyBowlDataCapture/jfrc_metadata_tools/src'));
else
  addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling/
  addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc/
  addpath ../FlyBowlDataCapture
  addpath(genpath('../FlyBowlDataCapture/jfrc_metadata_tools/src'));
end

%% data

%datatype = 'scratched_polycarbonate_CtraxTest20101118';
datatype = 'pBDPGAL4U_CtraxTest20101118';
%datatype = 'GAL4_CtraxTest20101118';

switch datatype,

  case 'scratched_polycarbonate_CtraxTest20101118',
  
    % scratched polycarbonate
    rootdatadir = '/groups/branson/bransonlab/tracking_data/olympiad/FlyBowl/polycarbonate_scratched';
    expdirs = {'DL-wildtype_TrpA_Rig1Plate01BowlC_20101112T152624'};
    resultsdir = '/groups/branson/home/bransonk/tracking/data/olympiad/FlyBowl/CtraxTest20101118';
    ctraxfilestr = 'fixerrors_results.mat';
    NOWRITEACCESS = false;
    ExpBGFGModelMatFile = '';
    
  case 'pBDPGAL4U_CtraxTest20101118',
    
    % controls from 10/19 - 10/21
    rootdatadir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data/';
    resultsdir = '/groups/branson/home/bransonk/tracking/data/olympiad/FlyBowl/CtraxTest20101118';
    expdirs = {...
      %'pBDPGAL4U_TrpA_Rig1Plate01BowlA_20101019T114325',...
      'pBDPGAL4U_TrpA_Rig1Plate01BowlA_20101020T141505',...
      'pBDPGAL4U_TrpA_Rig1Plate01BowlA_20101021T102920',...
      'pBDPGAL4U_TrpA_Rig1Plate01BowlA_20101021T133637',...
      'pBDPGAL4U_TrpA_Rig1Plate01BowlB_20101019T114334',...
      'pBDPGAL4U_TrpA_Rig1Plate01BowlB_20101020T101842',...
      'pBDPGAL4U_TrpA_Rig1Plate01BowlB_20101020T141510',...
      'pBDPGAL4U_TrpA_Rig1Plate01BowlC_20101020T141438',...
      'pBDPGAL4U_TrpA_Rig1Plate01BowlC_20101021T102851',...
      'pBDPGAL4U_TrpA_Rig1Plate01BowlC_20101021T133605',...
      'pBDPGAL4U_TrpA_Rig1Plate01BowlD_20101020T101817',...
      'pBDPGAL4U_TrpA_Rig1Plate01BowlD_20101020T141444',...
      'pBDPGAL4U_TrpA_Rig1Plate01BowlD_20101021T133611'};
    ctraxfilestr = 'fixerrors_results.mat';
    NOWRITEACCESS = true;    
    ExpBGFGModelMatFile = '/groups/branson/bransonlab/tracking_data/olympiad/FlyBowl/CtraxTest20101013/LearnCtraxParams/ExpBGFGModelResults20101029.mat';
    
  case 'GAL4_CtraxTest20101118',
    
    
    % GAL4 from 10/19 - 10/21
    rootdatadir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data/';
    resultsdir = '/groups/branson/home/bransonk/tracking/data/olympiad/FlyBowl/CtraxTest20101118';
    expdirs = {...
      'GMR_13F01_AE_01_TrpA_Rig1Plate01BowlC_20101020T092643'
      %'GMR_13F07_AE_01_TrpA_Rig1Plate01BowlB_20101019T111623'
      'GMR_14D12_AE_01_TrpA_Rig1Plate01BowlA_20101020T110814'
      'GMR_15D05_AE_01_TrpA_Rig1Plate01BowlB_20101021T142528'
      'GMR_15E07_AE_01_TrpA_Rig1Plate01BowlB_20101021T154143'
      %'GMR_15H01_AE_01_TrpA_Rig1Plate01BowlA_20101019T101741'
      'GMR_16B02_AE_01_TrpA_Rig1Plate01BowlD_20101021T131205'
      'GMR_16B04_AE_01_TrpA_Rig1Plate01BowlD_20101020T131236'
      'GMR_16B10_AE_01_TrpA_Rig1Plate01BowlC_20101020T144254'
      'GMR_16C09_AE_01_TrpA_Rig1Plate01BowlA_20101020T153925'
      %'GMR_16E02_AE_01_TrpA_Rig1Plate01BowlB_20101019T150721'
      %'GMR_16E08_AE_01_TrpA_Rig1Plate01BowlD_20101019T141236'
      %'GMR_16G08_AE_01_TrpA_Rig1Plate01BowlC_20101019T131208'
      'GMR_19F01_AE_01_TrpA_Rig1Plate01BowlA_20101021T091513'
      'GMR_22D03_AE_01_TrpA_Rig1Plate01BowlD_20101021T111833'
      %'GMR_42A08_AE_01_TrpA_Rig1Plate01BowlC_20101019T092128'
      };
    ctraxfilestr = 'ctrax_results.mat';
    NOWRITEACCESS = true;    
    ExpBGFGModelMatFile = '/groups/branson/bransonlab/tracking_data/olympiad/FlyBowl/CtraxTest20101013/LearnCtraxParams/ExpBGFGModelResults20101029.mat';
    
end

%% parameters

if ~isempty(regexpi(datatype,'polycarbonate')),
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
  detectregistrationparams.debug = true;
  detectregistrationparams.isBowlMarker = false;
  bowl2MarkerAngle = {};
else
  detectregistrationparams = struct(...
    'bkgdNSampleFrames',10,...
    'method','normcorr',...
    'crossFilterRadius',9,...
    'nRotations',20,...
    'minCrossFeatureStrength',.92,...
    'minFeatureStrengthHigh',30,...
    'minDistCenterFrac',.5,...
    'maxDistCenterFrac',.57,...
    'maxDistCornerFrac_BowlLabel',.2,...
    'nRegistrationPoints',8,...
    'featureRadius',25,...
    'maxDThetaMate',10*pi/180,...
    'pairDist_mm',133,...
    'maxDThetaBowlMarkerPair',pi/6,...
    'markerPairAngle_true',pi/6,...
    'debug',true,...
    'nr',1024,'nc',1024);
  bowl2MarkerAngle = {...
    'A',3*pi/4
    'B',pi/4
    'C',3*pi/4
    'D',pi/4};
  
  ExpBGFGModelMatFile = '/groups/branson/bransonlab/tracking_data/olympiad/FlyBowl/CtraxTest20101013/LearnCtraxParams/ExpBGFGModelResults20101029.mat';

end

%% initialize data

obj = CtraxDiagnostics('rootdatadir',rootdatadir,...
  'ctraxfilestr',ctraxfilestr,...
  'NOWRITEACCESS',NOWRITEACCESS,...
  'detectregistrationparams',detectregistrationparams,...
  'bowl2MarkerAngle',bowl2MarkerAngle,...
  'ExpBGFGModelMatFile',ExpBGFGModelMatFile);

%% add all the experiment directories

%for i = 1:numel(expdirs),
i = 1;
  obj.AddExpDir(expdirs{i});
%  drawnow;
%  tmp = fullfile(resultsdir,'figs',sprintf('registration_%s',obj.expdir_bases{i}));
%  savefig(tmp,i,'png');
%end

%% plot some trajectories

fliesplot = nan(2,obj.nexpdirs);
for n = 1:obj.nexpdirs,
  [~,i] = max([obj.trx(obj.movie2flies{n}).nframes]);
  fliesplot(:,n) = [obj.movie2flies{n}(1),obj.movie2flies{n}(end-1)];
end
fliesplot = reshape(fliesplot,[5,4]);
fliesplot = fliesplot(:)';
figure(1);
clf;
hax = createsubplots(4,5,[[.05,.01];[.05,.05]]);
hax = reshape(hax,[4,5])';
hax = hax(:)';
hax = obj.PlotTrajectories('flies',fliesplot,'hax',hax);
fliesplot = fliesplot';
for i = 1:numel(fliesplot),
  fly = fliesplot(i);
  n = obj.fly2movie(fly);
  title(hax(i),strrep(obj.expdir_bases{n},'_AE_01_TrpA_Rig1Plate01',''));
end

%% derive some more measurements
obj.ComputeSpeedMeasurements();
obj.ComputeLandmarkMeasurements();

%% plot histograms per experiment
minspeed = 1.5;
if ~isempty(regexpi(datatype,'GAL4')),
  for n = 1:obj.nexpdirs,
    expdirs = obj.expdir_bases(n);
    hfig = [n,100+n];
    obj.CenterPositionHeatmap('hfig',hfig,'doplot',true,'nbins',50,...
      'minspeed',minspeed,'fraclogscale',true,'plotperexp',true,'plotperfly',false,...
      'ploterrorbars','stderr','expdirs',expdirs);
    tmp = fullfile(resultsdir,'figs',sprintf('GAL4_CenterPositionHeatmap_%s',obj.expdir_bases{n}));
    savefig(tmp,hfig(1),'png');
    hfig = [200+n,300+n];
    obj.CenterPositionHeatmapPolar('hfig',hfig,'doplot',true,'minspeed',minspeed,...
      'fraclogscale',true,'plotperexp',false,'plotperfly',false,'ploterrorbars','stderr',...
      'nbins_r',10,'nbins_theta',20,'expdirs',expdirs);
    tmp = fullfile(resultsdir,'figs',sprintf('GAL4_CenterPositionHeatmapPolar_%s',obj.expdir_bases{n}));
    savefig(tmp,hfig(1),'png');
  end
end

%% where in the arena do the flies spend time?

minspeed = 1.5;
hfig = [1,2];
obj.CenterPositionHeatmap('hfig',hfig,'doplot',true,'nbins',50,...
  'minspeed',minspeed,'fraclogscale',true,'plotperexp',true,'plotperfly',false,...
  'ploterrorbars','stderr');
hfig = [3,4];
obj.CenterPositionHeatmapPolar('hfig',hfig,'doplot',true,'minspeed',minspeed,...
  'fraclogscale',true,'plotperexp',true,'plotperfly',false,'ploterrorbars','stderr',...
  'nbins_r',10,'nbins_theta',20);

%obj.CenterPositionHeatmapPolar('hfig',hfig,'doplot',true,'nbins_r',5,'nbins_theta',20);
%hfig = hfig+1;
%obj.CenterPositionBinEntriesPolar('hfig',hfig,'doplot',true,'nbins_r',5,'nbins_theta',20);

%% per bowl

Bowls = {'A','B','C','D'};
for i = 1:4,
  [ns,expdirscurr] = obj.Metadata2Exp('bowl',Bowls{i});
  obj.CenterPositionHeatmapPolar('doplot',true,'minspeed',minspeed,...
    'fraclogscale',true,'plotperexp',false,'plotperfly',false,'ploterrorbars','none',...
    'nbins_r',10,'nbins_theta',20,'expdirs',expdirscurr,'clim',[0,.03]);
end

%% is there a bias toward one side of the arena

hfig = 5;
Bowls = {'A','B','C','D'};
for i = 1:4,
  hfig = 4+i;
  [ns,expdirscurr] = obj.Metadata2Exp('bowl',Bowls{i});
  obj.BowlAngleBias('doplot',true,'hfig',hfig,'plotperexp',false,'ploterrorbars','stderr','expdirs',expdirscurr,'minspeed',minspeed);
  pause(1);
  tmp = fullfile(resultsdir,'figs',['BowlAngleBias_PerBowl_',Bowls{i}]);
  savefig(tmp,hfig,'png')
end

%% plot time series

fly = 1;
fnsplot = {'velmag','du_cor','absdv_cor','absdtheta','absyaw','dist2wall','ddist2wall','angle2wall','dangle2wall','dnose2ell','anglesub','veltoward_nose2ell','absthetadiff_center'};
obj.PlotTimeSeries(fly,fnsplot,'hfig',1,'tstart',400,'tend',500);



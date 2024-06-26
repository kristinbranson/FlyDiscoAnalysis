%% setup path
modpath

%% set parameters
% settingsdir = '/groups/branson/home/robiea/Code_versioned/FlyDiscoAnalysis/settings';
% analysis_protocol = '20210531_flybubble_LED_ARtestingperframe';
% settingsdir = '/groups/branson/home/robiea/Code_versioned/FlyDiscoAnalysis/settings-internal';
% analysis_protocol = '20230213_flybubble_LED_VNC2';


settingsdir = '/groups/branson/home/robiea/Code_versioned/FlyDiscoAnalysis';
analysis_protocol = 'settings-internal/20240507_flybubble_LED_VNC3';
params = {'settingsdir',settingsdir,...
  'analysis_protocol',analysis_protocol};

%%
% expdirs = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210713_perframefeatures/VNC_YNA_K_162984_RigB_20210602T135902'};

% expdirs = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_JRC_SS42702_RigC_20210505T153129', ...
% '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_JRC_SS43522_RigC_20210412T131721', ...
% '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_JRC_SS43660_RigA_20210414T135538', ...
% '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_JRC_SS44231_RigA_20210419T141517', ...
% '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_JRC_SS46301_RigC_20210512T131504', ...
% '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_JRC_SS48619_RigA_20210426T141945', ...
% '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_JRC_SS49172_RigA_20210504T145545', ...
% '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_JRC_SS49220_RigB_20210421T143507', ...
% '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_JRC_SS51893_RigA_20210506T150106', ...
% '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_JRC_SS52094_RigA_20210602T144112', ...
% '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_JRC_SS52799_RigD_20210602T145135', ...
% '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_JRC_SS52985_RigA_20210503T152107', ...
% '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_JRC_SS60240_RigB_20210420T141142', ...
% '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_JRC_SS62706_RigA_20210525T134347', ...
% '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_JRC_SS64175_RigB_20210518T145000', ...
% '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_JRC_SS64214_RigA_20210520T151504', ...
% '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_JRC_SS64218_RigA_20210527T162508', ...
% '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_JRC_SS64234_RigC_20210513T154200', ...
% '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_JRC_SS65692_RigC_20210427T135618', ...
% '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_JRC_SS65710_RigB_20210513T142817', ...
% '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_JRC_SS65765_RigA_20210525T140224', ...
% '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_JRC_SS66932_RigA_20210525T153325', ...
% '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_JRC_SS67323_RigA_20210513T145544', ...
% '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_YNA_K_162984_RigA_20210519T153103', ...
% '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_YNA_K_162984_RigB_20210602T135902', ...
% '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_YNA_K_162984_RigC_20210604T152231', ...
% '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_YNA_K_162984_RigD_20210526T155055'};

% expdirs = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_JRC_SS64175_RigB_20210518T145000'}

% expdirs = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20230509_autochecksERROR/VNC2_YNA_K_162984_RigD_20230503T135414'};

% expdir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20240508_VNC3/VNC3_YNA_K_162984_RigA_20240507T182023';
expdir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20240530_testACC/VNC3_JRC_SS57980_RigD_20240529T114756';

%% 
FlyDiscoAutomaticChecksComplete(expdir,params{:});
% success = [];
% msg = {};
% 
% for i = 1:numel(expdirs),
%   expdir = expdirs{i};
%   disp(expdir);
%   try 
%       FlyDiscoAutomaticChecksComplete(expdir,params{:})
%       success(i) = false ;  %#ok<SAGROW>
%       msg{i} = '' ;  %#ok<SAGROW>
%   catch me
%       success(i) = false ;  %#ok<SAGROW>
%       msg{i} = me.getReport() ;  %#ok<SAGROW>
%   end
% end

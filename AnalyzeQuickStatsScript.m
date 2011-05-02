%% set up path
if ispc,
  addpath ../JCtrax/filehandling/
  addpath ../JCtrax/misc/
  addpath ../FlyBowlDataCapture
  addpath E:\Code\SAGE\MATLABInterface\Trunk;
  rootdir = 'O:\Olympiad_Screen\fly_bowl\bowl_data';
  %addpath(genpath('../FlyBowlDataCapture/jfrc_metadata_tools/src'));
else
  addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling/
  addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc/
  addpath ../FlyBowlDataCapture
  addpath /groups/branson/bransonlab/projects/olympiad/SAGE/MATLABInterface/Trunk;
  rootdir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data';
  %addpath(genpath('../FlyBowlDataCapture/jfrc_metadata_tools/src'));
end

%% parameters

experiment_params = struct;
% type of data to analyze
experiment_params.experiment_protocol = {'*8*','*7*'};
experiment_params.rearing_protocol = '*8*';
% what dates should we analyze
experiment_params.daterange = {};
% name of ufmf diagnostics file
UFMFDiagnosticsFileStr = 'ufmf_diagnostics.txt';
% minimum fraction of experiment directories required to store stream
%ufmfstream_minfracexps = .9;
% mat file to save stats to 
savename = ['QuickStats_Stats_',datestr(now,30),'.mat'];

%% do it
quickstats_stats = AnalyzeQuickStats(...
  'experiment_params',struct2paramscell(experiment_params),...
  'UFMFDiagnosticsFileStr',UFMFDiagnosticsFileStr,...
  'savename',savename,...
  'rootdir',rootdir);
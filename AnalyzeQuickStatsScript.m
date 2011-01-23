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

%% parameters

experiment_params = struct;
% type of data to analyze
experiment_params.protocol = 'CtraxTest20110111';
% what dates should we analyze
experiment_params.daterange = cell(1,2);
% what lines
experiment_params.linename = 'pB*';
% whether the experiments did not start
experiment_params.notstarted = false;
% name of ufmf diagnostics file
UFMFDiagnosticsFileStr = 'ufmf_diagnostics.txt';
% minimum fraction of experiment directories required to store stream
ufmfstream_minfracexps = .9;
% mat file to save stats to 
savename = ['QuickStats_Stats_',datestr(now,30),'.mat'];

%% do it
quickstats_stats = AnalyzeQuickStats(...
  'experiment_params',struct2paramscell(experiment_params),...
  'UFMFDiagnosticsFileStr',UFMFDiagnosticsFileStr,...
  'ufmfstream_minfracexps',ufmfstream_minfracexps,...
  'savename',savename);
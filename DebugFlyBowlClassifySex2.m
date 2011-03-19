%% set up path

if ispc,
  addpath E:\Code\JCtrax\misc;
  addpath E:\Code\JCtrax\filehandling;
  addpath E:\Code\hmm;
  settingsdir = 'E:\Code\FlyBowlAnalysis\settings';
  rootdir = 'E:\Data\FlyBowl\bowl_data';
  expdir = 'E:\Data\FlyBowl\bowl_data\GMR_14B02_AE_01_TrpA_Rig2Plate14BowlD_20110204T094327';
else
  addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
  addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
  addpath /groups/branson/home/bransonk/tracking/code/lds/hmm;
  settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
  rootdir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data';
  expdir = fullfile(rootdir,'GMR_12E07_AE_01_TrpA_Rig2Plate14BowlA_20110209T134320');
end

%% parameters

analysis_protocol = '20110222';
params = {'settingsdir',settingsdir,...
  'analysis_protocol',analysis_protocol};

%% 

FlyBowlClassifySex2(expdir,params{:});

% %% 
% 
% nexpdirs_gmr_analyze = 10;
% nexpdirs_control_analyze = 8;
% rootwritedir = '/groups/branson/bransonlab/tracking_data/olympiad/FlyBowl/CtraxTest20110315';
% 
% %%
% 
% [~,expdirs_gmr] = getExperimentDirs('daterange',{'20110311T000000',''},...
%   'protocol',analysis_protocol,'linename','GMR*');
% [~,expdirs_control] = getExperimentDirs('daterange',{'20110311T000000',''},...
%   'protocol',analysis_protocol,'linename','pBD*');
% if numel(expdirs_gmr) > nexpdirs_gmr_analyze,
%   expdirs_gmr = expdirs(randperm(numel(expdirs_gmr)));
%   expdirs_gmr = expdirs_gmr(1:nexpdirs_gmr_analyze);
% end
% if numel(expdirs_control) > nexpdirs_control_analyze,
%   expdirs_control = expdirs_gmr(randperm(numel(expdirs_control)));
%   expdirs_control = expdirs_control(1:nexpdirs_control_analyze);
% end
% expdirs = [expdirs_gmr,expdirs_control];
% 
% %%
% 
% if false,
% 
% expdirs_write = cell(size(expdirs));
% for i = 1:numel(expdirs),
%   [~,expdir_base] = fileparts(expdirs{i});
%   expdirs_write{i} = fullfile(rootwritedir,expdir_base);
%   fprintf('Copying %s to %s\n',expdirs{i},expdirs_write{i});
%   if ~exist(expdirs_write{i},'file'),
%     copyfile(expdirs{i},expdirs_write{i});
%   end
% end
% expdirs = expdirs_write;
% else
  expdirs = dir(fullfile(rootwritedir,'*_*'));
  expdirs = cellfun(@(s) fullfile(rootwritedir,s),{expdirs([expdirs.isdir]).name},'UniformOutput',false);
% end
% 
% %%
% 
% for i = 1:numel(expdirs),
%   FlyBowlClassifySex2(expdirs{i},params{:},'dosave',true);
% end

%%

for i = 1:numel(expdirs),
  FlyBowlComputePerFrameFeatures(expdirs{i},params{:});
  FlyBowlComputePerFrameStats(expdirs{i},params{:});
  FlyBowlPlotPerFrameStats(expdirs{i},params{:});
end


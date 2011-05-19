%% set up path

if ispc,
  addpath E:\Code\JCtrax\misc;
  addpath E:\Code\JCtrax\filehandling;
  addpath('E:\Code\SAGE\MATLABInterface\Trunk\')
  settingsdir = 'E:\Code\FlyBowlAnalysis\settings';
  rootdatadir = 'O:\Olympiad_Screen\fly_bowl\bowl_data';
  rootoutdir = 'L:\projects\olympiad\data\flybowl\albins';
  
else
  addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
  addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
  addpath /groups/branson/bransonlab/projects/olympiad/SAGE/MATLABInterface/Trunk;
  settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
  rootdatadir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data';
  rootoutdir = '/groups/branson/bransonlab/projects/olympiad/data/flybowl/albins';
end

%% parameters

analysis_protocol = '20110407';
starved_handling_protocol = 'HP_flybowl_v005p0.xls';
fed_handling_protocol = 'HP_flybowl_v005p1.xls';
Trp_effector = 'UAS_dTrpA1_2_0002';
CS_effector = 'CTRL_CantonS_1101243_0016';
compare_name = 'fed_50H05_compare';

compare_queries = struct;
compare_queries.fed_50H05_Trp = {'line_name','GMR_50H05_AE_01','effector',Trp_effector,'handling_protocol',fed_handling_protocol};
compare_queries.fed_50H05_CS = {'line_name','GMR_50H05_AE_01','effector',CS_effector,'handling_protocol',fed_handling_protocol};

%% output locations per condition

fns = fieldnames(compare_queries);
expdirs = cell(1,numel(fns));
for i = 1:numel(fns),
  fn = fns{i};
  expdirs{i} = fullfile(rootoutdir,fn);
end

%% combine experiments per condition

fns = fieldnames(compare_queries);
for i = 1:numel(fns),
  fn = fns{i};
  FlyBowlCombineExperiments(rootdatadir,expdirs{i},...
    'analysis_protocol',analysis_protocol,...
    'settingsdir',settingsdir,...
    compare_queries.(fn){:});
end

%% create comparison histograms

outdir = fullfile(rootoutdir,compare_name);

CompareHistograms2(expdirs,outdir,...
  'analysis_protocol',analysis_protocol,...
  'settingsdir',settingsdir,...
  'visible','on');

  
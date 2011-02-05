function FlyBowlComputePerFrameFeatures(expdir,varargin)

[analysis_protocol,settingsdir,datalocparamsfilestr] = ...
  myparse(varargin,...
  'analysis_protocol','current',...
  'settingsdir','/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings',...
  'datalocparamsfilestr','dataloc_params.txt');

%% read in the data locations
datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
dataloc_params = ReadParams(datalocparamsfile);
      
%% read in the per-frame parameters

perframeparamsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.perframeparamsfilestr);
perframe_params = ReadParams(perframeparamsfile);

%% load the trx

trx = Trx(
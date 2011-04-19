function [success,msgs] = FlyBowlAutomaticChecks(expdir,varargin)

success = true;
msgs = {};

[analysis_protocol,settingsdir,datalocparamsfilestr] = ...
  myparse(varargin,...
  'analysis_protocol','current',...
  'settingsdir','/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings',...
  'datalocparamsfilestr','dataloc_params.txt');

%% parameters

datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
dataloc_params = ReadParams(datalocparamsfile);

%% check for Ctrax bugs: nan, inf in trajectories

ctraxfile = fullfile(expdir,dataloc_params.ctraxfilestr);
if ~exist(ctraxfile,'file'),
  msgs{end+1} = sprintf('Ctrax output mat file %s does not exist',ctraxfile);
  success = false;
else
  % name of annotation file
  annfile = fullfile(expdir,dataloc_params.annfilestr);
  
  % name of movie file
  moviefile = fullfile(expdir,dataloc_params.moviefilestr);
  
  % load trajectories
  [trx,~,succeeded,timestamps] = load_tracks(ctraxfile,moviefile,'annname',annfile); %#ok<NASGU>
  if ~succeeded,
    msgs{end+1} = sprintf('Could not load trajectories from file %s',ctraxfile);
    success = false;
  else
    
    for fly = 1:numel(trx),
      badidx = isnan(trx(fly).x) | ...
        isnan(trx(fly).y) | ...
        isnan(trx(fly).a) | ...
        isnan(trx(fly).b) | ...
        isnan(trx(fly).theta);
      if any(badidx),
        [starts,ends] = get_interval_ends(badidx);
        starts = starts - trx(fly).off;
        ends = ends - trx(fly).off;
        msgs{end+1} = [sprintf('Trajectory %d has NaNs in frames',fly),sprintf(' %d-%d',[starts,ends]')];
        success = false;
      end
      badidx = isinf(trx(fly).x) | ...
        isinf(trx(fly).y) | ...
        isinf(trx(fly).a) | ...
        isinf(trx(fly).b) | ...
        isinf(trx(fly).theta);
      if any(badidx),
        [starts,ends] = get_interval_ends(badidx);
        starts = starts - trx(fly).off;
        ends = ends - trx(fly).off;
        msgs{end+1} = [sprintf('Trajectory %d has Infs in frames',fly),sprintf(' %d-%d',[starts,ends]')];
        success = false;
      end
      
    end
  end
end
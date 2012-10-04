function success = FlyBowlFixPerFrameFeatures(expdir,varargin)

success = true;

traj_fns = {'x_mm','y_mm','a_mm','b_mm','theta_mm','sex',...
  'x','y','theta','a','b','timestamps'};
fnsrecompute = {'du_tail','dv_tail',...
  'magveldiff_nose2ell','magveldiff_anglesub',...
  'veltoward_nose2ell','veltoward_anglesub',...
  'closestfly_nose2ell_angle_min30to30','dnose2ell_nose2ell_angle_min30to30',...
  'closestfly_nose2ell_angle_min20to20','dnose2ell_nose2ell_angle_min20to20',...
  'closestfly_nose2ell_angle_30tomin30','dnose2ell_nose2ell_angle_30tomin30'};

[analysis_protocol,settingsdir,datalocparamsfilestr,fnsrecompute,DEBUG] = ...
  myparse(varargin,...
  'analysis_protocol','current',...
  'settingsdir','/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings',...
  'datalocparamsfilestr','dataloc_params.txt',...
  'fnsrecompute',fnsrecompute,...
  'debug',false...
  );

[~,experiment_name] = fileparts(expdir);

%% list of all per-frame functions

datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
dataloc_params = ReadParams(datalocparamsfile);
% perframefnsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.perframefnsfilestr);
% perframefns = importdata(perframefnsfile);
% perframefns = setdiff(perframefns,[traj_fns,{'dt'}]);
perframedir = fullfile(expdir,dataloc_params.perframedir);
perframefns = dir(fullfile(perframedir,'*.mat'));
perframefns = regexprep({perframefns.name},'\.mat$','');
isscoreorlabel = ~cellfun(@isempty,regexp(perframefns,'([sS]core)|([lL]abel)','once'));
perframefns(isscoreorlabel) = [];
perframefns = setdiff(perframefns,{'dt'});

if ~exist(perframedir,'dir'),
  fprintf('No per-frame directory for experiment %s\n',experiment_name);
  return;
end

%% load trx

trxfile = fullfile(expdir,dataloc_params.trxfilestr);
if ~exist(trxfile,'file'),
  fprintf('Trajectory file %s does not exist\n',trxfile);
  success = false;
  return;
end
try
  td = load(trxfile);
  if ~all(isfield(td,{'trx','timestamps'})),
    error('trx and/or timestamps fields missing.');
  end
catch ME,
  fprintf('Error loading trx and timestamp data from file %s: %s\n',trxfile,getReport(ME));
  success = false;
  return;
end
trx = td.trx;
timestamps = td.timestamps;



%% fix trajectory files

for i = 1:numel(traj_fns),
  fn = traj_fns{i};
  filename = fullfile(perframedir,[fn,'.mat']);
  if ~exist(filename,'file'),
    fprintf('Trajectory perframe file %s does not exist\n',filename);
    continue;
  end
  try
    pd = load(filename,'data','units');
  catch ME,
    fprintf('Error loading data and units from perframefile %s: %s\n',filename,getReport(ME));
    continue;
  end
  if ~iscell(pd.data),
    fprintf('Fixing trajectory perframe feature %s, experiment %s.\n',fn,experiment_name);
    pd.data = {trx.(fn)};
    if ~DEBUG,
      delete(filename);
      save(filename,'-struct','pd');
    end
  end
end

%% change dt in all per-frame features that used it

tofix = fnsrecompute;

if isfield(td,'dtfixed') && td.dtfixed,
  fprintf('dt already fixed for trajectory file, skipping fixing dt in experiment %s\n',experiment_name);
else
  
  tmp = diff(timestamps);
  tmp = tmp(~isnan(tmp));
  if isempty(tmp),
    error('No non-nan timestamps in %s',trxfile);
  end
  meddt = median(tmp(~isnan(tmp)));
  dt = {trx.dt};
  
  for fni = 1:numel(perframefns),
    fn = perframefns{fni};
    if ismember(fn,fnsrecompute),
      continue;
    end
    filename = fullfile(perframedir,[fn,'.mat']);
    if ~exist(filename,'file'),
      fprintf('Per-frame file %s does not exist\n',filename);
      continue;
    end
    try
      pd = load(filename);
      if ~all(isfield(pd,{'data','units'})),
        error('data and/or units missing from file');
      end
    catch ME,
      fprintf('Error loading data and units from perframefile %s: %s\n',filename,getReport(ME));
      continue;
    end
    if isfield(pd,'dtfixed') && pd.dtfixed,
      fprintf('dt already fixed for feature %s, experiment %s, skipping\n',fn,experiment_name);
      continue;
    end
    nden = numel(strcmp(pd.units.den,'s'));
    if nden == 0,
      continue;
    elseif nden > 1,
      fprintf('Can''t fix %s easily, skipping for now\n',fn);
      tofix{end+1} = fn; %#ok<AGROW>
      continue;
    end
    fprintf('Fixing dt for %s, experiment %s\n',fn,experiment_name);
    for i = 1:numel(pd.data),
      pd.data{i} = pd.data{i} .* dt{i} ./ meddt;
    end
    pd.dtfixed = true;
    if ~DEBUG,
      delete(filename);
      save(filename,'-struct','pd');
    end
    
  end
  
%% replace dt.mat

fn = 'dt';
pd = struct;
pd.data = cell(1,numel(trx));
pd.units = parseunits('s');
pd.dtfixed = true;
for i = 1:numel(trx),
  pd.data{i} = repmat(meddt,[1,trx(i).nframes-1]);
  trx(i).dt = pd.data{i};
end
filename = fullfile(perframedir,[fn,'.mat']);
fprintf('Fixing dt.mat, experiment %s\n',experiment_name);
if ~DEBUG,
  if exist(filename,'file'),
    delete(filename);
  end
  save(filename,'-struct','pd');
end

%% re-save trajectory file

fprintf('Fixing trx file %s\n',trxfile);
td.trx = trx;
td.dtfixed = true;
if ~DEBUG,
  if exist(trxfile,'file'),
    delete(trxfile);
  end
  save(trxfile,'-struct','td');
end

end

%% remove features that need to be recomputed

if ~isempty(tofix),
  traj = trx;
  trx = Trx('analysis_protocol',analysis_protocol,'settingsdir',settingsdir,...
    'datalocparamsfilestr',datalocparamsfilestr);
  trx.AddExpDir(expdir,'openmovie',false,...
    'dooverwrite',false,'traj',traj);
  if ~DEBUG,
    trx.CleanPerFrameData(tofix);
  end
%   for i = 1:numel(tofix),
%     fn = tofix{i};
%     fprintf('Recomputing %s for %s\n',fn,expdir);
%     trx.(fn);
%   end
end

%% compute all missing per-frame features
FlyBowlComputePerFrameFeatures(expdir,'analysis_protocol',analysis_protocol,...
  'settingsdir',settingsdir,'datalocparamsfilestr',datalocparamsfilestr,...
  'forcecompute',false);
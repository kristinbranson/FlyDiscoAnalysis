function [isconsistent,filesmissing,expdirs] = SAGEExpDirConsistencyCheckMany(varargin)

[analysis_protocol,settingsdir,datalocparamsfilestr,outfilename,rootdatadir,restartexperiment,experiments,screen_type,leftovers] = ...
  myparse_nocheck(varargin,...
  'analysis_protocol','current',...
  'settingsdir','/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings',...
  'datalocparamsfilestr','dataloc_params.txt',...
  'outfilename','',...
  'rootdatadir','/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data',...
  'restartexperiment','',...
  'experiments',[],...
  'screen_type','primary');

global DISABLESAGEWAITBAR;
DISABLESAGEWAITBAR = true;
datetimeformat = 'yyyymmddTHHMMSS';

%% get all experiment directories

if isempty(experiments),

fprintf('Getting all experiment directories...\n');
[expdirs,~,~,experiments] = getExperimentDirs('rootdir',rootdatadir,...
  'settingsdir',settingsdir,'protocol',analysis_protocol,...
  leftovers{:});

fprintf('Reading cross dates and screen types...\n');
% get cross date and only find experiments with the right screen_type
ignoreidx = false(1,numel(expdirs));
badidx = false(1,numel(expdirs));
for i = 1:numel(expdirs),
  metadatafile = fullfile(rootdatadir,expdirs{i},'Metadata.xml');
  if ~exist(metadatafile,'file'),
    badidx(i) = true;
    continue;
  end
  try
    metadata = ReadMetadataFile(metadatafile);
  catch ME,
    warning('Could not parse metadata file %s:\n%s',metadatafile,getReport(ME));
    badidx(i) = true;
    continue;
  end
  if ~isfield(metadata,'cross_date'),
    warning('Cross date not found in metadata file %s',metadatafile);
    continue;
  end
  if ~isfield(metadata,'screen_type'),
    warning('Screen type not found in metadata file %s',metadatafile);
  else
    experiments(i).screen_type = metadata.screen_type;
    if ~isempty(screen_type) && ~strcmp(screen_type,metadata.screen_type),
      ignoreidx(i) = true;
      continue;
    end
  end
  cross_date = metadata.cross_date;
  experiments(i).cross_date = cross_date;
  if isfield(metadata,'wish_list'),
    experiments(i).wish_list = metadata.wish_list;
  else
    experiments(i).wish_list = [];
  end
  try
    crossdatenum = datenum(cross_date,datetimeformat);
  catch %#ok<*CTCH>
    try
      crossdatenum = datenum(cross_date);
    catch
      crossdatenum = nan;
      warning('Could not parse cross date %s',cross_date);
     end
  end
  experiments(i).cross_datenum = crossdatenum;
end

fprintf('Removing %d experiment directories that have no metadata file...\n',nnz(badidx&~ignoreidx));
fprintf('Removing %d experiment directories based on the screen type constraint ...\n',nnz(ignoreidx));
experiments(ignoreidx|badidx) = [];
expdirs(ignoreidx|badidx) = [];

else
  expdirs = {experiments.file_system_path};

end

fprintf('Sorting by cross date...\n');
[~,order] = sort([experiments.cross_datenum]);
expdirs = expdirs(order(end:-1:1));
experiments = experiments(order(end:-1:1));

%% check each experiment

lastcrossdatenum = nan;

isconsistent = nan(1,numel(expdirs));
filesmissing = nan(1,numel(expdirs));

if ~isempty(restartexperiment),
  experiment_names = cell(1,numel(experiments));
  for i = 1:numel(experiments),
    [~,experiment_names{i}] = fileparts(experiments(i).file_system_path);
  end
  starti = find(strcmp(experiment_names,restartexperiment),1);
  if isempty(starti),
    error('Could not find experiment %s to restart.',restartexperiment);
  end
  %start_cross_datenum = experiments(i).cross_datenum;
  %starti = find([experiments.cross_datenum] == start_cross_datenum,1);
else
  fid = myfopen(outfilename,'w');
  myfclose(fid);
  starti = 1;
end
  
for i = starti:numel(expdirs),  

  expdir = experiments(i).file_system_path;
  experiment = experiments(i);
  if isnan(experiment.cross_datenum),
    cross_date = '?';
  else
    cross_date = datestr(experiment.cross_datenum,'yyyy-mm-dd');
  end
  if experiment.cross_datenum ~= lastcrossdatenum,
    fid = myfopen(outfilename,'a');
    if isempty(experiment.wish_list),
      fprintf(fid,'\n\n  *** CROSS DATE %s *** \n',cross_date);
    else
      fprintf(fid,'\n\n *** WISH LIST %d, CROSS DATE %s *** \n',experiment.wish_list,cross_date);
    end
    lastcrossdatenum = experiment.cross_datenum;
    myfclose(fid);
  end
  fprintf('%s (cross date %s) ... \n',expdirs{i},cross_date);
  [isconsistent(i),filesmissing(i)] = SAGEExpDirConsistencyCheck(expdir,...
    'analysis_protocol',analysis_protocol,...
    'settingsdir',settingsdir,...
    'datalocparamsfilestr',datalocparamsfilestr,...
    'docheckhist',false,...
    'outfilename',outfilename);
  if isconsistent(i),
    fprintf('    consistent ');
  else
    fprintf('    INCONSISTENT ');
  end
  if filesmissing(i),
    fprintf('    FILESMISSING\n');
  else
    fprintf('    filesexist\n');
  end
end

% %% grab out only the inconsistent experiments
% 
% idx = find(~isconsistent);
% fid = fopen(outfilename,'r');
% fid2 = fopen('SAGEInconsisteExps20110907.txt','w');
% for i = idx(:)',
%   expdir = experiments(i).file_system_path;
%   isfirst = true;
%   didfind = false;
%   while true,
%     s = fgetl(fid);
%     if ~ischar(s),
%       if isfirst,
%         fseek(fid,0,'bof');
%         isfirst = false;
%       else
%         break;
%       end
%     end
%     if strcmp(s,expdir),
%       didfind = true;
%       break;
%     end
%   end
%   fprintf(fid2,'\n%s\n',expdir);
%   while true,
%     s = fgetl(fid);
%     if ~ischar(s) || isempty(s),
%       break;
%     end
%     fprintf(fid2,'%s\n',s);
%   end
% end
% fclose(fid2);
% fclose(fid);
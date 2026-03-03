function propagateMetadataFixes(infilename,outfilename,varargin)

% stuff to read and where it is
property_names = {'rig','plate','top_plate','bowl','camera','computer','harddrive','apparatus_id','file_system_path','effector','exp_datetime'};
property_types = {'session','session','session','session','session','session','session','session','experiment','session','experiment'};
cv_tables = {'fly_olympiad_apparatus','fly_olympiad_apparatus','fly_olympiad_apparatus','fly_olympiad_apparatus','fly_olympiad_apparatus','fly_olympiad_apparatus','fly_olympiad_apparatus','fly_olympiad','fly_olympiad','fly_olympiad','fly_olympiad'};
property_isnumeric = [true,true,true,false,false,false,false,false,false,false,false];

%% parse arguments

sage_params_path = '/groups/branson/home/bransonk/SAGEConfig/SAGEReadParams.txt';

[sage_params_path,rootdatadir,DEBUG] = ...
  myparse(varargin,...
  'sage_params_path',sage_params_path,...
  'rootdatadir','/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data',...
  'debug',0);

if ischar(DEBUG),
  DEBUG = str2double(DEBUG);
end

%% read in the experiments to fix

infid = fopen(infilename,'r');
if infid < 0,
  error('Could not open %s to read',infilename);
end
experiment_names = {};
while true,
  s = fgetl(infid);
  if ~ischar(s),
    break;
  end
  s = strtrim(s);
  if isempty(s),
    continue;
  end
  if s(1) == '#',
    continue;
  end
  experiment_names{end+1} = s; %#ok<AGROW>
end
fclose(infid);

%% connect to sage

warning off MATLAB:javaclasspath:jarAlreadySpecified
classpath = getenv('CLASSPATH');
for path = regexp(classpath, ':', 'split')
  javaaddpath(path);
end

% Read the parameter file.
fid = fopen(sage_params_path);
line_synonyms = {};
try
  params = strtrim(textscan(fid, '%s %s'));
  host = params{2}{strmatch('host:', params{1})}; %#ok<*REMFF1>
  db_name = params{2}{strmatch('database:', params{1})};
  user_name = params{2}{strmatch('username:', params{1})};
  password = params{2}{strmatch('password:', params{1})};
  for index = strmatch('line-synonym:', params{1})'
    syn_count = size(line_synonyms, 1);
    parts = regexp(params{2}{index}, '=', 'split');
    line_synonyms{syn_count + 1, 1} = parts{1}; %#ok
    line_synonyms{syn_count + 1, 2} = parts{2}; %#ok
  end
  fclose(fid);
catch ME
  fclose(fid);
  rethrow(ME);
end

% connect
db = database(db_name, user_name, password, 'com.mysql.jdbc.Driver', ['jdbc:mysql://' host '/' db_name]);

if ~isempty(db.Message)
  if strcmp(db.Message, 'JDBC Driver Error: com.mysql.jdbc.Driver. Driver Not Found/Loaded.')
    message = [db.Message char(10) char(10) 'Make sure you add the path to the MySQL Connector JAR file is in the CLASSPATH.'];
  else
    message = db.Message;
  end
  error(['Could not create database connection: ' message]);
end
% Make sure the experiment is loaded atomically.
set(db, 'AutoCommit', 'off');

%% initialize outputs

doclose = false;
if DEBUG >= 2 || isempty(outfilename),
  outfid = 1;
else
  doclose = true;
  outfid = fopen(outfilename,'w');
end

%% loop over all experiments

%% read data from SAGE
for i = 1:numel(experiment_names),
  
  experiment_name = experiment_names{i};
  
  fprintf('Processing %s...\n',experiment_name);
  
  data = struct;
  data.experiment_name = experiment_name;
  
  % get experiment id
  query = sprintf('select id from experiment where name = ''%s''',experiment_name);
  try
    curs = exec(db, query);
    curs = fetch(curs);
    close(curs);
    data.experiment_id = curs.Data{1};
    if ischar(data.experiment_id) && strcmp(data.experiment_id,'No Data'),
      error('Could not find experiment %s\n',experiment_name);
    end
  catch ME,
    getReport(ME)
    continue;
  end
  
  % get experiment stuff
  % query = sprintf('select * from experiment
  
  % get session id
  query = sprintf('select id from session where experiment_id = %d',data.experiment_id);
  try
    curs = exec(db, query);
    curs = fetch(curs);
    close(curs);
    data.session_id = curs.Data{1};
    if ischar(data.session_id) && strcmp(data.session_id,'No Data'),
      error('Could not get session id for %s\n',experiment_name);
    end
  catch ME,
    getReport(ME)
    continue;
  end
  
  % get line
  query = sprintf('select line_id from session where id = %d',data.session_id);
  try
    curs = exec(db, query);
    curs = fetch(curs);
    close(curs);
    data.line_id = curs.Data{1};
    if ischar(data.line_id) && strcmp(data.line_id,'No Data'),
      error('Could not get line id for %s\n',experiment_name);
    end
  catch ME,
    getReport(ME)
    continue;
  end
  query = sprintf('select name from line_vw where id = %d',data.line_id);
  try
    curs = exec(db, query);
    curs = fetch(curs);
    close(curs);
    data.line_name = curs.Data{1};
    if ischar(data.line_name) && strcmp(data.line_name,'No Data'),
      error('Could not get line_name for %s\n',experiment_name);
    end
  catch ME,
    getReport(ME)
    continue;
  end
  
  try
    for j = 1:numel(property_names),
      if strcmp(property_types{j},'session'),
        id = data.session_id;
      else
        id = data.experiment_id;
      end
      query = sprintf('select value from %s_property where %s_id = %d and type_id = %s',...
        property_types{j},property_types{j},id,Sage_cv_term(property_names{j},cv_tables{j}));
      curs = exec(db, query);
      curs = fetch(curs);
      close(curs);
      data.(property_names{j}) = curs.Data{1};
      if ischar(data.(property_names{j})) && strcmp(data.(property_names{j}),'No Data'),
        error('Could not get %s for %s\n',property_names{j},experiment_name);
      end
      if property_isnumeric(j),
        data.(property_names{j}) = str2double(data.(property_names{j}));
      end
    end
  catch ME,
    getReport(ME)
    continue;
  end

  %% make apparatus_id match individual fields
  
  apparatus_id = sprintf('Rig%d__Plate%02d__Lid%02d__Bowl%s__Camera%s__Computer%s__HardDrive%s',...
    data.rig,data.plate,data.top_plate,data.bowl,data.camera,data.computer,data.harddrive);
  
  if ~strcmp(apparatus_id,data.apparatus_id),
    
    [~,~,success] = updateSageProperty(data,db,'apparatus_id',apparatus_id,'session',DEBUG);
    
    if ~success,
      warning('Could not update apparatus id');
    end
    
  end
  
  %% make file system path match individual fields
  
  expname = GetExperimentName(data);
  file_system_path = sprintf('%s/%s',rootdatadir,expname);
  if ~strcmp(file_system_path,data.file_system_path),
    [~,~,success] = updateSageProperty(data,db,'file_system_path',file_system_path,'experiment',DEBUG);
    if ~success,
      warning('Could not update file_system_path');
    end
    
    % rename the experiment directory if necessary
    curr_file_system_path = data.file_system_path;
    if exist(data.file_system_path,'file') && ~exist(file_system_path,'file'),
      cmd = sprintf('mv %s/ %s/',data.file_system_path,file_system_path);
      if DEBUG > 0,
        fprintf('** TODO: Execute:\n  %s\n',cmd);
      else
        unix(cmd);
        curr_file_system_path = file_system_path;
      end
    end
    
    % tell user to create the ctrax results movie for this file
    [~,oldexpname] = fileparts(data.file_system_path);
    oldmoviename = fullfile(curr_file_system_path,sprintf('ctrax_results_movie_%s.avi',oldexpname));
    newmoviename = fullfile(curr_file_system_path,sprintf('ctrax_results_movie_%s.avi',expname));
    if exist(oldmoviename,'file') && ~exist(newmoviename,'file'),
      fprintf(outfid,'%s\n',file_system_path);
    end
  end

  %% make experiment name match individual fields
  
  experiment_name = ['FlyBowl_',expname];
  if ~strcmp(data.experiment_name,experiment_name);
    
    [~,~,success] = updateSageExperimentName(data,db,experiment_name,DEBUG);
    if ~success,
      warning('Could not update experiment_name');
    end

  end
  
end

if doclose,
  fclose(outfid);
end

commit(db);
close(db);

function [query,comments,success] = updateSageProperty(data,db,property_name,property_value,property_type,DEBUG)

comments = '';
success = false;

id = data.(sprintf('%s_id',property_type));

query = sprintf('select id from %s_property where %s_id = %d and type_id = %s',property_type,property_type,id,Sage_cv_term(property_name));
if isempty(db),
  fprintf([query,'\n']);
  property_id = nan;
else
  curs = exec(db, query);
  curs = fetch(curs);
  close(curs);
  property_id = curs.Data{1};
end

if strcmp(property_id,'No Data'),
  warning('%s not yet set for experiment %s',property_name,data.experiment_name);
  query = '';
  return;
end
comments = sprintf('Update %s for experiment %s to %s',property_name,data.experiment_name,property_value);
query = sprintf('update %s_property set value = ''%s'' where id = %d',property_type,escape_string(property_value),property_id);
fprintf(['#',comments,'\n']);
fprintf([query,'\n']);
if DEBUG == 0,
  curs = exec(db, query);
  if ~isempty(curs.Message)
    error('Could not update %s for experiment %s to %s',property_name,data.experiment_name,property_value);
  end
  close(curs);
end
success = true;

function [query,comments,success] = updateSageExperimentName(data,db,property_value,DEBUG)

comments = sprintf('Update name for experiment %s to %s',data.experiment_name,property_value);
query = sprintf('update experiment set name = ''%s'' where id = %d',property_value,data.experiment_id);
fprintf(['#',comments,'\n']);
fprintf([query,'\n']);
if DEBUG == 0,
  curs = exec(db, query);
  if ~isempty(curs.Message)
    error('Could not update %s for experiment %s to %s',property_name,data.experiment_name,property_value);
  end
  close(curs);
end
success = true;


function expname = GetExperimentName(data)
  
effector_abbrs = {'UAS_dTrpA1_2_0002','TrpA'
  'EXT_CantonS_1101243','CS'
  'NoEffector_0_9999','None'};

i = find(strcmp(data.effector,effector_abbrs(:,1)),1);
if ~isempty(i),
  effector = effector_abbrs{i,2};
else
  effector = data.effector;
end
timestampstr = data.exp_datetime;

expname = sprintf('%s_%s_Rig%dPlate%02dBowl%s_%s',...
  data.line_name,...
  effector,...
  data.rig,data.plate,data.bowl,...
  timestampstr);


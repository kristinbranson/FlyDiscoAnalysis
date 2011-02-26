function db = ConnectToSage(sage_params_path)

% set up to connect to database
warning off MATLAB:javaclasspath:jarAlreadySpecified
classpath = getenv('CLASSPATH');
for path = regexp(classpath, ':', 'split')
  javaaddpath(path);
end

%% Read the parameter file.
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

%% Connect to SAGE
db = database(db_name, user_name, password, 'com.mysql.jdbc.Driver', ['jdbc:mysql://' host '/' db_name]);

if ~isempty(db.Message)
  if strcmp(db.Message, 'JDBC Driver Error: com.mysql.jdbc.Driver. Driver Not Found/Loaded.')
    message = [db.Message char(10) char(10) 'Make sure you the path to the MySQL Connector JAR file is in the CLASSPATH.'];
  else
    message = db.Message;
  end
  error(['Could not create database connection: ' message]);
end
% Make sure the experiment is loaded atomically.
set(db, 'AutoCommit', 'off');

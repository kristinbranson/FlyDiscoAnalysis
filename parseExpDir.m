function [res,success] = parseExpDir(expdir,sagecompatible)

if nargin < 2,
  sagecompatible = false;
end

effectors = {'TrpA','CTRL_CantonS_1101243','CTRL_CantonS_1500002_0029'};
expr = ['^(?<pathstr>(.*[/\\])*)(?<line>[^/^\\]+)_(?<effector>(',sprintf('|(%s)',effectors{:}),'))_Rig(?<rig>[0-9]+)Plate(?<plate>[0-9]+)Bowl(?<bowl>[A-Z]+)_(?<notstarted>(notstarted_)?)(?<date>\d{8}T\d{6}).*$'];
res = regexp(expdir,expr,'names','once');
success = ~isempty(res);

if ~success,
  expr = ['^(?<pathstr>(.*[/\\])*)(?<line>[^/^\\]+)_(?<effector>(([^_]*)',sprintf('|(%s)',effectors{:}),'))_Rig(?<rig>[0-9]+)Plate(?<plate>[0-9]+)Bowl(?<bowl>[A-Z]+)_(?<notstarted>(notstarted_)?)(?<date>\d{8}T\d{6}).*$'];
  res = regexp(expdir,expr,'names','once');
  success = ~isempty(res);
end
  
if ~success,
  expr = '^(?<pathstr>(.*[/\\])*)(?<screen_type>[a-zA-Z0-9]+)_(?<screen_reason>\w+)_Rig(?<rig>[0-9]+)Bowl(?<bowl>[A-Z]+)_(?<notstarted>(notstarted_)?)(?<date>\d{8}T\d{6})$';
  res = regexp(expdir,expr,'names','once');
  success = ~isempty(res);
  if success,
    res.type = 'condition2';
  end
end
if ~success,
  expr = '^(?<pathstr>(.*[/\\])*)(?<screen_type>[a-zA-Z0-9]+)_(?<screen_reason>\w+)_(?<notstarted>(notstarted_)?)(?<date>\d{8}T\d{6})$';
  res = regexp(expdir,expr,'names','once');
  success = ~isempty(res);
  if success,
    res.type = 'condition1';
  end
end

if success,
  res.notstarted = strcmp(res.notstarted,'notstarted_');
end

if sagecompatible,
  res1 = struct;
  res1.file_system_path = expdir;
  res1.line_name = res.line;
  res1.exp_datetime = res.date;
  [~,name] = fileparts(expdir);
  res1.experiment_name = ['FlyBowl_',name];
  res1.bowl = res.bowl;
  res1.rig = str2double(res.rig);
  res1.plate = str2double(res.plate);
  res1.flag_aborted = res.notstarted;
  res = res1;
end
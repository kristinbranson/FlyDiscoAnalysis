function [success,msgs] = CheckExperiment20120726(expdir,varargin)

success = true;
msgs = {};

[analysis_protocol,settingsdir,datalocparamsfilestr,...
  mindatestr_retrack,maxdatestr_retrack,...
  mintimestamp_retrack,mintimestamp_fix] = ...
  myparse(varargin,...
  'analysis_protocol','current',...
  'settingsdir','/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings',...
  'datalocparamsfilestr','dataloc_params.txt',...
  'mindatestr_retrack','20110601T000000',...
  'maxdatestr_retrack','20120211T000000',...
  'mintimestamp_retrack','20120301T000000',...
  'mintimestamp_fix','20120701T000000');

datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
dataloc_params = ReadParams(datalocparamsfile);

mindatenum_retrack = datenum(mindatestr_retrack,'yyyymmddTHHMMSS');
maxdatenum_retrack = datenum(maxdatestr_retrack,'yyyymmddTHHMMSS');
mintimestampnum_retrack = datenum(mintimestamp_retrack,'yyyymmddTHHMMSS');
mintimestampnum_fix = datenum(mintimestamp_fix,'yyyymmddTHHMMSS');

files_retrack = {'registered_trx.mat','ctrax_results.mat','movie.ufmf.ann',...
  'ctrax_diagnostics.txt','registrationdata.mat','registrationdata.txt',...
  'sexclassifier_diagnostics.txt'};
files_fix = {'registered_trx.mat'};

perframefns_fix = {...
  'closestfly_nose2ell_angle_min30to30','dnose2ell_nose2ell_angle_min30to30',...
  'closestfly_nose2ell_angle_min20to20','dnose2ell_nose2ell_angle_min20to20',...
  'closestfly_nose2ell_angle_30tomin30','dnose2ell_nose2ell_angle_30tomin30',...
  'absdangle2wall','absdtheta','absdv_cor','abssmoothdtheta','dangle2wall',...
  'danglesub','ddcenter','ddell2nose','ddist2wall','ddnose2ell','dtheta',...
  'du_cor','du_ctr','dv_cor','dv_ctr','flipdv_cor','velmag','velmag_ctr',...
  'magveldiff_nose2ell','magveldiff_anglesub','veltoward_nose2ell',...
  'veltoward_anglesub','da','db','darea','decc','dphi','du_tail','dv_tail',...
  'velmag_nose','velmag_tail','accmag'};

%% check for experiment

if ~exist(expdir,'dir'),
  msg = sprintf('Directory %s does not exist',expdir);
  fprintf([msg,'\n']);
  success = false;
  msgs{end+1} = msg;
  return;
end

%% check registered_trx

filename = fullfile(expdir,dataloc_params.trxfilestr);
if ~exist(filename,'file'),
  msg = sprintf('Trx file %s does not exist',filename);
  fprintf([msg,'\n']);
  success = false;
  msgs{end+1} = msg;
else
  
  tmp = load(filename);
  n = numel(unique([tmp.trx.dt]));
  if n > 1,
    msg = sprintf('More than one unique value for dt in %s',filename);
    fprintf([msg,'\n']);
    success = false;
    msgs{end+1} = msg;
  end
  
end

%% check dt

filename = fullfile(expdir,dataloc_params.perframedir,'dt.mat');
if ~exist(filename,'file'),
  msg = sprintf('dt file %s does not exist',filename);
  fprintf([msg,'\n']);
  success = false;
  msgs{end+1} = msg;
else

  tmp = load(filename);
  
  if ~iscell(tmp.data),
    msg = sprintf('dt data is not a cell');
    fprintf([msg,'\n']);
    success = false;
    msgs{end+1} = msg;
  else
    n = numel(unique([tmp.data{:}]));
    if n > 1,
      msg = sprintf('More than one unique value for dt in %s',filename);
      fprintf([msg,'\n']);
      success = false;
      msgs{end+1} = msg;
    end
  end
end

%% check that a_mm has enough elements

filename = fullfile(expdir,dataloc_params.perframedir,'dt.mat');
if ~exist(filename,'file'),
  msg = sprintf('a_mm file %s does not exist',filename);
  fprintf([msg,'\n']);
  success = false;
  msgs{end+1} = msg;
else
  tmp = load(filename);
  if ~iscell(tmp.data),
    msg = sprintf('a_mm data is not a cell');
    fprintf([msg,'\n']);
    success = false;
    msgs{end+1} = msg;
  end
end

%% check timestamps 

exp_datetime = regexp(expdir,'\d{8}T\d{6}$','once','match');
exp_datenum = datenum(exp_datetime,'yyyymmddTHHMMSS');

retrack = exp_datenum >= mindatenum_retrack & exp_datenum <= maxdatenum_retrack;
if retrack,
  files = files_retrack;
  mintimestampnum = mintimestampnum_retrack;
  mintimestamp = mintimestamp_retrack;
  
  for i = 1:numel(files),
    filename = fullfile(expdir,files{i});
    tmp = dir(filename);
    if isempty(tmp),
      msg = sprintf('File %s does not exist',filename);
      fprintf([msg,'\n']);
      success = false;
      msgs{end+1} = msg;
    else
      if tmp.datenum < mintimestampnum,
        msg = sprintf('File %s has timestamp %s < %s',filename,tmp.date,mintimestamp);
        fprintf([msg,'\n']);
        success = false;
        msgs{end+1} = msg;
      end
    end
  end

end

files = files_fix;
mintimestampnum = mintimestampnum_fix;
mintimestamp = mintimestamp_fix;

for i = 1:numel(files),
  filename = fullfile(expdir,files{i});
  tmp = dir(filename);
  if isempty(tmp),
    msg = sprintf('File %s does not exist',filename);
    fprintf([msg,'\n']);
    success = false;
    msgs{end+1} = msg;
  else
    if tmp.datenum < mintimestampnum,
      msg = sprintf('File %s has timestamp %s < %s',filename,tmp.date,mintimestamp);
      fprintf([msg,'\n']);
      success = false;
      msgs{end+1} = msg;
    end
  end
end

  
%% check that all the per-frame features are there

perframefnsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.perframefnsfilestr);
perframefns = importdata(perframefnsfile);
nfns = numel(perframefns);
isfixfn = ismember(perframefns,perframefns_fix);

for i = 1:nfns,
  filename = fullfile(expdir,dataloc_params.perframedir,[perframefns{i},'.mat']);
  tmp = dir(filename);
  if isempty(tmp),
    msg = sprintf('File %s does not exist',filename);
    fprintf([msg,'\n']);
    success = false;
    msgs{end+1} = msg;
  elseif isfixfn(i),
    if tmp.datenum < mintimestampnum_fix,
      msg = sprintf('File %s has timestamp %s < %s',filename,tmp.date,mintimestamp_fix);
      fprintf([msg,'\n']);
      success = false;
      msgs{end+1} = msg;
    end
  elseif retrack,
    if tmp.datenum < mintimestampnum_retrack,
      msg = sprintf('File %s has timestamp %s < %s',filename,tmp.date,mintimestamp_retrack);
      fprintf([msg,'\n']);
      success = false;
      msgs{end+1} = msg;
    end    
  end
end


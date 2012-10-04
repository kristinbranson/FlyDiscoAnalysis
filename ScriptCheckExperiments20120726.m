% ScriptCheckExperiments20120726

addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
addpath /groups/branson/bransonlab/projects/olympiad/SAGE/MATLABInterface/Trunk;
addpath /groups/branson/bransonlab/projects/olympiad/FlyBowlDataCapture;
settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
rootdatadir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data';

analysis_protocol = '20120706';
mintimestamp = '20110201T000000';


%% find experiments that are new enough

experiment_names = dir(fullfile(rootdatadir,'*201*'));
experiment_names = {experiment_names.name};
timestamps = regexp(experiment_names,'\d{8}T\d{6}$','once','match');

datenums = datenum(timestamps,'yyyymmddTHHMMSS');
idx = datenums>datenum(mintimestamp,'yyyymmddTHHMMSS');
experiment_names = experiment_names(idx);

%% remove anything where the screen type is not primary

nonprimary = SAGEListBowlExperiments('screen_type',{'*non*','*albin*','*quest*'},'checkflags',false,'removemissingdata',false);

nonprimary_names = regexprep({nonprimary.experiment_name},'FlyBowl_','');
experiment_names = setdiff(experiment_names,nonprimary_names);
timestamps = regexp(experiment_names,'\d{8}T\d{6}$','once','match');
[timestamps,order] = sort(timestamps);
experiment_names = experiment_names(order);
datenums = datenum(timestamps,'yyyymmddTHHMMSS');


%% initialize

issuccess = nan(1,numel(experiment_names));
msgs = cell(1,numel(experiment_names));
nproblems = 0;
ii = 0;

hfig = 1;
figure(hfig);
clf;

hdata = plot(datenums,issuccess,'k.-');
set(gca,'XLim',[min(datenums)-1,max(datenums)+1],'YLim',[-.25,1.25]);
datetick('x','mmmyy','keeplimits');
hti = title(sprintf('%d problems / %d checked',nproblems,ii));

order = randperm(numel(experiment_names));

%%

for ii = 1:npar:numel(experiment_names),
  
  %i = order(ii);
  
  idxcurr = order(ii:min(ii+npar-1,numel(experiment_names)));
  currissuccess = nan(1,numel(idxcurr));
  currmsgs = cell(1,numel(idxcurr));
  curr_experiment_names = experiment_names(idxcurr);
  
  parfor j = 1:numel(idxcurr),
    fprintf('***** Checking experiment %s...\n',curr_experiment_names{j});
    [currissuccess(j),currmsgs{j}] = CheckExperiment20120726(fullfile(rootdatadir,curr_experiment_names{j}),'analysis_protocol',analysis_protocol);
  end
  
  issuccess(idxcurr) = currissuccess;
  msgs(idxcurr) = currmsgs;
  nproblems = nnz(issuccess==0);
  
  set(hdata,'YData',issuccess);
  set(hti,'String',sprintf('%d problems / %d checked',nproblems,ii+numel(idxcurr)-1));
  
  drawnow;
  
  if mod(ii,100) < npar,
    save CheckExperiments20120726.mat issuccess msgs order ii experiment_names;
  end
  
end

%%

for i = idxworrisome,
  [issuccess(i),msgs{i}] = CheckExperiment20120726(fullfile(rootdatadir,experiment_names{i}),'analysis_protocol',analysis_protocol);
end

%% ignore things with automated_checks_incoming = F

isworrisome = issuccess == 0;
for i = find(issuccess==0),
  if ~isempty(regexp(experiment_names{i},'CantonS','once')),
    isworrisome(i) = false;
    continue;
  end
  expdir = fullfile(rootdatadir,experiment_names{i});
  if ~exist(fullfile(expdir,'movie.ufmf'),'file') || ~exist(fullfile(expdir,'Metadata.xml'),'file'),
    isworrisome(i) = false;
    continue;
  end
  filename = fullfile(expdir,'automatic_checks_incoming_results.txt');
  if ~exist(filename,'file'),
    isworrisome(i) = false;
    continue;
  end
  try
    params = ReadParams(filename);
    if strcmpi(params.automated_pf,'F'),
      isworrisome(i) = false;
      continue;
    end
  catch ME,
    warning(getReport(ME));
  end
  filename = fullfile(expdir,'ctrax_log.txt');
  [~,res] = unix(sprintf('grep ShortUFMFFileError %s',filename));
  if ~isempty(strtrim(res)),
    isworrisome(i) = false;
    continue;
  end
  
end

%%

tmp = find(isworrisome);
for i = tmp(randperm(numel(tmp))),
  fprintf('%s:\n',experiment_names{i});
  fprintf('%s\n',msgs{i}{:});
  input('');
end

%% 

idxworrisome = find(isworrisome);
[~,order] = sort(datenums(isworrisome));
idxworrisome = idxworrisome(order);

newexp = datenums(idxworrisome)' > datenum('20120701T000000','yyyymmddTHHMMSS');
badtimestamp = cellfun(@(x) any(~cellfun(@isempty,regexp(x,'has timestamp.*< 20120701T000000','once'))),msgs(idxworrisome));
badtrajfns = cellfun(@(x) any(~cellfun(@isempty,regexp(x,'((dt)|(a_mm)) data is not a cell','once'))),msgs(idxworrisome));
dtnotfixed = cellfun(@(x) any(~cellfun(@isempty,regexp(x,'More than one unique value for dt in ','once'))),msgs(idxworrisome));
notretracked = cellfun(@(x) any(~cellfun(@isempty,regexp(x,'has timestamp.*< 20120301T000000','once'))),msgs(idxworrisome));
missingtrx = cellfun(@(x) any(~cellfun(@isempty,regexp(x,'File.*registered_trx.*does not exist','once'))),msgs(idxworrisome));
missingotherfiles = ~missingtrx & cellfun(@(x) any(~cellfun(@isempty,regexp(x,'File.*does not exist','once'))),msgs(idxworrisome));

timestamps_ctrax = cell(1,numel(idxworrisome));
timestamps_reg = cell(1,numel(idxworrisome));
timestamps_sex = cell(1,numel(idxworrisome));
timestamps_velmag = cell(1,numel(idxworrisome));
timestamps_anglerange = cell(1,numel(idxworrisome));

for ii = 1:numel(idxworrisome),
  i = idxworrisome(ii);
 filename = fullfile(rootdatadir,experiment_names{i},'ctrax_results.mat');
  tmp = dir(filename);
  if isempty(tmp),
    timestamps_ctrax{ii} = 'MISSING';
  else
    timestamps_ctrax{ii} = datestr(tmp.datenum,'yyyymmdd');
  end
  
  filename = fullfile(rootdatadir,experiment_names{i},'registered_trx.mat');
  tmp = dir(filename);
  if isempty(tmp),
    timestamps_reg{ii} = 'MISSING';
  else
    timestamps_reg{ii} = datestr(tmp.datenum,'yyyymmdd');
  end
  
  filename = fullfile(rootdatadir,experiment_names{i},'sexclassifier.mat');
  tmp = dir(filename);
  if isempty(tmp),
    timestamps_sex{ii} = 'MISSING';
  else
    timestamps_sex{ii} = datestr(tmp.datenum,'yyyymmdd');
  end

  filename = fullfile(rootdatadir,experiment_names{i},'perframe','velmag.mat');
  tmp = dir(filename);
  if isempty(tmp),
    timestamps_velmag{ii} = 'MISSING';
  else
    timestamps_velmag{ii} = datestr(tmp.datenum,'yyyymmdd');
  end

  filename = fullfile(rootdatadir,experiment_names{i},'perframe','dnose2ell_angle_min20to20.mat');
  tmp = dir(filename);
  if isempty(tmp),
    timestamps_anglerange{ii} = 'MISSING';
  else
    timestamps_anglerange{ii} = datestr(tmp.datenum,'yyyymmdd');
  end

end

%%

experiments_worrisome = experiment_names(idxworrisome);

expdirs_retrack = importdata('expdirs_retrack_20120705.txt');
expdirs_unarchive = importdata('expdirs_movefromarchive_20120710.txt');
expdirs_fix = importdata('expdirs_fixperframefeatures_20120710.txt');

experiments_unarchive = cell(size(expdirs_unarchive));
for i = 1:numel(expdirs_unarchive),
  [~,experiments_unarchive{i}] = fileparts(expdirs_unarchive{i});
end
experiments_retrack = cell(size(expdirs_retrack));
for i = 1:numel(expdirs_retrack),
  [~,experiments_retrack{i}] = fileparts(expdirs_retrack{i});
end
experiments_fix = cell(size(expdirs_fix));
for i = 1:numel(expdirs_fix),
  [~,experiments_fix{i}] = fileparts(expdirs_fix{i});
end

isnonprimary = ismember(experiments_worrisome,nonprimary_names);
onarchivelist = ismember(experiments_worrisome,experiments_unarchive);
onretracklist = ismember(experiments_worrisome,experiments_retrack);
onfixlist = ismember(experiments_worrisome,experiments_fix);


%%

fid = fopen('CheckExperiments20120730.csv','w');
timestampsworrisome = cellstr(datestr(datenums(idxworrisome),'yyyymmddTHHMMSS'));

fprintf(fid,'Timestamp,Experiment name,Missing trx file,Not retracked,Not fixed,Missing other files,dt or a_mm bug,Bad timestamp,New experiment,Ctrax timestamp,Register timestamp,Sex class timestamp,Velmag timestamp,Anglerange timestamp,Non-primary,Unarchive list,Retrack list,Fix list,Notes\n');
for ii = 1:numel(idxworrisome),
  i = idxworrisome(ii);
  fprintf(fid,'%s,%s,%d,%d,%d,%d,%d,%d,%d,%s,%s,%s,%s,%s,%d,%d,%d,%d,',timestampsworrisome{ii},experiment_names{i},missingtrx(ii),notretracked(ii),dtnotfixed(ii),missingotherfiles(ii),badtrajfns(ii),badtimestamp(ii),newexp(ii),...
    timestamps_ctrax{ii},timestamps_reg{ii},timestamps_sex{ii},timestamps_velmag{ii},timestamps_anglerange{ii},...
    isnonprimary(ii),onarchivelist(ii),onretracklist(ii),onfixlist(ii));
  fprintf(fid,'%s   ',msgs{i}{:});
  fprintf(fid,'\n');
end  

fclose(fid);

%% check for no timestamps

badtimestampsintrx = false(1,numel(idxworrisome));
for ii = 1:numel(idxworrisome),
  i = idxworrisome(ii);
  filename = fullfile(rootdatadir,experiment_names{i},'registered_trx.mat');
  if ~exist(filename,'file'),
    fprintf('%s: no registered_trx.mat\n',experiment_names{i});
    continue;
  end
  td = load(filename,'timestamps');
  if isempty(td.timestamps),
    fprintf('%s: empty timestamps\n',experiment_names{i});
    badtimestampsintrx(ii) = true;
  end
end

%% 

idx_fixerror = ~isnonprimary & ~onarchivelist & onfixlist;
fprintf('%s\n',experiments_worrisome{idx_fixerror});

idx_retrackerror = ~isnonprimary & ~onarchivelist & onretracklist;
fprintf('%s\n',experiments_worrisome{idx_retrackerror});

tmp = importdata('List_RetrackMissing.txt');
tmpmetadata = SAGEListBowlExperiments('experiment_name',cellfun(@(x) ['FlyBowl_',x],tmp,'UniformOutput',false),'checkflags',false)

%%

fid = fopen('CheckExperiments20120726_unexplained.csv','w');
timestampsworrisome = cellstr(datestr(datenums(idxworrisome),'yyyymmddTHHMMSS'));

fprintf(fid,'Timestamp,Experiment name,Missing trx file,Not retracked,Not fixed,Missing other files,dt or a_mm bug,Bad timestamp,New experiment,Ctrax timestamp,Register timestamp,Sex class timestamp,Velmag timestamp,Anglerange timestamp,Non-primary,Unarchive list,Retrack list,Fix list,Notes\n');
for ii = find(idx_fixerror | idx_retrackerror),
  i = idxworrisome(ii);
  fprintf(fid,'%s,%s,%d,%d,%d,%d,%d,%d,%d,%s,%s,%s,%s,%s,%d,%d,%d,%d,',timestampsworrisome{ii},experiment_names{i},missingtrx(ii),notretracked(ii),dtnotfixed(ii),missingotherfiles(ii),badtrajfns(ii),badtimestamp(ii),newexp(ii),...
    timestamps_ctrax{ii},timestamps_reg{ii},timestamps_sex{ii},timestamps_velmag{ii},timestamps_anglerange{ii},...
    isnonprimary(ii),onarchivelist(ii),onretracklist(ii),onfixlist(ii));
  fprintf(fid,'%s   ',msgs{i}{:});
  fprintf(fid,'\n');
end  

fclose(fid);

%%

expdirs_retrack = importdata('expdirs_retrack_20120705.txt');
expdirs_unarchive = importdata('expdirs_movefromarchive_20120710.txt');
expdirs_fix = importdata('expdirs_fixperframefeatures_20120710.txt');

experiments_unarchive = cell(size(expdirs_unarchive));
for i = 1:numel(expdirs_unarchive),
  [~,experiments_unarchive{i}] = fileparts(expdirs_unarchive{i});
end
experiments_retrack = cell(size(expdirs_retrack));
for i = 1:numel(expdirs_retrack),
  [~,experiments_retrack{i}] = fileparts(expdirs_retrack{i});
end
experiments_fix = cell(size(expdirs_fix));
for i = 1:numel(expdirs_fix),
  [~,experiments_fix{i}] = fileparts(expdirs_fix{i});
end

experiments_worrisome = experiment_names(idxworrisome);

%% 

experiments_new = experiments_worrisome(~cellfun(@isempty,regexp(experiments_worrisome,'201207\d{2}T\d{6}$','once')));
tmp = setdiff(experiments_worrisome,experiments_unarchive);
tmp = setdiff(tmp,nonprimary_names);
tmp = setdiff(tmp,experiments_new);
tmp = intersect(tmp,union(experiments_fix,experiments_retrack))

%% choose experiments to analyze

analysis_protocols = cell(size(experiment_names));
for i = 1:numel(experiment_names),
  if mod(i,100) == 0,
    fprintf('%d / %d\n',i,numel(experiment_names));
  end
  filename = fullfile(rootdatadir,experiment_names{i},'analysis_protocol.txt');
  didread = false;
  if exist(filename,'file'),
    try
      fid = fopen(filename,'r');
      analysis_protocols{i} = fgetl(fid);
      didread = true;
      fclose(fid);
    catch ME
      error(getReport(ME));
    end
  end
  if ~didread,
    tmp = SAGEListBowlExperiments('experiment_name',['FlyBowl_',experiment_names{i}],'checkflags',false);
    if ~isempty(tmp),
      analysis_protocols{i} = tmp.analysis_protocol;
    else
      analysis_protocols{i} = '???';
    end
  end
end

nanalyze_complete_date = 100;
nanalyze_complete_perprotocol = 2;
nanalyze_incomplete_date = 500;
nanalyze_incomplete_perprotocol = 2;

tmp = find(issuccess);
idx_complete_date = tmp(round(linspace(1,numel(tmp),nanalyze_complete_date)));
[unique_protocols,~,idxprotocol] = unique(analysis_protocols);
idx_complete_protocol = [];
for i = 1:numel(unique_protocols),
  idxcurr = find(idxprotocol == i & issuccess);
  if isempty(idxcurr),
    continue;
  end
  if numel(idxcurr) <= nanalyze_complete_perprotocol,
    idx_complete_protocol = [idx_complete_protocol,idxcurr];
  else
    idx_complete_protocol = [idx_complete_protocol,idxcurr(round(linspace(1,numel(idxcurr),nanalyze_complete_perprotocol)))];
  end
end

idx_complete = union(idx_complete_protocol,idx_complete_date);

fid = fopen('expdirs_check_complete_20120810.txt','w');
for i = 1:numel(idx_complete),
  fprintf(fid,'%s\n',fullfile(rootdatadir,experiment_names{idx_complete(i)}));
end
fclose(fid);

tmp = setdiff(find(issuccess),idx_complete);
idx_incomplete_date = tmp(round(linspace(1,numel(tmp),nanalyze_incomplete_date)));
tmpallowed = issuccess & ~ismember(1:numel(issuccess),idx_complete);
idx_incomplete_protocol = [];
for i = 1:numel(unique_protocols),
  idxcurr = find(idxprotocol == i & tmpallowed);
  if isempty(idxcurr),
    continue;
  end
  if numel(idxcurr) <= nanalyze_incomplete_perprotocol,
    idx_incomplete_protocol = [idx_incomplete_protocol,idxcurr];
  else
    idx_incomplete_protocol = [idx_incomplete_protocol,idxcurr(round(linspace(1,numel(idxcurr),nanalyze_incomplete_perprotocol)))];
  end
end

idx_incomplete = union(idx_incomplete_date,idx_incomplete_protocol);

fid = fopen('expdirs_check_incomplete_20120810.txt','w');
for i = 1:numel(idx_incomplete),
  fprintf(fid,'%s\n',fullfile(rootdatadir,experiment_names{idx_incomplete(i)}));
end
fclose(fid);

%% clone directories

tmpdatadir = 'Check20120810';
EPSILON = .01;
if ~exist(tmpdatadir,'dir'),
  mkdir(tmpdatadir);
end
outexpdirs = cell(size(idx_complete));
for i = 1:numel(idx_complete),
  experiment_name = experiment_names{idx_complete(i)};
  expdir = fullfile(rootdatadir,experiment_name);
  SymbolicCopyFlyBowlDataCaptureFiles(expdir,tmpdatadir);
  outexpdirs{i} = fullfile(tmpdatadir,experiment_name);
  cmd = sprintf('ln -s %s %s',fullfile(expdir,'movie.ufmf.ann'),fullfile(outexpdirs{i},'movie.ufmf.ann'));
  unix(cmd);
  cmd = sprintf('ln -s %s %s',fullfile(expdir,'ctrax_results.mat'),fullfile(outexpdirs{i},'ctrax_results.mat'));
  unix(cmd);
end
%%
isregmismatch = false(1,numel(idx_complete));
regnotes = cell(1,numel(idx_complete));
fnscheck = {'offX','offY','offTheta','scale',...
  'seconds_crop_start','seconds_crop_end','start_frame','end_frame'};
trajfns = {'x_mm','y_mm','a_mm','b_mm','theta_mm','firstframe','endframe','timestamps','dt'};
maxerrs = struct;
for j = 1:numel(trajfns),
  maxerrs.(trajfns{j}) = 0;
end
for j = 1:numel(fnscheck),
  maxerrs.(fnscheck{j}) = 0;
end
maxerrs = repmat(maxerrs,[1,numel(idx_complete)]);
rd_eps = struct('offX',2,'offY',2,'offTheta',.1,'scale',.01,'start_frame',0,'end_frame',0,...
  'seconds_crop_start',0,'seconds_crop_end',0);

order = randperm(numel(idx_complete));

for ii = 1:numel(idx_complete),
  
  i = order(ii);
  
  if ii > 1,
    iprev = order(ii-1);
    fprintf('%s: %d',experiment_names{idx_complete(iprev)},isregmismatch(iprev));
    if isregmismatch(iprev),
      fprintf(', %s',regnotes{iprev});
    end
    fprintf('\n');
  end
    
  try
    trx = FlyBowlRegisterTrx(outexpdirs{i},'analysis_protocol',analysis_protocol);
  catch ME
    warning('Registration failed: %s',getReport(ME));
    isregmismatch(i) = true;
    regnotes{i} = 'registration failed';
    continue;
  end

  expdir = fullfile(rootdatadir,experiment_names{idx_complete(i)});
  td = load(fullfile(expdir,'registered_trx.mat'));
  if numel(td.trx) ~= numel(trx),
    isregmismatch(i) = true;
    regnotes{i} = 'trx structures are different sizes';
    continue;
  end
  traj_eps = struct;
  for j = 1:numel(trajfns),
    maxerrs(i).(trajfns{j}) = 0;
    tmp = [trx.(trajfns{j})];
    traj_eps.(trajfns{j}) = (max(tmp)-min(tmp))*EPSILON;
  end
  for j = 1:numel(fnscheck),
    maxerrs(i).(fnscheck{j}) = 0;
  end
  for fly = 1:numel(trx),
    for j = 1:numel(trajfns),
      trajfn = trajfns{j};
      if numel(trx(fly).(trajfn)) ~= numel(td.trx(fly).(trajfn)),
        isregmismatch(i) = true;
        regnotes{i} = sprintf('%s sizes do not match for fly %d',trajfn,fly);
        break;
      end
      if strcmp(trajfn,'theta_mm'),
        errcurr = abs(modrange(trx(fly).(trajfn)-td.trx(fly).(trajfn),-pi,pi));
      else
        errcurr = abs(trx(fly).(trajfn)-td.trx(fly).(trajfn));
      end
      maxerrs(i).(trajfn) = max(maxerrs(i).(trajfn),max(errcurr));
    end
  end
  if isregmismatch(i),
    continue;
  end
  for j = 1:numel(trajfns),
    trajfn = trajfns{j};
    if maxerrs(i).(trajfn) > traj_eps.(trajfn),
      isregmismatch(i) = true;
      regnotes{i} = sprintf('max error between %s = %f > %f',trajfn,maxerrs(i).(trajfn),traj_eps.(trajfn));
      break;
    end
  end
  if isregmismatch(i),
    continue;
  end
  % check registration data
  tmp1 = load(fullfile(expdir,'registrationdata.mat'));
  tmp2 = load(fullfile(outexpdirs{i},'registrationdata.mat'));
  for j = 1:numel(fnscheck),
    fn = fnscheck{j};
    if strcmp(fn,'offTheta'),
      errcurr = abs(modrange(tmp1.(fn)-tmp2.(fn),-pi,pi));
    else
      errcurr = abs(tmp1.(fn)-tmp2.(fn));
    end
    maxerrs(i).(fn) = errcurr;
    if errcurr > rd_eps.(fn),
      isregmismatch(i) = true;
      regnotes{i} = sprintf('max error between %s = %f > %f',fn,errcurr,rd_eps.(fn));
      break;
    end
  end
end

%% copy over registraion, sex-classification results

for i = 1:numel(idx_complete),
  experiment_name = experiment_names{idx_complete(i)};
  expdir = fullfile(rootdatadir,experiment_name);
  outexpdir = fullfile(tmpdatadir,experiment_name);
  cmd = sprintf('rm %s',fullfile(outexpdir,'registered_trx.mat'));
  unix(cmd);
  cmd = sprintf('ln -s %s %s',fullfile(expdir,'registered_trx.mat'),fullfile(outexpdirs{i},'registered_trx.mat'));
  unix(cmd);
end

%% test out per-frame features

perframechecksuccess = false(1,numel(idx_complete));

tmpexperiment_names = experiment_names(idx_complete);

parfor i = 1:numel(idx_complete),
  try
    perframechecksuccess(i) = FlyBowlCheckPerFrameFeatures(...
      fullfile(rootdatadir,tmpexperiment_names{i}),...
      'analysis_protocol',analysis_protocol,...
      'outfilename',fullfile(tmpdatadir,['CompletePerFrameCheck_',tmpexperiment_names{i}]));
    fprintf('%s: %d\n',tmpexperiment_names{i},success);
  catch ME,
    warning(getReport(ME));
  end
end

%% copy over per-frame features

for i = 1:numel(idx_complete),
  experiment_name = experiment_names{idx_complete(i)};
  expdir = fullfile(rootdatadir,experiment_name);
  outexpdir = fullfile(tmpdatadir,experiment_name);
  if exist(fullfile(outexpdirs{i},'perframe'),'dir'),
    warning('perframe directory exists for %d %s, skipping',i,experiment_name);
    continue;
  end
  cmd = sprintf('ln -s %s %s',fullfile(expdir,'perframe'),fullfile(outexpdirs{i},'perframe'));
  unix(cmd);
end

%% check per-frame stats
stats_params = ReadStatsPerFrameFeatures2(fullfile('settings',analysis_protocol,'stats_perframefeatures.txt'));

fns = cell(1,numel(stats_params));
for i = 1:numel(stats_params),
  fns{i} = sprintf('%s_fly%s_frame%s',stats_params(i).field,stats_params(i).flycondition,stats_params(i).framecondition);
end
sagefns = unique({stats_params.field});

statsnotes = cell(1,numel(idx_complete));
isstatsmismatch = false(1,numel(idx_complete));
EPSILON = .01;
for i = 1:numel(idx_complete),
  
  if i > 1,
    iprev = i-1;
    fprintf('%d %s: %d',iprev,experiment_names{idx_complete(iprev)},isstatsmismatch(iprev));
    if isstatsmismatch(iprev),
      fprintf(', %s',statsnotes{iprev});
    end
    fprintf('\n');
  end
  
  %FlyBowlComputePerFrameStats2(outexpdirs{i},'analysis_protocol',analysis_protocol,'verbose',false);
  
  stats = load(fullfile(outexpdirs{i},'stats_perframe.mat'));
  %oldstats = load(fullfile(rootdatadir,tmpexperiment_names{i},'stats_perframe.mat'));
  
  
  
  
%   missingfns = setdiff(fns,fieldnames(oldstats.statsperexp));
%   extrafns = setdiff(fieldnames(oldstats.statsperexp),fns);
%   if ~isempty(missingfns),
%     statsnotes{i} = ['Missing stats:',sprintf(' %s',missingfns{:})];
%     isstatsmismatch(i) = true;
%     continue;
%   end
%   if ~isempty(extrafns),
%     statsnotes{i} = ['Extra stats:',sprintf(' %s',extrafns{:})];
%     isstatsmismatch(i) = true;
%     continue;
%   end
%   
%   for j = 1:numel(fns),
%     err = abs(stats.statsperexp.(fns{j}).meanmean-oldstats.statsperexp.(fns{j}).meanmean);
%     if err > EPSILON,
%       statsnotes{i} = sprintf('%s mismatch by %f',fns{j},err);
%       isstatsmismatch(i) = true;
%       break;
%     end
%   end
%   
%   if isstatsmismatch(i),
%     continue;
%   end
%   
%   for j = 1:numel(fns),
%     err = abs(stats.statsperexp.(fns{j}).stdmean-oldstats.statsperexp.(fns{j}).stdmean);
%     if err > EPSILON,
%       statsnotes{i} = sprintf('std of %s mismatch by %f',fns{j},err);
%       isstatsmismatch(i) = true;
%       break;
%     end
%   end
%   
%   if isstatsmismatch(i),
%     continue;
%   end

  % try loading from SAGE
  sagedata = SAGEGetBowlData('experiment_name',['FlyBowl_',experiment_names{idx_complete(i)}],...
    'checkflags',false,'removemissingdata',false);
  
  sagefns_curr = fieldnames(sagedata);
  sagefns_curr(cellfun(@isempty,regexp(sagefns_curr,'^stats_perframe','match','once'))) = [];
  sagefns_curr = regexprep(sagefns_curr,'^stats_perframe_','');
  missingfns = setdiff(sagefns,sagefns_curr);
  extrafns = setdiff(sagefns_curr,sagefns);
  
  if ~isempty(missingfns),
    statsnotes{i} = ['Missing SAGE stats:',sprintf(' %s',missingfns{:})];
    isstatsmismatch(i) = true;
    continue;
  end
  if ~isempty(extrafns),
    statsnotes{i} = ['Extra SAGE stats:',sprintf(' %s',extrafns{:})];
    isstatsmismatch(i) = true;
    continue;
  end
  
  for j = 1:numel(fns),
    v1 = stats.statsperexp.(fns{j}).meanmean;
    fn1 = ['stats_perframe_',stats_params(j).field];
    fn2 = ['fly',stats_params(j).flycondition,'_frame',stats_params(j).framecondition];
    if ~isfield(sagedata.(fn1).meanmean_perexp,fn2),
      statsnotes{i} = sprintf('Missing SAGE stat: %s.%s',fn1,fn2);
      isstatsmismatch(i) = true;
      break;
    end
      
    v2 = sagedata.(fn1).meanmean_perexp.(fn2);
      
    err = abs(v2-v1);
    if err > EPSILON,
      statsnotes{i} = sprintf('SAGE %s mismatch by %f',fns{j},err);
      isstatsmismatch(i) = true;
      break;
    end
  end
  
  if isstatsmismatch(i),
    continue;
  end
  
  for j = 1:numel(fns),
    v1 = stats.statsperexp.(fns{j}).stdmean;
    fn1 = ['stats_perframe_',stats_params(j).field];
    fn2 = ['fly',stats_params(j).flycondition,'_frame',stats_params(j).framecondition];
    if ~isfield(sagedata.(fn1).stdmean_perexp,fn2),
      statsnotes{i} = sprintf('Missing SAGE std stat: %s.%s',fn1,fn2);
      isstatsmismatch(i) = true;
      break;
    end
      
    v2 = sagedata.(fn1).stdmean_perexp.(fn2);
      
    err = abs(v2-v1);
    if err > EPSILON,
      statsnotes{i} = sprintf('SAGE std %s mismatch by %f',fns{j},err);
      isstatsmismatch(i) = true;
      break;
    end
  end
  
  
end
    

%% check for extra/missing fns in all directories

statsok = nan(1,numel(experiment_names));

parfor i = 1:numel(experiment_names),
  
  filename = fullfile(rootdatadir,experiment_names{i},'stats_perframe.mat');
  if ~exist(filename,'file'),
    continue;
  end
  oldstats = load(fullfile(rootdatadir,experiment_names{i},'stats_perframe.mat'),'statsperexp');
  missingfns = setdiff(fns,fieldnames(oldstats.statsperexp));
  extrafns = setdiff(fieldnames(oldstats.statsperexp),fns);
  statsok(i) = isempty(missingfns) && isempty(extrafns);
  fprintf('%d %s: %d\n',i,experiment_names{i},statsok(i));
end

save tmp.mat

%%

badexps = experiment_names(statsok==0);
tmp = SAGEListBowlExperiments('experiment_name',...
  cellfun(@(x) ['FlyBowl_',x],badexps,'UniformOutput',false),...
  'checkflags',false);
tmp2 = cellfun(@(x) x(9:end),{tmp.experiment_name},'UniformOutput',false);

for i = 1:numel(badexps),

%   if ~ismember(badexps{i},sus),
%     continue;
%   end
  
  fprintf('%s\t',badexps{i});
  if ismember(badexps{i},experiments_unarchive),
    fprintf('unarchive\t');
  elseif ismember(badexps{i},experiments_retrack),
    fprintf('retrack\t');
  elseif ismember(badexps{i},experiments_fix),
    fprintf('fix\t');
  else
    fprintf('unrequested\t');
  end
  
  j = find(strcmp(badexps{i},tmp2),1);
  if isempty(j),
    filename = fullfile(rootdatadir,badexps{i},'Metadata.xml');
    if ~exist(filename,'file'),
      fprintf('no metadata file\n');
    else
      try
        tmpmetadata = ReadMetadataFile(filename);
        fprintf('Not in SAGE\n');
      catch ME
        fprintf('Could not read metadata file\n');
      end
    end

  elseif strcmp(tmp(j).automated_pf,'F'),
    fprintf('automated_pf = F, automated_pf_category = %s\n',tmp(j).automated_pf_category);
  elseif strcmp(tmp(j).manual_pf,'F'),
    fprintf('manual_pf = F\n');
  elseif ~strcmp(tmp(j).screen_type,'primary'),
    fprintf('screen_type = %s\n',tmp(j).screen_type);
  else
    fprintf('no error\n');
  end
end

%% check for extra/missing perframe mat files in all directories

perframeok = true(1,numel(experiment_names));
perframefns = importdata('settings/20120706/perframefns.txt');

parfor i = 1:numel(experiment_names),
  
  perframedir = fullfile(rootdatadir,experiment_names{i},'perframe');
  for j = 1:numel(perframefns),
    if ~exist(fullfile(perframedir,[perframefns{j},'.mat']),'file'),
      perframeok(i) = false;
      break;
    end
  end
  fprintf('%d %s: %d\n',i,experiment_names{i},perframeok(i));
end

save tmp.mat

%%

for i = 1:numel(idx_complete),
    
  try
    FlyBowlPlotPerFrameStats2(outexpdirs{i},'analysis_protocol',analysis_protocol);
  catch ME,
    warning(getReport(ME));
  end
  
end
  

%% check if plotting worked

hist_params = ReadStatsPerFrameFeatures2(fullfile('settings',analysis_protocol,'hist_perframefeatures.txt'));
histfns = cell(size(hist_params));
for i = 1:numel(hist_params),
  if strcmp(hist_params(i).framecondition,'any'),
    histfns{i} = hist_params(i).field;
  else
    histfns{i} = [hist_params(i).field,'_',hist_params(i).framecondition];
  end
end
histfns = unique(histfns);

plotsok = true(1,numel(idx_complete));
for i = 1:numel(idx_complete),
    
  analysisdir = fullfile(outexpdirs{i},'analysis_plots');
  if ~exist(analysisdir,'dir'),
    fprintf('%s does not exist\n',analysisdir);
    plotsok(i) = false;
    continue;
  end

  if ~exist(fullfile(analysisdir,'stats.png'),'file'),
    fprintf('%s does not exist\n',fullfile(analysisdir,'stats.png'));
    plotsok(i) = false;
    continue;
  end
  
  mindatenum = inf;
  for j = 1:numel(histfns),
    
    filename = fullfile(analysisdir,['hist_',histfns{j},'.png']);
    tmp = dir(filename);
    if isempty(tmp),
      fprintf('%s does not exist\n',fullfile(analysisdir,['hist_',histfns{j},'.png']));
      plotsok(i) = false;
      break;
    end
    mindatenum = min(mindatenum,tmp.datenum);
    
  end

  if plotsok(i),
    fprintf('%s first timestamp: %s\n',experiment_names{idx_complete(i)},datestr(mindatenum,'yyyymmddTHHMMSS'));
  end
  
end

%% 
fid = fopen('expdirs_20120819.txt','w');
for i = 1:numel(idx_complete),
  fprintf(fid,'%s\n',fullfile(pwd,outexpdirs{i}));
end

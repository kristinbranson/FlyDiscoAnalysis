function [datamerge,experiment_ids,iswarning] = SAGEGetBowlData(varargin)

iswarning = false;
experiment_ids = [];
currdir = which('SAGEGetBowlData');
[currdir,~] = fileparts(currdir);

%% parse inputs
[docheckflags,daterange,SAGEpath,removemissingdata,dataset,rootdir,MAX_SET_TIMERANGE,MAX_EXPS_PER_SET,unflatten,...
  analysis_protocol,settingsdir,datalocparamsfilestr, CIRCLECENTERX,CIRCLECENTERY,readfromfile,leftovers] = ...
  myparse_nocheck(varargin,...
  'checkflags',true,...
  'daterange',[],...
  'SAGEpath',fullfile(currdir,'..','SAGE','MATLABInterface','Trunk'),...
  'removemissingdata',true,...
  'dataset','data',...
  'rootdir',0,...
  'MAX_SET_TIMERANGE',10/(24*60),...
  'MAX_EXPS_PER_SET',4,...
  'unflatten',true,...
  'analysis_protocol','current',...
  'settingsdir','/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings',...
  'datalocparamsfilestr','dataloc_params.txt',...
  'CIRCLECENTERX',1025/2,'CIRCLECENTERY',1025/2,...
  'readfromfile',false);

%% add SAGE to path
if ~exist('SAGE.Lab','class'),
  addpath(SAGEpath);
end

%% add input queries
allqueries = {};
if mod(numel(leftovers),2) ~= 0,
  error('Even number of parameters expected');
end
for i = 1:2:numel(leftovers)-1,
  fn = leftovers{i};
  val = leftovers{i+1};
  
  % query arguments are added raw
  if strcmpi(fn,'query'),
    if iscell(val),
      allqueries = [allqueries,val];
    else
      allqueries{end+1} = val;
    end
    continue;
  end

  if ischar(val) || numel(val) == 1,
    % single possible value
    allqueries{end+1} = SAGE.Query.Compare(fn,'=',val);
  elseif isempty(val),
    % empty set means any value allowed -- ignore
    continue;
  elseif iscell(val),
    % set of possible values
    queriescurr = cell(size(val));
    for j = 1:numel(val),
      queriescurr{j} = SAGE.Query.Compare(fn,'=',val{j});
    end
    allqueries{end+1} = SAGE.Query.Any(queriescurr{:}); %#ok<*AGROW>
  elseif isnumeric(val) && numel(val) == 2,
    % range of possible values
    if ~isnan(val(1)),
      % lower bound
      allqueries{end+1} = SAGE.Query.Compare(fn,'>',num2str(val(1)));
    end
    if ~isnan(val(2)),
      % upper bound
      allqueries{end+1} = SAGE.Query.Compare(fn,'<',num2str(val(2)));
    end
  else % not handled
    error('Not implemented: value for field %s must be a cell array of possible values, a single possible value, or an upper and lower bound.',fn);
  end
end

%% add queries for checking flag
if docheckflags,
  allqueries{end+1} = SAGE.Query.Compare('flag_aborted','=','0');
  allqueries{end+1} = SAGE.Query.Compare('automated_pf','=','P');
  allqueries{end+1} = SAGE.Query.Compare('flag_redo','=','0');
  % TODO: this should be = 'P' someday
  allqueries{end+1} = SAGE.Query.Any(SAGE.Query.Compare('manual_pf','=','P'),...
    SAGE.Query.Compare('manual_pf','=','U'));
end

%% add date range queries
if ~isempty(daterange),
  if ischar(daterange),
    daterange = {daterange};
  end
  if ~isempty(daterange{1}),
    allqueries{end+1} = SAGE.Query.Compare('exp_datetime','>',daterange{1});
  end
  if numel(daterange) > 1 && ~isempty(daterange{2}),
    allqueries{end+1} = SAGE.Query.Compare('exp_datetime','<',daterange{2});
  end
end

%% grab data: put everything into as few queries as possible for speed

tmp = SAGE.Lab('olympiad');
tmp = tmp.assay('bowl');
bowlAssay_data = tmp.dataSet(dataset);

if isempty(allqueries),
  data = bowlAssay_data.findData();
else
  query = SAGE.Query.All(allqueries{:});
  % get data
  data = bowlAssay_data.findData(query);
end

% convert some fields to numbers
data = convert2numeric(data);

%% put the data together

if isempty(data),
  datamerge = data;
  return;
end

[experiment_ids,uniqueidx,idxperrow] = unique([data.experiment_id]);
datamerge = rmfield(data(uniqueidx),intersect(fieldnames(data),{'data','data_type'}));
if isfield(data,'data_type'),
  for i = 1:numel(data),
    j = idxperrow(i);
    fn = data(i).data_type;
    % need to shorten some names
    fn = regexprep(fn,'^hist_perframe','hist');
    fn = regexprep(fn,'^stats_perframe','stats');
    datamerge(j).(fn) = data(i).data;
  end
end

% convert file system path
if ischar(rootdir),
  isfilesystempath = isfield(datamerge,'file_system_path');
  for i = 1:numel(datamerge),
    if isfilesystempath && ~strcmpi(datamerge(i).file_system_path,'NULL'),
      [~,basename] = fileparts(datamerge(i).file_system_path);
    else
      basename = regexprep(datamerge(i).experiment_name,'^FlyBowl_','');
    end
    datamerge(i).file_system_path = fullfile(rootdir,basename);
  end
end

%% read missing data from file system

if readfromfile,

  filestrs = {
    'bias_diagnostics.mat'
    'bkgd_diagnostics.mat'
    'registrationdata.mat'
    'video_diagnostics.mat'
    'ctrax_diagnostics.txt'
    'sexclassifier_diagnostics.txt'
    'temperature_diagnostics.txt'
    };
  prefixes = {
    'bias_diagnostics_'
    'bkgd_diagnostics_'
    'registrationdata_'
    'video_diagnostics_'
    'ctrax_diagnostics_'
    'sexclassifier_diagnostics_'
    'temperature_diagnostics_'
    };

  isread = zeros(1,numel(filestrs));
  fnsread = {};
  fileidx = [];
  for i = 1:numel(datamerge),
    
    for j = find(~isread),
      filename = fullfile(datamerge(i).file_system_path,filestrs{j});
      if exist(filename,'file'),
        [~,~,ext] = fileparts(filestrs{j});
        if strcmp(ext,'.mat'),
          s = load(filename);
        else
          s = ReadParams(filename);
        end
        fprintf('Read fields for file %s from experiment %s\n',filestrs{j},datamerge(i).experiment_name);
        fnscurr = fieldnames(s)';
        fnscurr = cellfun(@(x) [prefixes{j},x],fnscurr,'UniformOutput',false);
        fnsread = [fnsread,fnscurr];
        fileidx(end+1:end+numel(fnscurr)) = j;
        isread(j) = i;
      end
      if all(isread),
        break;
      end
    end    
  end

  % choose those that match data_type
  i = find(strcmp(varargin(1:2:end-1),'data_type'));
  if ~isempty(i),
    i = i*2;
    data_types_requested = varargin{i};
    idxkeep = false(size(fileidx));
    for i = 1:numel(data_types_requested),
      s = data_types_requested{i};
      s = ['^',strrep(s,'*','.*'),'$'];
      idxkeep = idxkeep | ~cellfun(@isempty,regexp(fnsread,s,'once'));
    end    
  else
    idxkeep = true(size(fileidx));
  end
  
  for i = find(idxkeep),
    if ~isfield(datamerge,fnsread{i}),
      [datamerge.(fnsread{i})] = deal([]);
    end
  end
  
  fileidxread = unique(fileidx(idxkeep));
  fns = fieldnames(datamerge);
  fnsmetadata = fieldnames(data);
  ignoreidx = ismember(fns,fnsmetadata);
  for i = 1:numel(datamerge),
    ismissing = ~ignoreidx & structfun(@isempty,datamerge(i));
    if ~any(ismissing),
      continue;
    end
    for j = fileidxread,
      filename = fullfile(datamerge(i).file_system_path,filestrs{j});
      if ~exist(filename,'file'),
        continue;
      end
      [~,~,ext] = fileparts(filestrs{j});
      if strcmp(ext,'.mat'),
        s = load(filename);
      else
        s = ReadParams(filename);
      end
      fnscurr = fieldnames(s)';
      fnscurr1 = cellfun(@(x) [prefixes{j},x],fnscurr,'UniformOutput',false);
      idxcurr = find(ismember(fnscurr1,fns));
      for k = idxcurr,
        fn = fnscurr{k};
        fn1 = fnscurr1{k};
        if ~isempty(datamerge(i).(fn1)),
          fprintf('%s, %s: replacing %s with %s\n',datamerge(i).experiment_name,fn1,mat2str(datamerge(i).(fn1)),mat2str(s.(fn)));
        else
          if (~ischar(s.(fn)) && numel(s.(fn)) > 100) || ndims(s.(fn)) > 2 || iscell(s.(fn)) || isstruct(s.(fn)),
            fprintf('%s, %s: read from file, size = %s\n',datamerge(i).experiment_name,fn1,mat2str(size(s.(fn))));
          else
            fprintf('%s, %s: read %s from file\n',datamerge(i).experiment_name,fn1,mat2str(s.(fn)));
          end
        end
        datamerge(i).(fn1) = s.(fn);
      end
    end
    
  end
  
  
end

%% add rig x bowl
if isfield(datamerge,'rig') && isfield(datamerge,'bowl'),
  for i = 1:numel(datamerge),
    datamerge(i).rig_bowl = sprintf('%d%s',datamerge(i).rig,datamerge(i).bowl);
  end
end

%% add plate x bowl
if isfield(datamerge,'plate') && isfield(datamerge,'bowl'),
  for i = 1:numel(datamerge),
    datamerge(i).plate_bowl = sprintf('%d%s',datamerge(i).plate,datamerge(i).bowl);
  end
end


%% add line__effector
if isfield(datamerge,'line_name') && isfield(datamerge,'effector'),
  for i = 1:numel(datamerge),
    datamerge(i).line__effector = sprintf('%s__%s',datamerge(i).line_name,datamerge(i).effector);
  end
end

%% add ctrax_diagnostics_nframes_not_tracked
if isfield(datamerge,'ctrax_diagnostics_nframes_analyzed') && ...
    isfield(datamerge,'ufmf_diagnostics_summary_nFrames'),
  for i = 1:numel(datamerge),
    datamerge(i).ctrax_diagnostics_nframes_not_tracked = ...
      datamerge(i).ufmf_diagnostics_summary_nFrames - ...
      datamerge(i).ctrax_diagnostics_nframes_analyzed;
  end
end

%% add ctrax_diagnostics_mean_nsplit

if isfield(datamerge,'ctrax_diagnostics_sum_nsplit') && ...
    isfield(datamerge,'ctrax_diagnostics_nlarge_split'),
  for i = 1:numel(datamerge),
    datamerge(i).ctrax_diagnostics_mean_nsplit = ...
      datamerge(i).ctrax_diagnostics_sum_nsplit ./ ...
      datamerge(i).ctrax_diagnostics_nlarge_split;
  end
end

%% add reg_pxpermm
if isfield(datamerge,'registrationdata_scale'),
  for i = 1:numel(datamerge),
    datamerge(i).registrationdata_pxpermm = 1./datamerge(i).registrationdata_scale;
  end
end

%% add "set" -- super-experiment

[line_names,~,lineidx] = unique({datamerge.line_name});
sets = nan(1,numel(datamerge));
seti = 0;
for linei = 1:numel(line_names),
  expidx1 = find(lineidx==linei);
  [rigs,~,rigidx] = unique([datamerge(expidx1).rig]);
  for rigi = 1:numel(rigs),
    expidx2 = expidx1(rigidx==rigi);

    % sort by datetime
    [exp_datenum,order] = sort(datenum({datamerge(expidx2).exp_datetime},'yyyymmddTHHMMSS'));
    expidx2 = expidx2(order);
    min_set_time = 0;
    nperset = 0;
  
    for i = 1:numel(expidx2),
      
      % start a new set?
      if exp_datenum(i)-min_set_time > MAX_SET_TIMERANGE || nperset >= MAX_EXPS_PER_SET,
        min_set_time = inf;
        seti = seti+1;
        nperset = 0;
      end
      sets(expidx2(i)) = seti;
      nperset = nperset+1;
      min_set_time = min(min_set_time,exp_datenum(i));
      
    end
  end
end

for seti = 1:max(sets),
  expidx = find(sets==seti);
  %fprintf('Experiments in set %d:\n',seti);
  %fprintf('%s\n',datamerge(expidx).experiment_name);
  [~,order] = sort({datamerge(expidx).exp_datetime});
  min_datetime = datamerge(expidx(order(1))).exp_datetime;
  set_name = sprintf('%s__Rig%d__%s',datamerge(expidx(1)).line_name,...
    datamerge(expidx(1)).rig,...
    min_datetime);
  for i = expidx(:)',
    datamerge(i).set = set_name;
  end
end

%% add in normalized bkgd stats

fns_bkgd = {'bkgd_diagnostics_mean_bkgdcenter','bkgd_diagnostics_mean_bkgdcenter_llr',...
  'registrationdata_circleCenterX','registrationdata_circleCenterY','registrationdata_bowlMarkerTheta',...
  'temperature_diagnostics_mean','temperature_diagnostics_max'};

if any(isfield(datamerge,fns_bkgd)) && isfield(datamerge,'plate_bowl'),
  
  % read bkgd stats from file
  datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
  dataloc_params = ReadParams(datalocparamsfile);
  platebowlstatsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.statsperplatebowlfilestr);
  try
    platebowlstats = ReadParams(platebowlstatsfile);
  catch ME,
    getReport(ME)
    % if reading from file fails, just set to be empty
    platebowlstats = struct;
    platebowlstats.platebowls = {};
    for i = 1:numel(fns_bkgd),
      platebowlstats.(fns_bkgd{i}) = [];
    end
  end
  % find platebowls not in the file
  [ism,platebowlidx] = ismember({datamerge.plate_bowl},platebowlstats.platebowls);
  if any(~ism),
    % compute stats from this data if no stat in file
    new_plate_bowls = setdiff({datamerge.plate_bowl},platebowlstats.platebowls);
    n0 = numel(platebowlstats.platebowls);
    for i = 1:numel(fns_bkgd),
      fn = fns_bkgd{i};
      if ~isfield(datamerge,fn);
        continue;
      end
      y_bkgd = nan(1,numel(datamerge));
      goodidx = ~cellfun(@isempty,{datamerge.(fn)});
      y_bkgd(goodidx) = [datamerge(goodidx).(fn)];
      for j = 1:numel(new_plate_bowls),
        idx = strcmp({datamerge.plate_bowl},new_plate_bowls{j}) & goodidx;
        platebowlstats.platebowls{n0+j} = new_plate_bowls{j};
        platebowlidx(idx) = n0+j;
        platebowlstats.(fn)(n0+j) = median(y_bkgd(idx));
      end
    end
  end
  for i = 1:numel(fns_bkgd),
    fn = fns_bkgd{i};
    if isfield(datamerge,fn) && isfield(platebowlstats,fn),
      fnz = [fn,'_norm_platebowl'];
      for j = 1:numel(datamerge),
        datamerge(j).(fnz) = datamerge(j).(fn) - platebowlstats.(fn)(platebowlidx(j));
      end
    end
  end
end

if all(isfield(datamerge,{'registrationdata_circleCenterX','registrationdata_circleCenterY'})),
  for i = 1:numel(datamerge),
    datamerge(i).registrationdata_circleCenterMaxXY = ...   
      max(abs(datamerge(i).registrationdata_circleCenterX - CIRCLECENTERX),...
      abs(datamerge(i).registrationdata_circleCenterY - CIRCLECENTERY));
  end
end

%% format stats, hist into substructs

if unflatten,
  
  statfns = {'nflies_analyzed'
    'flies_analyzed'
    'mean_perfly'
    'std_perfly'
    'Z_perfly'
    'fracframesanalyzed_perfly'
    'prctiles_perfly'
    'meanmean_perexp'
    'stdmean_perexp'
    'meanstd_perexp'
    'Z_perexp'
    'meanZ_perexp'
    'meanprctiles_perexp'
    'stdprctiles_perexp'};
  [datamerge,iswarning1] = format2substructs(datamerge,'stats',statfns,'stats_perframe',removemissingdata);
  iswarning = iswarning || iswarning1;
  statfns = {'nflies_analyzed'
    'flies_analyzed'
    'frac_linear_perfly'
    'frac_log_perfly'
    'Z_perfly'
    'fracframesanalyzed_perfly'
    'mean_frac_linear_perexp'
    'std_frac_linear_perexp'
    'mean_frac_log_perexp'
    'std_frac_log_perexp'
    'Z_perexp'
    'meanZ_perexp'};
  [datamerge,iswarning1] = format2substructs(datamerge,'hist',statfns,'hist_perframe',removemissingdata);
  iswarning = iswarning || iswarning1;
  
end

function in = convert2numeric(in)

numeric_fields = {...
  'experiment_id'...
  'flag_aborted'...
  'hours_sorted'...
  'hours_starved'...
  'humidity'...
  'line_id'...
  'num_flies'...
  'num_flies_dead'...
  'plate'...
  'rearing_incubator'...
  'rig'...
  'seconds_fliesloaded'...
  'seconds_shiftflytemp'...
  'session_id'...
  'temperature'...
  'top_plate'...
  'flag_redo'...
  'flag_review'...
  };

% flags were once strings but are now numeric
flag_fields = {...
  'flag_redo'...
  'flag_review'...
};

for j = 1:numel(numeric_fields),
  fn = numeric_fields{j};
  if isfield(in,fn),
    for i = 1:numel(in),
      if ischar(in(i).(fn)),
        v = str2double(in(i).(fn));
        if isnan(v) && ismember(fn,flag_fields),
          v = ~(isempty(in(i).(fn)) || strcmpi(in(i).(fn),'None'));
        end          
        in(i).(fn) = v;
      end
    end
  end
end

function [datamerge,iswarning] = format2substructs(datamerge,prefix,statfns,newprefix,removemissingdata)

iswarning = false;

% get all perframe fns
expr = ['^',prefix,'_(?<perframefn>.+)_(?<statfn>',...
  statfns{1},sprintf('|%s',statfns{2:end}),')$'];
fns = fieldnames(datamerge);
m = regexp(fns,expr,'names','once');
statidx = ~cellfun(@isempty,m);
m(~statidx) = [];
fns(~statidx) = [];
m = cell2mat(m);
if isempty(m), return; end
perframefns = {m.perframefn};
[perframefns,~,idx] = unique(perframefns);
statfns = {m.statfn};

allmissingdata = false(size(datamerge));
for i = 1:numel(perframefns),
  perframefn = perframefns{i};
  idxcurr = find(idx == i);
  conditionfn = sprintf('%s_%s_conditions',prefix,perframefn);
  newfn = sprintf('%s_%s',newprefix,perframefn);
  nfliesanalyzedfn = sprintf('%s_%s_nflies_analyzed',prefix,perframefn);
  if ~isfield(datamerge,conditionfn),
    fprintf('Condition data %s does not exist, skipping reformatting for %s\n',conditionfn,perframefn);
    continue;
  end
  missingdata = false(size(datamerge));
  allconditions = {};
  for j = 1:numel(datamerge),
    try
    datamerge(j).(newfn) = struct;
    conditions = datamerge(j).(conditionfn);
    nconditions = numel(conditions);
    nfliesanalyzed = datamerge(j).(nfliesanalyzedfn);
    nfliesanalyzedtotal = sum(nfliesanalyzed);
    for kk = 1:numel(idxcurr),
      k = idxcurr(kk);
      statfn = statfns{k};
      fn = sprintf('%s_%s_%s',prefix,perframefn,statfn);
      datamerge(j).(newfn).(statfn) = struct;
      if nconditions == 0,
        if isempty(datamerge(j).(fn)),
          missingdata(j) = true;
        else
          error('For %s, nconditions = 0 but non-empty data',fn); 
        end
      else
        allconditions = union(allconditions,conditions);
        if strcmp(statfn,'flies_analyzed') || ...
            (~isempty(regexp(statfn,'perfly$','once')) && ...
            ~(strcmp(prefix,'hist') && ismember(statfn,{'Z_perfly','fracframesanalyzed_perfly'}))),
          n = numel(datamerge(j).(fn))/nfliesanalyzedtotal;
          if n ~= round(n),
            error('Error splitting %s into different conditions: size of array not a multiple of nflies = %d',fn,nfliesanalyzedtotal); 
          end
          datamerge(j).(fn) = mat2cell(reshape(datamerge(j).(fn),[n,nfliesanalyzedtotal]),n,nfliesanalyzed);
          for l = 1:nconditions,
            datamerge(j).(newfn).(statfn).(conditions{l}) = datamerge(j).(fn){l};
          end
        else
          n = numel(datamerge(j).(fn)) / nconditions;
          if n ~= round(n),
            error('Error splitting %s into different conditions: size of array not a multiple of nconditions = %d',fn,nconditions); 
          end
          datamerge(j).(fn) = reshape(datamerge(j).(fn),[n,nconditions]);
          for l = 1:nconditions,
            datamerge(j).(newfn).(statfn).(conditions{l}) = datamerge(j).(fn)(:,l)';
          end
        end
      end
    end
    catch ME,
      warning('Error splitting stats for experiment %s:\n%s',datamerge(j).experiment_name,getReport(ME));
      missingdata(j) = true;
      iswarning = true;
    end
  end
  % create empty structs for missing data
  if removemissingdata,
    allmissingdata = allmissingdata | missingdata;
  else
    conditions = allconditions;
    for j = find(missingdata),
      for kk = 1:numel(idxcurr),
        k = idxcurr(kk);
        statfn = statfns{k};
        warning('Missing data for %s, %s for experiment %s',newfn,statfn,datamerge(j).experiment_name);
        iswarning = true;
        for l = 1:nconditions,
          datamerge(j).(newfn).(statfn).(conditions{l}) = [];
        end
      end
    end
  end
  datamerge = rmfield(datamerge,[{conditionfn};fns(idxcurr)]);
end

if removemissingdata,
  if any(allmissingdata),
    fprintf('Removing the following experiments, which were missing some data:\n');
    fprintf('%s\n',datamerge(allmissingdata).experiment_name);
    datamerge(allmissingdata) = [];
  end
end
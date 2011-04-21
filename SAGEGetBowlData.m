function [datamerge,experiment_ids] = SAGEGetBowlData(varargin)

currdir = which('SAGEGetBowlData');
[currdir,~] = fileparts(currdir);

%% parse inputs
[docheckflags,daterange,SAGEpath,removemissingdata,dataset,rootdir,leftovers] = ...
  myparse_nocheck(varargin,...
  'checkflags',true,...
  'daterange',[],...
  'SAGEpath',fullfile(currdir,'..','SAGE','MATLABInterface','Trunk'),...
  'removemissingdata',true,...
  'dataset','data',...
  'rootdir',0);

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
datamerge = rmfield(data(uniqueidx),{'data','data_type'});
for i = 1:numel(data),
  j = idxperrow(i);
  fn = data(i).data_type;
  % need to shorten some names
  fn = regexprep(fn,'^hist_perframe','hist');
  fn = regexprep(fn,'^stats_perframe','stats');
  datamerge(j).(fn) = data(i).data;
end

% convert file system path
if ischar(rootdir),
  isfilesystempath = isfield(datamerge,'file_system_path');
  for i = 1:numel(datamerge),
    if isfilesystempath,
      [~,basename] = fileparts(datamerge(i).file_system_path);
    else
      basename = regexprep(datamerge(i).experiment_name,'^FlyBowl_','');
    end
    datamerge(i).file_system_path = fullfile(rootdir,basename);
  end
end


%% format stats, hist into substructs

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
datamerge = format2substructs(datamerge,'stats',statfns,'stats_perframe',removemissingdata);
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
datamerge = format2substructs(datamerge,'hist',statfns,'hist_perframe',removemissingdata);

function in = convert2numeric(in)

numeric_fields = {...
  'experiment_id'...
  'environmental_chamber'...
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

function datamerge = format2substructs(datamerge,prefix,statfns,newprefix,removemissingdata)

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
          error('For %s, nconditions = 0 but non-empty data',fn); %#ok<SPERR>
        end
      else
        allconditions = union(allconditions,conditions);
        if strcmp(statfn,'flies_analyzed') || ...
            (~isempty(regexp(statfn,'perfly$','once')) && ...
            ~(strcmp(prefix,'hist') && ismember(statfn,{'Z_perfly','fracframesanalyzed_perfly'}))),
          n = numel(datamerge(j).(fn))/nfliesanalyzedtotal;
          if n ~= round(n),
            error('Error splitting %s into different conditions: size of array not a multiple of nflies = %d',fn,nfliesanalyzedtotal); %#ok<SPERR>
          end
          datamerge(j).(fn) = mat2cell(reshape(datamerge(j).(fn),[n,nfliesanalyzedtotal]),n,nfliesanalyzed);
          for l = 1:nconditions,
            datamerge(j).(newfn).(statfn).(conditions{l}) = datamerge(j).(fn){l};
          end
        else
          n = numel(datamerge(j).(fn)) / nconditions;
          if n ~= round(n),
            error('Error splitting %s into different conditions: size of array not a multiple of nconditions = %d',fn,nconditions); %#ok<SPERR>
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
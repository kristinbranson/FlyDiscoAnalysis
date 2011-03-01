function [datamerge,experiment_ids] = SAGEGetBowlData(varargin)

currdir = which('SAGEGetBowlData');
[currdir,~] = fileparts(currdir);

%% parse inputs
[docheckflags,daterange,SAGEpath,leftovers] = ...
  myparse_nocheck(varargin,...
  'checkflags',true,...
  'daterange',[],...
  'SAGEpath',fullfile(currdir,'..','SAGE','MATLABInterface','Trunk'));

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
  allqueries{end+1} = SAGE.Query.Compare('flag_redo','=','None');
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

bowlAssay_data = SAGE.Lab('olympiad').assay('bowl').dataSet('data');

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

[experiment_ids,uniqueidx,idxperrow] = unique([data.experiment_id]);
datamerge = rmfield(data(uniqueidx),{'data','data_type'});
for i = 1:numel(data),
  j = idxperrow(i);
  % need to shorten some names
  fn = data(i).data_type;
  fn = regexprep(fn,'^hist_perframe','hist');
  fn = regexprep(fn,'^stats_perframe','stats');
  datamerge(j).(fn) = data(i).data;
end

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
  };

for j = 1:numel(numeric_fields),
  fn = numeric_fields{j};
  if isfield(in,fn),
    for i = 1:numel(in),
      if ischar(in(i).(fn)),
        in(i).(fn) = str2double(in(i).(fn));
      end
    end
  end
end

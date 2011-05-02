function data = pullBowlData(qObjs,varargin)
%pullBowlData Pull Bowl Data from SAGE
%   pullBowlData is a SAGE data-pulling function for use with
%   OlyDat.DataSelector. qObjs is a cell array of query objects generated
%   by the DataSelector. data is a single data struct that results from
%   pulling data from the flybowl dataSets, cleaning and merging that data, and
%   computing additional statistics of interest.
%
%   See also OlyDat.DataSelector.

% FORMING QUERIES USING THE DATASELECTOR FOR THE BOWL:
%

% STAGES OF DATA PULLING:
% * Stage 1: actually pull data from env and AR. group/clean/merge. Put
% into standard format.
% * Stage 2: compute addn'l statistics. Right now this uses the Analyzer
% and i) computes env stats and ii) reduces multi-dim bstats. QC can take
% place here as well-- practically, the only QC check that actually goes
% off at the moment is for missing vibrational peaks. This stage isn't
% really sorted out yet. Does the QC really belong here?

assert(all(cellfun(@(x)isa(x,'SAGE.Query.Clause'),qObjs)));

% allow dataset to be specified
[dataset,rootdir,removemissingdata,ignore_missingdata_fns] = myparse(varargin,...
  'dataset','data',...
  'rootdir',0,...
  'removemissingdata',true,...
  'ignore_missingdata_fns',{'temperature_stream'});

%% grab data
dDS = SAGE.Lab('olympiad').assay('bowl').dataSet(dataset);

% hardcoded behavioral stats of interest. don't put anything, so we should
% grab everytrhing for now
bStatQry = SAGE.Query.Compare('experiment_name','=','FlyBowl_*');%,...

if isempty(bStatQry),
  dQueryObjs = qObjs(:);
else
  dQueryObjs = [qObjs(:);{bStatQry}];
end
dData = zlclApplyQueryObjsToDataSet(dQueryObjs,dDS);
if isempty(dData)
    warning('pullBowlData:noDData','No matching records in %s table. Aborting.',dataset);
    data = [];
    return;
end

%% group, clean, merge

dData = datagroupclean(dData,rootdir,removemissingdata,ignore_missingdata_fns);
% no data to merge
data = dData;

%% compute addn'l statistics

% add in exp_datenum
for i = 1:numel(data),
  data(i).edt = datenum(data(i).exp_datetime,'yyyymmddTHHMMSS');
end

end

function data = zlclApplyQueryObjsToDataSet(qObjs,dataSet)
flds = dataSet.fields;
fieldNames = {flds.name}';

% pick out parts of query that are relevant to this dataSet
relevantQObjs = cell(0,1);
for c = 1:numel(qObjs)
    fns = zlclGetAllFieldsInSAGEQuery(qObjs{c});
    tf = ismember(fns,fieldNames);
    if all(tf)
        % all fieldnames in this query object are in this dataSet
        relevantQObjs{end+1,1} = qObjs{c}; %#ok<AGROW>
    elseif ~any(tf)
        % no fieldnames in this query object are in this dataSet
    else
        % mixed situation
        warning('pullBowlData:zlclApplyQueryObjsToDataSet',...
            'Query clause with a mix of fields that are and are not part of dataSet ''%s''. Ignoring this clause.',dataSet.displayName);
    end
end

% pull data
qry = SAGE.Query.All(relevantQObjs{:});
data = dataSet.findData(qry);

if isempty(data)
    % Create a 0x1 struct with the right fields
    
    % when you pull data from SAGE, it automatically removes these three
    % fields
    reducedFieldNames = setdiff(fieldNames,{'data_columns';'data_format';'data_rows'}); 
    
    data = cell2struct(cell(0,numel(reducedFieldNames)),reducedFieldNames,2);
end
end

function fieldNames = zlclGetAllFieldsInSAGEQuery(q)
assert(isa(q,'SAGE.Query.Clause')&&isscalar(q));
switch q.type
    case 'Compare'
        fieldNames = {q.fieldName};
    case {'Any' 'All'}
        v = q.value;
        fieldNames = cell(0,1);
        for c = 1:numel(v)
            fieldNames = [fieldNames;zlclGetAllFieldsInSAGEQuery(v{c})]; %#ok<AGROW>
        end
    otherwise
        assert(false,'Unhandled query type.');
end
fieldNames = unique(fieldNames);
end

function datamerge = datagroupclean(data,rootdir,removemissingdata,ignore_missingdata_fns)

data = convert2numeric(data);
if isempty(data),
  datamerge = data;
  return;
end

% combine rows corresponding to the same experiment into one row
[~,uniqueidx,idxperrow] = unique([data.experiment_id]');
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

% format stats, hist into substructs
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
datamerge = format2substructs(datamerge,'stats',statfns,'stats',removemissingdata);
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
datamerge = format2substructs(datamerge,'hist',statfns,'hist',removemissingdata);

% check for missing data
ismissingdata = false(size(datamerge));
if removemissingdata,
  fns = setdiff(fieldnames(datamerge),ignore_missingdata_fns);
  for i = 1:numel(fns),
    fn = fns{i};
    ismissingdatacurr = cellfun(@isempty,{datamerge.(fn)});
    if ~all(ismissingdatacurr),
      ismissingdata = ismissingdata | ismissingdatacurr;
    end
    if any(ismissingdatacurr),
      warning(['The following experiments are missing data for %s:',sprintf('\n%s',datamerge(ismissingdatacurr).experiment_name)],fn);
    end
  end
  datamerge(ismissingdata) = [];  
end


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
  nfliesanalyzedfn = sprintf('%s_%s_nflies_analyzed',prefix,perframefn);
  if ~isfield(datamerge,conditionfn),
    fprintf('Condition data %s does not exist, skipping reformatting for %s\n',conditionfn,perframefn);
    continue;
  end
  missingdata = false(size(datamerge));
  allconditions = {};
  for j = 1:numel(datamerge),
    try
    conditions = datamerge(j).(conditionfn);
    nconditions = numel(conditions);
    nfliesanalyzed = datamerge(j).(nfliesanalyzedfn);
    nfliesanalyzedtotal = sum(nfliesanalyzed);
    for kk = 1:numel(idxcurr),
      k = idxcurr(kk);
      statfn = statfns{k};
      fn = sprintf('%s_%s_%s',prefix,perframefn,statfn);
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
            conditionfn1 = strrep(conditions{l},'flyany_frame','');
            newfn = sprintf('%s_%s_%s_%s',newprefix,perframefn,statfn,conditionfn1);
            datamerge(j).(newfn) = datamerge(j).(fn){l};
          end
        else
          n = numel(datamerge(j).(fn)) / nconditions;
          if n ~= round(n),
            error('Error splitting %s into different conditions: size of array not a multiple of nconditions = %d',fn,nconditions); %#ok<SPERR>
          end
          datamerge(j).(fn) = reshape(datamerge(j).(fn),[n,nconditions]);
          for l = 1:nconditions,
            conditionfn1 = strrep(conditions{l},'flyany_frame','');
            newfn = sprintf('%s_%s_%s_%s',newprefix,perframefn,statfn,conditionfn1);
            datamerge(j).(newfn) = datamerge(j).(fn)(:,l)';
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
          conditionfn1 = strrep(conditions{l},'flyany_frame','');
          newfn = sprintf('%s_%s_%s_%s',newprefix,perframefn,statfn,conditionfn1);
          datamerge(j).(newfn) = [];
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
end
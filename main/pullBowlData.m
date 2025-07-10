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
[dataset,rootdir,removemissingdata,ignore_missingdata_fns,...
  MAX_SET_TIMERANGE,MAX_EXPS_PER_SET,...
  analysis_protocol,settingsdir,datalocparamsfilestr,...
  CIRCLECENTERX,CIRCLECENTERY] = myparse(varargin,...
  'dataset','data',...
  'rootdir',0,...
  'removemissingdata',true,...
  'ignore_missingdata_fns',{'temperature_stream'},...
  'MAX_SET_TIMERANGE',10/(24*60),...
  'MAX_EXPS_PER_SET',4,...
  'analysis_protocol','current',...
  'settingsdir','/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings',...
  'datalocparamsfilestr','dataloc_params.txt',...
  'CIRCLECENTERX',1025/2,'CIRCLECENTERY',1025/2);

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
  
  % numeric date
  data(i).edt = datenum(data(i).exp_datetime,'yyyymmddTHHMMSS');
  
  % rig x bowl
  if isfield(data,'rig') && isfield(data,'bowl'),
    data(i).rig_bowl = sprintf('%d%s',data(i).rig,data(i).bowl);
  end

  % plate x bowl
  if isfield(data,'plate') && isfield(data,'bowl'),
    data(i).plate_bowl = sprintf('%d%s',data(i).plate,data(i).bowl);
  end
  
  % numerical date, rounded to day
  data(i).edt_day = floor(data(i).edt);
  
  % numerical date, rounded to week
  data(i).edt_week = floor((data(i).edt-2)/7)*7+2;
  
  % numerical date, rounded to month
  data(i).edt_month = datenum(datestr(data(i).edt,'yyyymm'),'yyyymm');
  
  % part of day
  hrofday = data(i).edt - floor(data(i).edt);
  data(i).part_of_day = double(hrofday>=.5);
  
  % hour of day
  data(i).hour_of_day = floor(hrofday*24);
  
end


%% add line__effector
if isfield(data,'line_name') && isfield(data,'effector'),
  for i = 1:numel(data),
    data(i).line__effector = sprintf('%s__%s',data(i).line_name,data(i).effector);
  end
end

%% add ctrax_diagnostics_nframes_not_tracked
if isfield(data,'ctrax_diagnostics_nframes_analyzed') && ...
    isfield(data,'ufmf_diagnostics_summary_nFrames'),
  for i = 1:numel(data),
    data(i).ctrax_diagnostics_nframes_not_tracked = ...
      data(i).ufmf_diagnostics_summary_nFrames - ...
      data(i).ctrax_diagnostics_nframes_analyzed;
  end
end

%% add ctrax_diagnostics_mean_nsplit

if isfield(data,'ctrax_diagnostics_sum_nsplit') && ...
    isfield(data,'ctrax_diagnostics_nlarge_split'),
  for i = 1:numel(data),
    data(i).ctrax_diagnostics_mean_nsplit = ...
      data(i).ctrax_diagnostics_sum_nsplit ./ ...
      data(i).ctrax_diagnostics_nlarge_split;
  end
end

%% add reg_pxpermm
if isfield(data,'registrationdata_scale'),
  for i = 1:numel(data),
    data(i).registrationdata_pxpermm = 1./data(i).registrationdata_scale;
  end
end

%% add "set" -- super-experiment

[line_names,~,lineidx] = unique({data.line_name});
sets = nan(1,numel(data));
seti = 0;
for linei = 1:numel(line_names),
  expidx1 = find(lineidx==linei);
  [rigs,~,rigidx] = unique([data(expidx1).rig]);
  for rigi = 1:numel(rigs),
    expidx2 = expidx1(rigidx==rigi);

    % sort by datetime
    [exp_datenum,order] = sort(datenum({data(expidx2).exp_datetime},'yyyymmddTHHMMSS'));
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
  %fprintf('%s\n',data(expidx).experiment_name);
  [~,order] = sort({data(expidx).exp_datetime});
  min_datetime = data(expidx(order(1))).exp_datetime;
  set_name = sprintf('%s__Rig%d__%s',data(expidx(1)).line_name,...
    data(expidx(1)).rig,...
    min_datetime);
  for i = expidx(:)',
    data(i).set = set_name;
  end
end

%% add in normalized bkgd stats

fns_bkgd = {'bkgd_diagnostics_mean_bkgdcenter','bkgd_diagnostics_mean_bkgdcenter_llr',...
  'registrationdata_circleCenterX','registrationdata_circleCenterY','registrationdata_bowlMarkerTheta',...
  'temperature_diagnostics_mean','temperature_diagnostics_max'};

if any(isfield(data,fns_bkgd)) && isfield(data,'plate_bowl'),
  
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
  [ism,platebowlidx] = ismember({data.plate_bowl},platebowlstats.platebowls);
  if any(~ism),
    % compute stats from this data if no stat in file
    new_plate_bowls = setdiff({data.plate_bowl},platebowlstats.platebowls);
    n0 = numel(platebowlstats.platebowls);
    for i = 1:numel(fns_bkgd),
      fn = fns_bkgd{i};
      if ~isfield(data,fn);
        continue;
      end
      y_bkgd = nan(1,numel(data));
      goodidx = ~cellfun(@isempty,{data.(fn)});
      y_bkgd(goodidx) = [data(goodidx).(fn)];
      for j = 1:numel(new_plate_bowls),
        idx = strcmp({data.plate_bowl},new_plate_bowls{j}) & goodidx;
        platebowlstats.platebowls{n0+j} = new_plate_bowls{j};
        platebowlidx(idx) = n0+j;
        platebowlstats.(fn)(n0+j) = median(y_bkgd(idx));
      end
    end
  end
  for i = 1:numel(fns_bkgd),
    fn = fns_bkgd{i};
    if isfield(data,fn) && isfield(platebowlstats,fn);
      fnz = [fn,'_norm_platebowl'];
      for j = 1:numel(data),
        data(j).(fnz) = data(j).(fn) - platebowlstats.(fn)(platebowlidx(j));
      end
    end
  end
end

if all(isfield(data,{'registrationdata_circleCenterX','registrationdata_circleCenterY'})),
  for i = 1:numel(data),
    data(i).registrationdata_circleCenterMaxXY = ...   
      max(abs(data(i).registrationdata_circleCenterX - CIRCLECENTERX),...
      abs(data(i).registrationdata_circleCenterY - CIRCLECENTERY));
  end
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
    if isfilesystempath && ~strcmpi(datamerge(i).file_system_path,'NULL'),
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
if removemissingdata,
  ismissingdata = false(size(datamerge));
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
else
  % replace missing fields
  fns = setdiff(fieldnames(datamerge),ignore_missingdata_fns);
  for i = 1:numel(fns),
    fn = fns{i};
    ismissingdatacurr = cellfun(@isempty,{datamerge.(fn)});
    if ~any(ismissingdatacurr), continue; end
    isnumericdatacurr = cellfun(@isnumeric,{datamerge.(fn)});
    isscalardatacurr = cellfun(@numel,{datamerge.(fn)}) <= 1;
    if all( (isnumericdatacurr&isscalardatacurr) | ismissingdatacurr),
      missingidx = find(ismissingdatacurr);
      for j = missingidx(:)',
        datamerge(j).(fn) = nan;
      end
      warning(['The following experiments are missing data for %s:',sprintf('\n%s',datamerge(ismissingdatacurr).experiment_name)],fn);
    end
  end
  
end


end

function in = convert2numeric(in)

numeric_fields = {...
  'experiment_id'...
  'flag_aborted'...
  'flag_review'...
  'flag_redo'...
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
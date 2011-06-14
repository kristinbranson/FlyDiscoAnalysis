function FlyBowlMetadataCheck(varargin)

[outfilename,...
  data_type,max_maxdiff_sorting_time_minutes,max_maxdiff_starvation_time_minutes,...
  max_maxdiff_exp_datetime_minutes,...
  min_mindt_exp_datetime_diff_sets,...
  leftovers] = ...
  myparse_nocheck(...
  'outfilename','',...
  'rootdatadir','/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data',...
  'data_type','QuickStats_BackSubStats_meanNConnComps',...
  'max_maxdiff_sorting_time_minutes',5,...
  'max_maxdiff_exp_datetime_minutes',10,...
  'max_maxdiff_starvation_time_minutes',5,...
  'min_mindt_exp_datetime_diff_sets',10);

%% constants

allowed_bowls = {'A','B','C','D'};
nallowed_bowls = numel(allowed_bowls);
allowed_rigs = [1,2];
allowed_plates = [10,14];
allowed_top_plates = [1,2];

%% get metadata
data = SAGEGetBowlData('data_type',data_type,leftovers{:});

%% within-set checks

[sets,~,setidx] = unique({data.set});
nsets = numel(sets);
msgs = cell(1,nsets);
iserror = false(1,nsets);

for seti = 1:nsets,

  msgs{seti} = {};
  
  idxcurr = setidx==seti;
  
  % check if bowls are unique
  bowls = {data(idxcurr).bowl};
  [isallowedbowl,bowlidx] = ismember(bowls,allowed_bowls);
  bowl_counts = hist(bowlidx,1:nallowed_bowls);
  if any(bowl_counts > 1),
    msgs{seti}{end+1} = ['The following bowls are repeated:',sprintf(' %s',allowed_bowls{bowl_counts>1})];
    iserror(seti) = true;
  end
  % check if illegal bowl values
  if ~all(isallowedbowl),
    msgs{seti}{end+1} = ['The following bowl values are not allowed:',sprintf(' %s',bowls{~isallowedbowl})];
    iserror(seti) = true;
  end
  % warning if bowls are missing
  if any(bowl_counts < 1),
    msgs{seti}{end+1} = ['The following bowls are not present:',sprintf(' %s',allowed_bowls{bowl_counts==0})];
    iserror(seti) = true;
  end

  % check if rigs are all the same
  rigs = [data(idxcurr).rig];
  if numel(unique(rigs)) ~= 1,
    msgs{seti}{end+1} = ['Multiple rigs recorded:',sprintf(' %d',unique(rigs))];
    iserror(seti) = true;
  end
  % check if illegal rig values
  [isallowedrig,idx_rig] = ismember(rigs,allowed_rigs);
  if ~all(isallowedrig),
    msgs{seti}{end+1} = ['The following rig values are not allowed:',sprintf(' %s',rigs{~isallowedrig})];
    iserror(seti) = true;
  end
  
  % check if plates are all the same
  plates = [data(idxcurr).plate];
  if numel(unique(plates)) ~= 1,
    msgs{seti}{end+1} = ['Multiple plates recorded:',sprintf(' %d',unique(plates))];
    iserror(seti) = true;
  end
  % check if illegal plate values
  [isallowedplate,idx_plate] = ismember(plates,allowed_plates);
  if ~all(isallowedplate),
    msgs{seti}{end+1} = ['The following plate values are not allowed:',sprintf(' %s',plates{~isallowedplate})];
    iserror(seti) = true;
  end
  
  % check if top_plates are all the same
  top_plates = [data(idxcurr).top_plate];
  if numel(unique(top_plates)) ~= 1,
    msgs{seti}{end+1} = ['Multiple top_plates recorded:',sprintf(' %d',unique(top_plates))];
    iserror(seti) = true;
  end
  % check if illegal top_plate values
  [isallowedtop_plate,idx_top_plate] = ismember(top_plates,allowed_top_plates);
  if ~all(isallowedtop_plate),
    msgs{seti}{end+1} = ['The following top_plate values are not allowed:',sprintf(' %s',top_plates{~isallowedtop_plate})];
    iserror(seti) = true;
  end
  
  % check if rig, plate, top_plate match
  if any(idx_rig(isallowedrig&isallowedplate) ~= idx_plate(isallowedrig&isallowedplate)),
    msgs{seti}{end+1} = 'Rigs and plates don''t match.';
    iserror(seti) = true;
  end
  if any(idx_rig(isallowedrig&isallowedtop_plate) ~= idx_top_plate(isallowedrig&isallowedtop_plate)),
    msgs{seti}{end+1} = 'Rigs and top_plates don''t match.';
    iserror(seti) = true;
  end
  if any(idx_plate(isallowedplate&isallowedtop_plate) ~= idx_top_plate(isallowedplate&isallowedtop_plate)),
    msgs{seti}{end+1} = 'Plates and top_plates don''t match.';
    iserror(seti) = true;
  end

  % check if experimenter is the same
  experimenters = unique({data(idxcurr).experimenter});
  if numel(experimenters) ~= 1,
    msgs{seti}{end+1} = ['Multiple experimenters recorded:',sprintf(' %s',experimenters{:})];
    iserror(seti) = true;
  end
  
  % check if handler_sorting is the same
  handler_sortings = unique({data(idxcurr).handler_sorting});
  if numel(handler_sortings) ~= 1,
    msgs{seti}{end+1} = ['Multiple handler_sortings recorded:',sprintf(' %s',handler_sortings{:})];
    iserror(seti) = true;
  end
  
  % check if handler_cross is the same
  handler_crosses = unique({data(idxcurr).handler_cross});
  if numel(handler_crosses) ~= 1,
    msgs{seti}{end+1} = ['Multiple handler_crosses recorded:',sprintf(' %s',handler_crosses{:})];
    iserror(seti) = true;
  end
  
  % check if handler_starvation is the same
  handler_starvations = unique({data(idxcurr).handler_starvation});
  if numel(handler_starvations) ~= 1,
    msgs{seti}{end+1} = ['Multiple handler_starvations recorded:',sprintf(' %s',handler_starvations{:})];
    iserror(seti) = true;
  end
  
  % check if cross_date is the same
  cross_dates = unique({data(idxcurr).cross_date});
  if numel(cross_dates) ~= 1,
    msgs{seti}{end+1} = ['Multiple cross_dates recorded:',sprintf(' %s',cross_dates{:})];
    iserror(seti) = true;
  end
  
  % check if flip_date is the same
  flip_dates = unique({data(idxcurr).flip_date});
  if numel(flip_dates) ~= 1,
    msgs{seti}{end+1} = ['Multiple flip_dates recorded:',sprintf(' %s',flip_dates{:})];
    iserror(seti) = true;
  end
  
  % sorting_time should be about the same
  hours_sorted = [data(idxcurr).hours_sorted];
  sorting_datenum = datenum({data(idxcurr).exp_datetime},'yyyymmddTHHMMSS')' - hours_sorted/24;
  maxdiff_sorting_time_minutes = (max(sorting_datenum) - min(sorting_datenum))*24*60;
  if maxdiff_sorting_time_minutes > maxdiff_sorting_time_minutes,
    msgs{seti}{end+1} = sprintf('Difference between first and last sorting time = %f > %f minutes',...
      maxdiff_sorting_time_minutes,max_maxdiff_sorting_time_minutes);
    iserror(seti) = true;
  end
  
  % starvation_time should be about the same
  hours_starved = [data(idxcurr).hours_starved];
  starvation_datenum = datenum({data(idxcurr).exp_datetime},'yyyymmddTHHMMSS')' - hours_starved/24;
  maxdiff_starvation_time_minutes = (max(starvation_datenum) - min(starvation_datenum))*24*60;
  if maxdiff_starvation_time_minutes > maxdiff_starvation_time_minutes,
    msgs{seti}{end+1} = sprintf('Difference between first and last starvation time = %f > %f minutes',...
      maxdiff_starvation_time_minutes,max_maxdiff_starvation_time_minutes);
    iserror(seti) = true;
  end
  
  % exp_datetime should be about the same
  exp_datenum = datenum({data(idxcurr).exp_datetime},'yyyymmddTHHMMSS')';
  maxdiff_exp_datetime_minutes = (max(exp_datenum) - min(exp_datenum))*24*60;
  if maxdiff_exp_datetime_minutes > maxdiff_exp_datetime_minutes,
    msgs{seti}{end+1} = sprintf('Difference between first and last experiment starts = %f > %f minutes',...
      maxdiff_exp_datetime_minutes,max_maxdiff_exp_datetime_minutes);
    iserror(seti) = true;
  end
  
  % wish list should be the same
  wish_lists = str2double({data(idxcurr).wish_list});
  if numel(unique(wish_lists)) ~= 1,
    msgs{seti}{end+1} = ['Multiple wish_lists recorded:',sprintf(' %d',unique(wish_lists))];
    iserror(seti) = true;
  end

  % barcode not currently in flattened view
  % barcode == -1 and not control
  barcodes = str2double({data(idxcurr).cross_barcode});
  isolympiad = ~strcmp({data(idxcurr).screen_type},'non-olympiad');
  iscontrol = strcmp({data(idxcurr).screen_reason},'control');
  bad_barcode = barcodes(isolympiad & ~iscontrol) <= 0;
  if any(bad_barcode),
    msgs{seti}{end+1} = sprintf('Barcode is not set for %d experiments.',nnz(bad_barcode));
    iserror(seti) = true;
  end
  
end

%% across-set checks

% look for different sets on the same rig at about the same time
for seti = 1:nsets-1,
  idxi = setidx==seti;
  same_rig_time = [];
  rigsi = unique([data(idxi).rig]);
  exp_datenumi = datenum({data(idxi).exp_datetime},'yyyymmddTHHMMSS')';
  for setj = seti+1:nsets,
    idxj = setidx==setj;
    rigsj = unique([data(idxj).rig]);
    exp_datenumj = datenum({data(idxj).exp_datetime},'yyyymmddTHHMMSS')';
    mindt = min(pdist2(exp_datenumi',exp_datenumj','euclidean','smallest',1));
    mindt_minutes = mindt*24*60;
    if mindt_minutes < min_mindt_exp_datetime_diff_sets && ~isempty(intersect(rigsi,rigsj)),
      same_rig_time(end+1) = setj; %#ok<AGROW>
    end
  end
  if ~isempty(same_rig_time),
    msgs{seti}{end+1} = ['The following sets occur at approximately the same time on the same rig:',...
      sprintf(' %s',sets{same_rig_time})];
    iserror(seti) = true;
  end
end

%% print results

if isempty(outfilename),
  fid = 1;
else
  fid = fopen(outfilename,'w');
end

sorted_exp_datetimes = sort({data.exp_datetime});
fprintf(fid,'Metadata check for the following %d experiments collected between %s and %s:\n',...
  numel(data),sorted_exp_datetimes{1},sorted_exp_datetimes{end});
fprintf(fid,'%s\n',data.experiment_name);
fprintf(fid,'\n');

fprintf(fid,'%d within-set errors:\n',nnz(iserror));
if ~any(iserror),
  fprintf(fid,'**NONE**\n');
else
  for seti = find(iserror),
    fprintf(fid,'Set %s:\n',sets{seti});
    fprintf(fid,'  %s\n',msgs{seti}{:});
  end
end

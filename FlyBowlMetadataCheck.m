function [data,iserror,msgs,iserror_exp,msgs_exp,iserror_date,msgs_date] = ...
  FlyBowlMetadataCheck(varargin)

[outfilename,...
  data_type,max_maxdiff_sorting_time_minutes,max_maxdiff_starvation_time_minutes,...
  max_maxdiff_exp_datetime_minutes,...
  min_mindt_exp_datetime_diff_sets,...
  first_barcode_datetime,...
  max_maxdiff_sorting_time_perday_days,...
  dosetchecks,doexpchecks,dodatechecks,...
  leftovers] = ...
  myparse_nocheck(varargin,...
  'outfilename','',...
  'data_type','QuickStats_BackSubStats_meanNConnComps',...
  'max_maxdiff_sorting_time_minutes',5,...
  'max_maxdiff_exp_datetime_minutes',10,...
  'max_maxdiff_starvation_time_minutes',5,...
  'min_mindt_exp_datetime_diff_sets',10,...
  'first_barcode_datetime','20110413T000000',...
  'max_maxdiff_sorting_time_perday_days',4.5,...
  'dosetchecks',true,...
  'doexpchecks',true,...
  'dodatechecks',true);

%% constants

datetime_format = 'yyyymmddTHHMMSS';
allowed_bowls = {'A','B','C','D'};
nallowed_bowls = numel(allowed_bowls);
allowed_rigs = [1,2];
allowed_plates_per_rig = {[10,11,15],[14,17]};
allowed_top_plates_per_rig = {[1,5],[2,3,4]};
allowed_plates = unique([allowed_plates_per_rig{:}]);
allowed_top_plates = unique([allowed_top_plates_per_rig{:}]);
idx_plate_to_rig = nan(size(allowed_plates));
for i = 1:numel(allowed_rigs),
  for j = 1:numel(allowed_plates_per_rig{i}),
    k = allowed_plates == allowed_plates_per_rig{i}(j);
    idx_plate_to_rig(k) = allowed_rigs(i);
  end
end
idx_top_plate_to_rig = nan(size(allowed_top_plates));
for i = 1:numel(allowed_rigs),
  for j = 1:numel(allowed_top_plates_per_rig{i}),
    k = allowed_top_plates == allowed_top_plates_per_rig{i}(j);
    idx_top_plate_to_rig(k) = allowed_rigs(i);
  end
end
%allowed_plates = [10,14,15,16,17,18];
%allowed_top_plates = [1,2,3,4,5,6];
first_barcode_datenum = datenum(first_barcode_datetime,datetime_format);
start_of_day = 7/24; % 7am
end_of_day = 20/24; % 8pm


%% get metadata
data = SAGEGetBowlData('data_type',data_type,'dataset','score',leftovers{:});

% sort by date
[~,order] = sort({data.exp_datetime});
data = data(order);

fprintf('Data pulled ... checking ...\n');

exp_datenums = datenum({data.exp_datetime},datetime_format)';
cross_datenums = datenum({data.cross_date},datetime_format)';
% flip date seems to not be set all the time -- set to nan here
flip_datenums = nan(1,numel(data));
for i = 1:numel(data),
  try
    flip_datenums(i) = datenum(data(i).flip_date,datetime_format);
  catch %#ok<CTCH>
  end
end

hours_sorted = [data.hours_sorted];
sorting_datenums = exp_datenums - hours_sorted/24;

hours_starved = [data.hours_starved];
starvation_datenums = exp_datenums - hours_starved/24;

%% within-experiment checks

msgs_exp = cell(1,numel(data));
iserror_exp = false(1,numel(data));

if doexpchecks,

for i = 1:numel(data),

  msgs_exp{i} = {};
  % apparatus name matches individual fields
  
  % parse apparatus name
  % Rig2__Plate14__Lid02__BowlD__Camera0053300063A94001__Computerbransonlab-ww4__HardDriveInternal_D
  m = regexp(data(i).apparatus_id,...
    '^Rig(?<rig>\s*\d+\s*)__Plate(?<plate>\s*\d+\s*)__Lid(?<lid>\s*\d+\s*)__Bowl(?<bowl>.+)__Camera(?<camera>.+)__Computer(?<computer>.+)__HardDrive(?<harddrive>.+)$','names');
  if isempty(m),
    msgs_exp{i}{end+1} = sprintf('Could not parse apparatus %s',data(i).apparatus_id);
    iserror_exp(i) = true;
  else
    if ~strcmp(sprintf('%d',data(i).rig),m.rig),
      msgs_exp{i}{end+1} = sprintf('Rig %d does not match apparatus %s',data(i).rig,data(i).apparatus_id);
      iserror_exp(i) = true;
    end
    if ~strcmp(sprintf('%02d',data(i).plate),m.plate),
      msgs_exp{i}{end+1} = sprintf('plate %d does not match apparatus %s',data(i).plate,data(i).apparatus_id);
      iserror_exp(i) = true;
    end
    if ~strcmp(sprintf('%02d',data(i).top_plate),m.lid),
      msgs_exp{i}{end+1} = sprintf('lid %d does not match apparatus %s',data(i).top_plate,data(i).apparatus_id);
      iserror_exp(i) = true;
    end
    if ~strcmp(data(i).bowl,m.bowl),
      msgs_exp{i}{end+1} = sprintf('bowl %s does not match apparatus %s',data(i).bowl,data(i).apparatus_id);
      iserror_exp(i) = true;
    end
    if ~strcmp(data(i).camera,m.camera),
      msgs_exp{i}{end+1} = sprintf('camera %s does not match apparatus %s',data(i).camera,data(i).apparatus_id);
      iserror_exp(i) = true;
    end
    if ~strcmp(data(i).computer,m.computer),
      msgs_exp{i}{end+1} = sprintf('computer %s does not match apparatus %s',data(i).computer,data(i).apparatus_id);
      iserror_exp(i) = true;
    end
    if ~strcmp(data(i).harddrive,m.harddrive),
      msgs_exp{i}{end+1} = sprintf('harddrive %s does not match apparatus %s',data(i).harddrive,data(i).apparatus_id);
      iserror_exp(i) = true;
    end
  end

  % file system path matches individual fields

  % parse file system path
  [m,success1] = parseExpDir(data(i).file_system_path);
  if ~success1,
    msgs_exp{i}{end+1} = sprintf('Error parsing file system path %s',data(i).file_system_path);
    iserror_exp(i) = true;
  else
    if ~strcmp(sprintf('%d',data(i).rig),m.rig),
      msgs_exp{i}{end+1} = sprintf('Rig %d does not match file_system_path %s',data(i).rig,data(i).file_system_path);
      iserror_exp(i) = true;
    end
    if ~strcmp(sprintf('%02d',data(i).plate),m.plate),
      msgs_exp{i}{end+1} = sprintf('plate %d does not match file_system_path %s',data(i).plate,data(i).file_system_path);
      iserror_exp(i) = true;
    end
    if ~strcmp(data(i).bowl,m.bowl),
      msgs_exp{i}{end+1} = sprintf('bowl %s does not match file_system_path %s',data(i).bowl,data(i).file_system_path);
      iserror_exp(i) = true;
    end
    if ~strcmp(data(i).line_name,m.line),
      msgs_exp{i}{end+1} = sprintf('line %s does not match file_system_path %s',data(i).line_name,data(i).file_system_path);
      iserror_exp(i) = true;
    end
  end

  % experiment name matches file_system_path
  if ~strcmp(data(i).experiment_name(1:numel('FlyBowl_')),'FlyBowl_'),
    msgs_exp{i}{end+1} = sprintf('experiment name %s does not begin with FlyBowl_',data(i).experiment_name);
    iserror_exp(i) = true;
  else
    [m,success1] = parseExpDir(data(i).experiment_name(numel('FlyBowl_')+1:end));
    
    if ~success1,
      msgs_exp{i}{end+1} = sprintf('Error parsing experiment_name %s',data(i).experiment_name);
      iserror_exp(i) = true;
    else
      if ~strcmp(sprintf('%d',data(i).rig),m.rig),
        msgs_exp{i}{end+1} = sprintf('Rig %d does not match experiment_name %s',data(i).rig,data(i).experiment_name);
        iserror_exp(i) = true;
      end
      if ~strcmp(sprintf('%02d',data(i).plate),m.plate),
        msgs_exp{i}{end+1} = sprintf('plate %d does not match experiment_name %s',data(i).plate,data(i).experiment_name);
        iserror_exp(i) = true;
      end
      if ~strcmp(data(i).bowl,m.bowl),
        msgs_exp{i}{end+1} = sprintf('bowl %s does not match experiment_name %s',data(i).bowl,data(i).experiment_name);
        iserror_exp(i) = true;
      end
      if ~strcmp(data(i).line_name,m.line),
        msgs_exp{i}{end+1} = sprintf('line %s does not match experiment_name %s',data(i).line_name,data(i).experiment_name);
        iserror_exp(i) = true;
      end
    end
  end
  
end


end


%% within-set checks

[sets,~,setidx] = unique({data.set});
nsets = numel(sets);
msgs = cell(1,nsets);
iserror = false(1,nsets);

if dosetchecks,
  
for seti = 1:nsets,

  msgs{seti} = {};
  
  idxcurr = setidx==seti;
  idxcurr1 = find(idxcurr);  
  
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
    msgs{seti}{end+1} = ['The following plate values are not allowed:',sprintf(' %d',plates(~isallowedplate))];
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
    msgs{seti}{end+1} = ['The following top_plate values are not allowed:',sprintf(' %d',top_plates(~isallowedtop_plate))];
    iserror(seti) = true;
  end
  
  
  % check if rig, plate, top_plate match
  if any(rigs(isallowedrig&isallowedplate) ~= idx_plate_to_rig(idx_plate(isallowedrig&isallowedplate))),
    msgs{seti}{end+1} = sprintf('Rigs and plates don''t match. Plates: %s, Rigs: %s',...
      sprintf('%d ',plates(isallowedrig&isallowedplate)),...
      sprintf('%d ',rigs(isallowedrig&isallowedplate)));
    iserror(seti) = true;
  end
  if any(rigs(isallowedrig&isallowedtop_plate) ~= idx_top_plate_to_rig(idx_top_plate(isallowedrig&isallowedtop_plate))),
    msgs{seti}{end+1} = sprintf('Rigs and top_plates don''t match. Top plates: %s, Rigs: %s',...
      sprintf('%d ',top_plates(isallowedrig&isallowedplate)),...
      sprintf('%d ',rigs(isallowedrig&isallowedplate)));
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
  % make sure flip dates are valid
  if any(isnan(flip_datenums(idxcurr))),
    msgs{seti}{end+1} = 'Flip dates could not be parsed.';
    iserror(seti) = true;
  end

  % sorting_time should be about the same
  maxdiff_sorting_time_minutes = (max(sorting_datenums(idxcurr)) - min(sorting_datenums(idxcurr)))*24*60;
  if maxdiff_sorting_time_minutes > maxdiff_sorting_time_minutes,
    msgs{seti}{end+1} = sprintf('Difference between first and last sorting time = %f > %f minutes',...
      maxdiff_sorting_time_minutes,max_maxdiff_sorting_time_minutes);
    iserror(seti) = true;
  end
  
  % starvation_time should be about the same
  maxdiff_starvation_time_minutes = (max(starvation_datenums(idxcurr)) - min(starvation_datenums(idxcurr)))*24*60;
  if maxdiff_starvation_time_minutes > max_maxdiff_starvation_time_minutes,
    msgs{seti}{end+1} = sprintf('Difference between first and last starvation time = %f > %f minutes',...
      maxdiff_starvation_time_minutes,max_maxdiff_starvation_time_minutes);
    iserror(seti) = true;
  end
  
  % time of day sorted should be reasonable
  sorting_tods = mod(sorting_datenums(idxcurr),1);
  badidx = sorting_tods < start_of_day | sorting_tods > end_of_day;
  if any(badidx),
    msgs{seti}{end+1} = 'The following experiments have sorting times not during the work day:';
    for j = find(badidx),
      msgs{seti}{end} = [msgs{seti}{end},sprintf('\n%s: %s',data(idxcurr1(j)).experiment_name,datestr(sorting_tods(j),'HH:MM:SS AM'))];
    end
    iserror(seti) = true;
  end
  
  % time of day starved should be reasonable
  starvation_tods = mod(starvation_datenums(idxcurr),1);
  badidx = starvation_tods < start_of_day | starvation_tods > end_of_day;
  if any(badidx),
    msgs{seti}{end+1} = 'The following experiments have starvation times not during the work day:';
    for j = find(badidx),
      msgs{seti}{end} = [msgs{seti}{end},sprintf('\n%s: %s',data(idxcurr1(j)).experiment_name,datestr(starvation_tods(j),'HH:MM:SS AM'))];
    end
    iserror(seti) = true;
  end
  
  % exp_datetime should be about the same
  exp_datenum = exp_datenums(idxcurr);
  maxdiff_exp_datetime_minutes = (max(exp_datenum) - min(exp_datenum))*24*60;
  if maxdiff_exp_datetime_minutes > maxdiff_exp_datetime_minutes,
    msgs{seti}{end+1} = sprintf('Difference between first and last experiment starts = %f > %f minutes',...
      maxdiff_exp_datetime_minutes,max_maxdiff_exp_datetime_minutes);
    iserror(seti) = true;
  end
  
  % wish list should be the same
  wish_lists = str2double({data(idxcurr).wish_list});
  if any(isnan(wish_lists)),
    msgs{seti}{end+1} = 'Wish list is nan';
    iserror(seti) = true;
  elseif numel(unique(wish_lists)) ~= 1,
    msgs{seti}{end+1} = ['Multiple wish_lists recorded:',sprintf(' %d',unique(wish_lists))];
    iserror(seti) = true;
  end

  % barcode not currently in flattened view
  % barcode == -1 and not control
  barcodes = str2double({data(idxcurr).cross_barcode});
  isolympiad = ~strcmp({data(idxcurr).screen_type},'non-olympiad');
  iscontrol = strcmp({data(idxcurr).screen_reason},'control');
  islateenough = isempty(first_barcode_datetime) | ...
    exp_datenums(idxcurr) >= first_barcode_datenum;
  bad_barcode = barcodes(isolympiad & ~iscontrol & islateenough) <= 0;
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
  exp_datenumi = exp_datenums(idxi);
  for setj = seti+1:nsets,
    idxj = setidx==setj;
    rigsj = unique([data(idxj).rig]);
    exp_datenumj = exp_datenums(idxj);
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

end

%% experiments on the same day should have some similar metadata

exp_dates = floor(exp_datenums);
[unique_exp_dates,~,exp_date_idx] = unique(exp_dates);
ndates = numel(unique_exp_dates);
iserror_date = false(1,ndates);
msgs_date = cell(1,ndates);

if dodatechecks,

for datei = 1:numel(unique_exp_dates),
  msgs_date{datei} = {};
  idxcurr = find(exp_date_idx == datei);
  
  % cross dates should be the same on the same experiment day
  cross_dates_curr = floor(cross_datenums(idxcurr));
  [unique_cross_dates,~,cross_date_idx] = unique(cross_dates_curr);
  if numel(unique_cross_dates) > 1,
    iserror_date(datei) = true;
    mode_cross_date = unique_cross_dates(mode(cross_date_idx));
    msgs_date{datei}{end+1} = sprintf('%d / %d experiments on %s have cross date %s. The following experiments have different cross dates:',...
      nnz(cross_dates_curr == mode_cross_date),...
      numel(idxcurr),...
      datestr(exp_dates(datei),'yyyy-mm-dd'),...
      datestr(mode_cross_date,'yyyy-mm-dd'));
    for j = find(cross_dates_curr ~= mode_cross_date),
      msgs_date{datei}{end} = [msgs_date{datei}{end},sprintf('\n%s: %s',data(idxcurr(j)).experiment_name,data(idxcurr(j)).cross_date)];
    end
  end
  
  % flip dates should be the same on the same experiment day
  flip_dates_curr = floor(flip_datenums(idxcurr));
  goodidx = find(~isnan(flip_dates_curr));
  flip_dates_curr = flip_dates_curr(goodidx);
  [unique_flip_dates,~,flip_date_idx] = unique(flip_dates_curr);
  if numel(unique_flip_dates) > 1,
    iserror_date(datei) = true;
    mode_flip_date = unique_flip_dates(mode(flip_date_idx));
    msgs_date{datei}{end+1} = sprintf('%d / %d experiments on %s have flip date %s. The following experiments have different flip dates:',...
      nnz(flip_dates_curr == mode_flip_date),...
      numel(idxcurr),...
      datestr(exp_dates(datei),'yyyy-mm-dd'),...
      datestr(mode_flip_date,'yyyy-mm-dd'));
    for j = goodidx(flip_dates_curr ~= mode_flip_date),
      msgs_date{datei}{end} = [msgs_date{datei}{end},sprintf('\n%s: %s',data(idxcurr(j)).experiment_name,data(idxcurr(j)).flip_date)];
    end
  end
  
  % sorting should take place within a span of 2 days
  maxdiff_sorting = max(sorting_datenums(idxcurr)) - min(sorting_datenums(idxcurr));
  if maxdiff_sorting > max_maxdiff_sorting_time_perday_days,
    msgs_date{datei}{end+1} = sprintf('Sorting for experiments on this day took place over a span of %f days.',maxdiff_sorting);
    iserror_date(datei) = true;
  end
  
  % starvation should take place on the same day
  maxdiff_starvation = max(starvation_datenums(idxcurr)) - min(starvation_datenums(idxcurr));
  if maxdiff_starvation > 1,
    msgs_date{datei}{end+1} = sprintf('Starvation for experiments on this day took place over a span of %f days.',maxdiff_starvation);
    iserror_date(datei) = true;
  end
  
end

end

%% print results

if isempty(outfilename),
  fid = 1;
else
  fid = fopen(outfilename,'w');
end

sorted_exp_datetimes = sort({data.exp_datetime});
fprintf(fid,'Metadata check for %d experiments collected between %s and %s:\n\n',...
  numel(data),sorted_exp_datetimes{1},sorted_exp_datetimes{end});
%fprintf(fid,'%s\n',data.experiment_name);
%fprintf(fid,'\n');

if doexpchecks,

fprintf(fid,'\n%d/%d experiment consistency errors:\n',nnz(iserror_exp),numel(iserror_exp));
if ~any(iserror_exp),
  fprintf(fid,'**NONE**\n');
else
  for expi = find(iserror_exp),
    fprintf(fid,'Experiment %s:\n',data(expi).experiment_name);
    fprintf(fid,'  %s\n',msgs_exp{expi}{:});
  end
end

end

if dosetchecks,

fprintf(fid,'\n%d/%d set errors:\n',nnz(iserror),numel(iserror));
if ~any(iserror),
  fprintf(fid,'**NONE**\n');
else
  for seti = find(iserror),
    fprintf(fid,'Set %s:\n',sets{seti});
    fprintf(fid,'  %s\n',msgs{seti}{:});
  end
end

end

if dodatechecks,
  
fprintf(fid,'\n%d/%d day errors:\n',nnz(iserror_date),numel(iserror_date));
for datei = find(iserror_date),
  fprintf(fid,'Experiment date %s:\n',datestr(unique_exp_dates(datei),'yyyy-mm-dd'));
  fprintf(fid,'  %s\n',msgs_date{datei}{:});
end
  
end

if fid > 1,
  fclose(fid);
end
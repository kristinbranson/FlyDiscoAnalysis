function [success, iserror, error_or_warning_messages] = check_for_nans_in_tracking(trx_file_name, has_wing_info, success, iserror, error_or_warning_messages)  

if ~exist(trx_file_name,'file'),
  error_or_warning_messages{end+1} = sprintf('File %s does not exist', trx_file_name);
  success = false;
  iserror(category2idx.missing_registration_files) = true;
  return
end

% load trajectories
[trx,~,succeeded] = load_tracks(trx_file_name) ;
if ~succeeded,
  error_or_warning_messages{end+1} = sprintf('Could not load trajectories from file %s',trx_file_name);
  success = false;
  iserror(category2idx.completed_checks_other) = true;
  return
end

% Check loaded trajectories for nan's
for fly = 1:numel(trx),
  if has_wing_info ,
    badidx = ...
      isnan(trx(fly).x) | ...
      isnan(trx(fly).y) | ...
      isnan(trx(fly).a) | ...
      isnan(trx(fly).b) | ...
      isnan(trx(fly).theta) | ...
      isnan(trx(fly).wing_anglel) | ...
      isnan(trx(fly).wing_angler) ;
  else
    badidx = ...
      isnan(trx(fly).x) | ...
      isnan(trx(fly).y) | ...
      isnan(trx(fly).a) | ...
      isnan(trx(fly).b) | ...
      isnan(trx(fly).theta) ;
  end
  if any(badidx),
    [starts,ends] = get_interval_ends(badidx);
    starts = starts - trx(fly).off;
    ends = ends - trx(fly).off;
    for se = 1:numel(starts)
      error_or_warning_messages{end+1} = [sprintf('Trajectory %d has NaNs in frames',fly),sprintf(' %d-%d',[starts(se),ends(se)]')]; %#ok<AGROW>
    end
    success = false;
    iserror(category2idx.flytracker_nans) = true;
  end
end

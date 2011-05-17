function data = AddExtraData(data,statnames,registration_params)

nexpdirs = numel(data);

% add in mean_nsplit
if ismember('ctrax_diagnostics_mean_nsplit',statnames),  
  for i = 1:nexpdirs,
    if isfield(data,'ctrax_diagnostics_sum_nsplit') && isfield(data,'ctrax_diagnostics_nlarge_split'),
      data(i).ctrax_diagnostics_mean_nsplit = ...
        data(i).ctrax_diagnostics_sum_nsplit / data(i).ctrax_diagnostics_nlarge_split;
    end
  end
end

% add in nframes_not_tracked
if ismember('ctrax_diagnostics_nframes_not_tracked',statnames),
  for i = 1:nexpdirs,
    if isfield(data,'ufmf_diagnostics_summary_nFrames') && isfield(data,'ctrax_diagnostics_nframes_analyzed'),
      data(i).ctrax_diagnostics_nframes_not_tracked = ...
        data(i).ufmf_diagnostics_summary_nFrames - data(i).ctrax_diagnostics_nframes_analyzed;
    end
  end
end

% mean x-position
if ismember('mean_x_mm',statnames),
  if isfield(data,'stats_perframe_x_mm'),
    for i = 1:nexpdirs,
      data(i).mean_x_mm = data(i).stats_perframe_x_mm.meanmean_perexp.flyany_frameany;
    end
  end
end

% mean y-position
if ismember('mean_x_mm',statnames),
  if isfield(data,'stats_perframe_y_mm'),
    for i = 1:nexpdirs,
      data(i).mean_y_mm = data(i).stats_perframe_y_mm.meanmean_perexp.flyany_frameany;
    end
  end
end

% pxpermm
if ismember('registrationdata_pxpermm',statnames),
  if isfield(data,'registrationdata_circleRadius') && ...
      exist('registration_params','var'),
    for i = 1:nexpdirs,
      data(i).registrationdata_pxpermm = data(i).registrationdata_circleRadius / registration_params.circleRadius_mm;
    end
  end
end

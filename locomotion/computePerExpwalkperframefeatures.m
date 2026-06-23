function [perexpperframefeatures] = computePerExpwalkperframefeatures(perwalk_metrics,perfly_metrics,fieldname)
% frm_mean_exp = mean of all valid frms across the movie - concatenate from
% perwalk_metrics - compute mean, std, n from frame data

% walk_mean_exp = mean of valid walks from across frames
% from perwalk_metrics - compute mean, std, n across walk means

% mean_frm_mean_fly = mean of mean_frm_fly data from perfly_metrics
% mean_walk_mean_fly = mean of mean_walk_fly data from perflymetrics


% frm_mean_exp (+ co-indexed smoothed speed for speed-binning)
curr_all = {};
spd_all = {};
data =[];
have_speed = isfield(perwalk_metrics, 'smoothvelmag_ctr');
for w = 1:numel(perwalk_metrics)
    curr_all{w} = [perwalk_metrics(w).(fieldname).data];
    if have_speed
        spd_all{w} = [perwalk_metrics(w).smoothvelmag_ctr.data];
    end
end
curr_all = horzcat(curr_all{:});
mask = ~isnan(curr_all);
data = curr_all(mask);
perexpperframefeatures.(fieldname).frm_data_exp = data;
perexpperframefeatures.(fieldname).frm_mean_exp = mean(data);
perexpperframefeatures.(fieldname).frm_n_exp = numel(data);
perexpperframefeatures.(fieldname).frm_std_exp = std(data);
% co-indexed smoothed speed (same NaN mask as frm_data_exp); [] if unavailable
% or per-walk length mismatch (then binning is skipped downstream).
if have_speed
    spd_all = horzcat(spd_all{:});
    if numel(spd_all) == numel(curr_all)
        perexpperframefeatures.(fieldname).frm_speed_exp = spd_all(mask);
    else
        perexpperframefeatures.(fieldname).frm_speed_exp = [];
    end
else
    perexpperframefeatures.(fieldname).frm_speed_exp = [];
end

% walk_mean_exp
curr_all = {};
data =[];
for w = 1:numel(perwalk_metrics)
    curr_all{w} = [perwalk_metrics(w).(fieldname).mean];
end
curr_all = horzcat(curr_all{:});
data = curr_all(~isnan(curr_all));
perexpperframefeatures.(fieldname).walk_data_exp = data;
perexpperframefeatures.(fieldname).walk_mean_exp = mean(data);
perexpperframefeatures.(fieldname).walk_n_exp = numel(data);
perexpperframefeatures.(fieldname).walk_std_exp = std(data);


% mean_frm_mean_fly
curr_all = {};
data =[];
for f = 1:numel(perfly_metrics)
    curr_all{f} = [perfly_metrics(f).(fieldname).frm_mean_fly];
end
curr_all = horzcat(curr_all{:});
data = curr_all(~isnan(curr_all));
perexpperframefeatures.(fieldname).data_frm_mean_fly = data;
perexpperframefeatures.(fieldname).mean_frm_mean_fly = mean(data);
perexpperframefeatures.(fieldname).n_frm_mean_fly = numel(data);
perexpperframefeatures.(fieldname).std_frm_mean_fly = std(data);

% mean_walk_mean_fly
curr_all = {};
data =[];
for f = 1:numel(perfly_metrics)
    curr_all{f} = [perfly_metrics(f).(fieldname).walk_mean_fly];   
end
curr_all = horzcat(curr_all{:});
data = curr_all(~isnan(curr_all));
perexpperframefeatures.(fieldname).data_walk_mean_fly = data;
perexpperframefeatures.(fieldname).mean_walk_mean_fly = mean(data);
perexpperframefeatures.(fieldname).n_walk_mean_fly = numel(data);
perexpperframefeatures.(fieldname).std_walk_mean_fly = std(data);





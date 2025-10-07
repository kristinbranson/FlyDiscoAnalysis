function [perexpperframefeatures] = computePerExpwalkperframefeatures(perwalk_metrics,perfly_metrics,fieldname)
% frm_mean_exp = mean of all valid frms across the movie - concatenate from
% perwalk_metrics - compute mean, std, n from frame data

% walk_mean_exp = mean of valid walks from across frames
% from perwalk_metrics - compute mean, std, n across walk means

% mean_frm_mean_fly = mean of mean_frm_fly data from perfly_metrics
% mean_walk_mean_fly = mean of mean_walk_fly data from perflymetrics


% frm_mean_exp
curr_all = {};
data =[];
for w = 1:numel(perwalk_metrics)
    curr_all{w} = [perwalk_metrics(w).(fieldname).data];
end
curr_all = horzcat(curr_all{:});
data = curr_all(~isnan(curr_all));
perexpperframefeatures.(fieldname).frm_data_exp = data;
perexpperframefeatures.(fieldname).frm_mean_exp = mean(data);
perexpperframefeatures.(fieldname).frm_n_exp = numel(data);
perexpperframefeatures.(fieldname).frm_std_exp = std(data);

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





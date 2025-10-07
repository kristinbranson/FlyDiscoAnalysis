function [perexpphasefeatures] = computePerExpphasefeatures(perwalk_metrics,perfly_metrics,fieldname)
% frm_mean_exp = mean of all valid frms across the movie - concatenate from
% perwalk_metrics - compute mean, std, n from frame data

% walk_mean_exp = mean of valid walks from across frames
% from perwalk_metrics - compute mean, std, n across walk means

% mean_frm_mean_fly = mean of mean_frm_fly data from perfly_metrics
% mean_walk_mean_fly = mean of mean_walk_fly data from perflymetrics

subfieldnames = fields(perwalk_metrics(1).(fieldname));

for sf = 1:numel(subfieldnames)

    if strcmp(subfieldnames{sf},'phasedata')
        curr_all = {};
        for w = 1:numel(perwalk_metrics(w))
            curr_all{w} = [perwalk_metrics(w).(fieldname).(subfieldnames{sf})];
        end
        perexpphasefeatures.(fieldname).(subfieldnames{sf}) = horzcat(curr_all{:});
    else
        

        % frm_mean_exp
        curr_all = {};
        data =[];
        for w = 1:numel(perwalk_metrics)
            curr_all{w} = [perwalk_metrics(w).(fieldname).(subfieldnames{sf}).data];
        end
        curr_all = horzcat(curr_all{:});
        data = curr_all(~isnan(curr_all));
        perexpphasefeatures.(fieldname).frm_data_exp = data;
        perexpphasefeatures.(fieldname).frm_mean_exp = mean(data);
        perexpphasefeatures.(fieldname).frm_n_exp = numel(data);
        perexpphasefeatures.(fieldname).frm_std_exp = std(data);

        % walk_mean_exp
        curr_all = {};
        data =[];
        for w = 1:numel(perwalk_metrics)
            curr_all{w} = [perwalk_metrics(w).(fieldname).(subfieldnames{sf}).mean];
        end
        curr_all = horzcat(curr_all{:});
        data = curr_all(~isnan(curr_all));
        perexpphasefeatures.(fieldname).walk_data_exp = data;
        perexpphasefeatures.(fieldname).walk_mean_exp = mean(data);
        perexpphasefeatures.(fieldname).walk_n_exp = numel(data);
        perexpphasefeatures.(fieldname).walk_std_exp = std(data);


        % mean_frm_mean_fly
        curr_all = {};
        data =[];
        for f = 1:numel(perfly_metrics)
            curr_all{f} = [perfly_metrics(f).(fieldname).(subfieldnames{sf}).frm_mean_fly];
        end
        curr_all = horzcat(curr_all{:});
        data = curr_all(~isnan(curr_all));
        perexpphasefeatures.(fieldname).data_frm_mean_fly = data;
        perexpphasefeatures.(fieldname).mean_frm_mean_fly = mean(data);
        perexpphasefeatures.(fieldname).n_frm_mean_fly = numel(data);
        perexpphasefeatures.(fieldname).std_frm_mean_fly = std(data);


        % mean_walk_mean_fly
        curr_all = {};
        data =[];
        for f = 1:numel(perfly_metrics)
            curr_all{f} = [perfly_metrics(f).(fieldname).(subfieldnames{sf}).walk_mean_fly];
        end
        curr_all = horzcat(curr_all{:});
        data = curr_all(~isnan(curr_all));
        perexpphasefeatures.(fieldname).data_walk_mean_fly = data;
        perexpphasefeatures.(fieldname).mean_walk_mean_fly = mean(data);
        perexpphasefeatures.(fieldname).n_walk_mean_fly = numel(data);
        perexpphasefeatures.(fieldname).std_walk_mean_fly = std(data);

    end

end




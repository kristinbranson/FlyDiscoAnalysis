function [perexpphasefeatures] = computePerExpphasefeatures(perwalk_metrics,perfly_metrics,fieldname)
% frm_mean_exp = mean of all valid frms across the movie - concatenate from
% perwalk_metrics - compute mean, std, n from frame data

% walk_mean_exp = mean of valid walks from across frames
% from perwalk_metrics - compute mean, std, n across walk means

% mean_frm_mean_fly = mean of mean_frm_fly data from perfly_metrics
% mean_walk_mean_fly = mean of mean_walk_fly data from perflymetrics

subfieldnames = fields(perwalk_metrics(1).(fieldname));
perexpphasefeatures = struct;

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
        perexpphasefeatures.(fieldname).(subfieldnames{sf}).frm_data_exp = data;
        perexpphasefeatures.(fieldname).(subfieldnames{sf}).frm_mean_exp = circ_mean(data);
        perexpphasefeatures.(fieldname).(subfieldnames{sf}).frm_n_exp = numel(data);
        perexpphasefeatures.(fieldname).(subfieldnames{sf}).frm_std_exp = circ_std(data);
       

        % walk_mean_exp
        curr_all = {};
        data =[];
        for w = 1:numel(perwalk_metrics)
            curr_all{w} = [perwalk_metrics(w).(fieldname).(subfieldnames{sf}).mean];
        end
        curr_all = horzcat(curr_all{:});
        data = curr_all(~isnan(curr_all));
        perexpphasefeatures.(fieldname).(subfieldnames{sf}).walk_data_exp = data;
        perexpphasefeatures.(fieldname).(subfieldnames{sf}).walk_mean_exp = circ_mean(data);
        perexpphasefeatures.(fieldname).(subfieldnames{sf}).walk_n_exp = numel(data);
        perexpphasefeatures.(fieldname).(subfieldnames{sf}).walk_std_exp = circ_std(data);


        % mean_frm_mean_fly
        curr_all = {};
        data =[];
        for f = 1:numel(perfly_metrics)
            curr_all{f} = [perfly_metrics(f).(fieldname).(subfieldnames{sf}).frm_mean_fly];
        end
        curr_all = horzcat(curr_all{:});
        data = curr_all(~isnan(curr_all));
        perexpphasefeatures.(fieldname).(subfieldnames{sf}).data_frm_mean_fly = data;
        perexpphasefeatures.(fieldname).(subfieldnames{sf}).mean_frm_mean_fly = circ_mean(data);
        perexpphasefeatures.(fieldname).(subfieldnames{sf}).n_frm_mean_fly = numel(data);
        perexpphasefeatures.(fieldname).(subfieldnames{sf}).std_frm_mean_fly = circ_std(data);


        % mean_walk_mean_fly
        curr_all = {};
        data =[];
        for f = 1:numel(perfly_metrics)
            curr_all{f} = [perfly_metrics(f).(fieldname).(subfieldnames{sf}).walk_mean_fly];
        end
        curr_all = horzcat(curr_all{:});
        data = curr_all(~isnan(curr_all));
        perexpphasefeatures.(fieldname).(subfieldnames{sf}).data_walk_mean_fly = data;
        perexpphasefeatures.(fieldname).(subfieldnames{sf}).mean_walk_mean_fly = circ_mean(data);
        perexpphasefeatures.(fieldname).(subfieldnames{sf}).n_walk_mean_fly = numel(data);
        perexpphasefeatures.(fieldname).(subfieldnames{sf}).std_walk_mean_fly = circ_std(data);

    end

end




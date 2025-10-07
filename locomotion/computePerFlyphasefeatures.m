function [perflyphasefeatures] = computePerFlyphasefeatures(walkfeaturestruct,fieldname)

subfieldnames = fields(walkfeaturestruct(1).(fieldname));
% subfieldnames = {'data','mean'};
for sf = 1:numel(subfieldnames)

    if strcmp(subfieldnames{sf},'phasedata')
        curr_all = {};
        for w = 1:numel(walkfeaturestruct)
            curr_all{w} = [walkfeaturestruct(w).(fieldname).(subfieldnames{sf})];
        end
        perflyphasefeatures.(fieldname).(subfieldnames{sf}) = horzcat(curr_all{:});
    else
        
        % frm_data
        curr_frm_all = {};
        for w = 1:numel(walkfeaturestruct)
            curr_frm_all{w} = [walkfeaturestruct(w).(fieldname).(subfieldnames{sf}).data];
        end
        currfrmdata = horzcat(curr_frm_all{:});
        currdata = currfrmdata(~isnan(currfrmdata));
        
        perflyphasefeatures.(fieldname).(subfieldnames{sf}).frm_data_fly = currdata;
        perflyphasefeatures.(fieldname).(subfieldnames{sf}).frm_mean_fly = circ_mean(currdata);        
        perflyphasefeatures.(fieldname).(subfieldnames{sf}).frm_n_fly = numel(currdata);
        perflyphasefeatures.(fieldname).(subfieldnames{sf}).frm_std_fly = circ_std(currdata);

        % walk_data

        curr_walk = {};
        for w = 1:numel(walkfeaturestruct)
            curr_walk{w} = [walkfeaturestruct(w).(fieldname).(subfieldnames{sf}).mean];
        end
        currwalkdata = horzcat(curr_walk{:});
        currdatawalks = currwalkdata(~isnan(currwalkdata));

        perflyphasefeatures.(fieldname).(subfieldnames{sf}).walk_data_fly = currdatawalks;
        perflyphasefeatures.(fieldname).(subfieldnames{sf}).walk_mean_fly = circ_mean(currdatawalks);        
        perflyphasefeatures.(fieldname).(subfieldnames{sf}).walk_n_fly = numel(currdatawalks);
        perflyphasefeatures.(fieldname).(subfieldnames{sf}).walk_std_fly = circ_std(currdatawalks);

    end
end     
        


end




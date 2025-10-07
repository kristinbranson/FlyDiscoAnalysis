function [perflyperframefeatures] = computePerFlywalkperframefeatures(walkfeaturestruct,fieldname)

% subfieldnames = fields(walkfeaturestruct(1).(fieldname));
% % subfieldname = {'data','mean'};
% for sf = 1:numel(subfieldnames)
%     curr_all = {};
%     for w = 1:numel(walkfeaturestruct)
%         curr_all{w} = [walkfeaturestruct(w).(fieldname).(subfieldnames{sf})];
%     end
%     if strcmp(subfieldnames{sf},'data')
%         perflyperframefeatures.(fieldname).(subfieldnames{sf}) = horzcat(curr_all{:});
%     elseif strcmp(subfieldnames{sf},'n')
%         perflyperframefeatures.(fieldname).('n_walks') = sum(numel(curr_all));
%     elseif strcmp(subfieldnames{sf},'pffname')
%         perflyperframefeatures.(fieldname).(subfieldnames{sf}) = walkfeaturestruct(1).(fieldname).(subfieldnames{sf});
%     else
%         perflyperframefeatures.(fieldname).(['mean_',subfieldnames{sf}]) = mean(horzcat(curr_all{:}),'omitnan');     
% 
%     end
% end
% perflyperframefeatures.(fieldname).mean = mean(perflyperframefeatures.(fieldname).data,'omitnan');
% perflyperframefeatures.(fieldname).std = std(perflyperframefeatures.(fieldname).data,'omitnan');
% perflyperframefeatures.(fieldname).n = numel(perflyperframefeatures.(fieldname).data);


% frm_data = data from all walking frames
curr_all = {};
for w = 1:numel(walkfeaturestruct)
    curr_all{w} = [walkfeaturestruct(w).(fieldname).data];
end
currdata = horzcat(curr_all{:});
perflyperframefeatures.(fieldname).frm_data_fly = currdata(~isnan(currdata));
perflyperframefeatures.(fieldname).frm_mean_fly = mean(currdata,'omitnan');
perflyperframefeatures.(fieldname).frm_n_fly = numel(currdata(~isnan(currdata)));
perflyperframefeatures.(fieldname).frm_std_fly = std(currdata,'omitnan');

% walk_data = average over walks for this fly
curr_means = [];
for w = 1:numel(walkfeaturestruct)
    curr_means{w} = [walkfeaturestruct(w).(fieldname).mean];
end
currmeandata = horzcat(curr_means{:});
perflyperframefeatures.(fieldname).walk_data_fly = currdata(~isnan(currmeandata));
perflyperframefeatures.(fieldname).walk_mean_fly = mean(currmeandata,'omitnan');
perflyperframefeatures.(fieldname).walk_n_fly = numel(currmeandata(~isnan(currmeandata)));
perflyperframefeatures.(fieldname).walk_std_fly = std(currmeandata,'omitnan');


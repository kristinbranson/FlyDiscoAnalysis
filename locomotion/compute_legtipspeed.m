function [mean_speed,std_speed, min_speed, max_speed] = compute_legtipspeed(fly,position_data,limb_idx,start_indices,end_indices,currfly_timestamps,pxpermm)


mean_speed = nan(1,numel(start_indices));
std_speed = nan(1,numel(start_indices));
min_speed = nan(1,numel(start_indices));
max_speed = nan(1,numel(start_indices));
if ~isempty(start_indices)
    for i = 1:numel(start_indices)

        t = currfly_timestamps(start_indices(i):end_indices(i));
        dt = diff(t);

        x = squeeze(position_data{fly}(limb_idx,1,start_indices(i):end_indices(i)))';
        y = squeeze(position_data{fly}(limb_idx,2,start_indices(i):end_indices(i)))';
        x_mm = x./pxpermm;
        y_mm = y./pxpermm;
        % diff
        dx_mm = diff(x_mm);
        dy_mm = diff(y_mm);

        speed_during_bout = sqrt(dx_mm.^2+dy_mm.^2)./dt;
        mean_speed(i) = mean(speed_during_bout);
        std_speed(i) = std(speed_during_bout);
        min_speed(i) = min(speed_during_bout);
        max_speed(i) = max(speed_during_bout);
    end
end
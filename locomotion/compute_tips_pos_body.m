function [tips_pos_body] = compute_tips_pos_body(aptdata,apt_pts_4_center,apt_pt_4_theta,legtip_landmarknums)

[body_aligned_pts] = convert_global2body(aptdata,apt_pts_4_center,apt_pt_4_theta);
tips_pos_body = cell(1,numel(aptdata));
for fly = 1:numel(aptdata)
  tips_pos_body{fly}  = body_aligned_pts{fly}(legtip_landmarknums,:,:);
end

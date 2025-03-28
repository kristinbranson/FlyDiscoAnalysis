function [body_aligned_pts] = convert_global2body(pts,ctr_KPnums, theta_KPnum)
% 
% input
% pts = double for 1 fly or cell for all flies in a movie
% aptdata for keypoints to be centered and rotated LD x d X T (21 x 2 x T landmarks)

% ctrKPs = which keypoints to use to compute centroid (keypoints 4 and 5, anterior thorax)
% thetarefKP = key point to use as reference angle for theta (keypoint 6, notum)
% output
% either double or cell per fly 21 x 2 x T

if iscell(pts)
    body_aligned_pts = cell(1,numel(pts));
    for fly = 1:numel(pts)
        curr_pts = pts{fly};
        body_aligned_pts{fly} = convert_global2body_perfly(curr_pts,ctr_KPnums, theta_KPnum);
    end

else
    body_aligned_pts = convert_global2body_perfly(pts,ctr_KPnums, theta_KPnum);

end

end


% code runs on single fly
function [body_aligned_pts_perfly] = convert_global2body_perfly(pts,ctr_KPnums, theta_KPnum)

[~,d,T] = size(pts);
assert(d==2);

% compute fly center
ctrx = mean(pts(ctr_KPnums,1,:),1);
ctry = mean(pts(ctr_KPnums,2,:),1);

% substract off center
% assert(d==2);
mu = cat(1,ctrx,ctry);
pts2 = permute(pts,[2,1,3]);
pts2 = pts2-mu;

% compute fly theta
cx = pts2(1,theta_KPnum,:);
cy = pts2(2,theta_KPnum,:);

theta = atan2(cy,cx);

% negative because theta should point to head + bias to get head to point
% up
b = pi/2 * ones(1,1,T);
costheta = -(cos(theta-b));
sintheta = -(sin(theta-b)) ;


% align to postive y axis (in image coordinates = negative y axis)
body_aligned_pts_perfly = [costheta.*pts2(1,:,:) + sintheta.*pts2(2,:,:)
    -sintheta.*pts2(1,:,:) + costheta.*pts2(2,:,:)];
body_aligned_pts_perfly = permute(body_aligned_pts_perfly,[2,1,3]);

end

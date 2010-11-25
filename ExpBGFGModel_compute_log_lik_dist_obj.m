function [log_lik_dist_obj_given_fore,log_lik_dist_obj_given_back] = ExpBGFGModel_compute_log_lik_dist_obj(im,model,c0,c1)
% ExpBGFGModel_compute_log_lik_dist_obj()
        
% computes the log likelihood of the distance to the nearest object detection
% for each pixel in the image self.im given that the pixel is in foreground and
% given that the pixel is in background.

nc = size(im,2);

if nargin < 3,
  c0 = 1;
end
if nargin < 4,
  c1 = nc;
end
        
% detect objects in the current image
isobj = obj_detect(c0,c1);

% compute distances to detections
obj_detection_dist = bwdist(isobj);

% find which bin each distance falls in
[~,idx] = hist(obj_detection_dist,model.obj_detection_dist_centers);

% lookup probabilities for each bin idx
log_lik_dist_obj_given_fore = log(model.obj_detection_dist_frac_fg(idx));
log_lik_dist_obj_given_back = log(model.obj_detection_dist_frac_bg(idx));

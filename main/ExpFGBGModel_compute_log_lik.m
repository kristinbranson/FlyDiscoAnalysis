function [log_lik_given_fore,log_lik_given_back,...
  log_lik_appearance_given_fore,log_lik_appearance_given_back,...
  log_lik_dist_obj_given_fore,log_lik_dist_obj_given_back] = ...
  ExpFGBGModel_compute_log_lik(im,model,c0,c1)

nc = size(im,2);

if nargin < 3,
  c0 = 1;
end
if nargin < 4,
  c1 = nc;
end

log_lik_appearance_given_fore = compute_log_lik_appearance_given_fore(im,model,c0,c1);
log_lik_appearance_given_back = compute_log_lik_appearance_given_back(im,model,c0,c1);
[log_lik_dist_obj_given_fore,log_lik_dist_obj_given_back] = ...
  compute_log_lik_dist_obj(im,model,c0,c1);
log_lik_given_fore = log_lik_appearance_given_fore + log_lik_dist_obj_given_fore;
log_lik_given_back =log_lik_appearance_given_back + log_lik_dist_obj_given_back;
function log_lik_appearance_given_back = ExpBGFGModel_compute_log_lik_appearance_given_back(im,model,c0,c1)

% log_lik_appearance_given_back = compute_log_lik_appearance_given_back(im,model,c0=1,c1=nc)
% computes the log likelihood of each pixel appearance in
% im given that the pixel is foreground.

nc = size(im,2);

if nargin < 3,
  c0 = 1;
end
if nargin < 4,
  c1 = nc;
end

d2 = (im(:,c0:c1,:) - model.bg_mu(:,c0:c1,:)).^2 ./ (2*model.bg_sigma(:,c0:c1,:).^2);
log_lik_appearance_given_back = -d2 - model.bg_log_Z(:,c0:c1,:);

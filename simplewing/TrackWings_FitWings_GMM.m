function wingtrx = TrackWings_FitWings_GMM(dthetawing,fore2dbkgd,wingtrxprev,params)

gmmweights = max(0,fore2dbkgd(idxcurr)-params.mindwing_high);

% or try gmm-based approach
aic = nan(1,2);
%bic = nan(1,2);
mu = cell(1,2);
S = cell(1,2);
priors = cell(1,2);

k = 1;
[mu{k},S{k}] = weighted_mean_cov(dthetawing,gmmweights);
priors{k} = 1; %#ok<NASGU>
nllcurr = -sum(log(normpdf(dthetawing,mu{k},sqrt(S{k}))).*gmmweights);
aic(k) = 2*nllcurr+2*3*k;
%bic(k) = 2*nllcurr+3*k*log(nwingpx);

k = 2;
if isempty(wingtrxprev),
  [mucurr,Scurr,priorscurr,~,nllcurr] = mygmm(dthetawing,k,'Replicates',10,'weights',gmmweights); %#ok<ASGLU>
else
  startmu = [wingtrxprev.wing_anglel;wingtrxprev.wing_angler];
  dmu = abs(startmu(1)-startmu(2));
  if dmu < .02;
    startmu(1) = startmu(1) - (.01-dmu);
    startmu(2) = startmu(2) + (.01-dmu);
  end
  [mucurr,Scurr,priorscurr,~,nllcurr] = mygmm(dthetawing,k,'Start',startmu,'weights',gmmweights); %#ok<ASGLU>
end
mu{k} = mucurr; %S{k} = Scurr; priors{k} = priorscurr;
aic(k) = 2*nllcurr+2*3*k;
%bic(k) = 2*nllcurr+3*k*log(sum(gmmweights));

%       if min(priors{2}) < wing_min_prior,
%         k = 1;
%       else
%         [~,k] = min(bic);
%       end
[~,k] = min(aic);
wing_angles_curr = mu{k};
if numel(wing_angles_curr) == 1,
  wing_angles_curr = repmat(wing_angles_curr,[1,2]);
end

wingtrx = struct('wing_anglel',[],...
  'wing_angler',[],...
  'nwingsdetected',[],...
  'wing_areal',[],...
  'wing_arear',[],...
  'wing_trough_angle',[]);

wingtrx.wing_anglel = wing_angles_curr(1);
wingtrx.wing_angler = wing_angles_curr(2);
% TODO: compute these other stats
% wingtrx.nwingsdetected = npeaks;
% wingtrx.wing_areal = area(1);
% wingtrx.wing_arear = area(2);
% wingtrx.wing_trough_angle = wing_trough_angle;

  
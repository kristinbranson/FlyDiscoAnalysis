function [bias_diagnostics,res] = BowlBiasDiagnostics(expdir,varargin)

bias_diagnostics = struct;

%% parse parameters
[analysis_protocol,settingsdir,datalocparamsfilestr] = ...
  myparse(varargin,...
  'analysis_protocol','current',...
  'settingsdir','/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings',...
  'datalocparamsfilestr','dataloc_params.txt');

%% read experiment trx

trx = Trx('analysis_protocol',analysis_protocol,'settingsdir',settingsdir,...
  'datalocparamsfilestr',datalocparamsfilestr);

fprintf('Loading trajectories for %s...\n',expdir);

trx.AddExpDir(expdir);

%% read parameters

biasdiagnosticsparamsfile = trx.getDataLoc('biasdiagnosticsparamsfilestr');
params = ReadParams(biasdiagnosticsparamsfile);

%% set up figure

res = struct;
hfig = 1;
figure(hfig);
clf(hfig,'reset');
set(hfig,'Units','Pixels','Position',params.figpos);
nax_r = 4;
nax_c = 4;
hax = createsubplots(nax_r,nax_c,.04,hfig);
hax = reshape(hax(1:nax_r*nax_c),[nax_r,nax_c]);
res.hax = hax;
res.hfig = hfig;

%% collect data
r = [];
theta = [];
for fly = 1:trx.nflies,
  ismoving = trx(fly).velmag >= params.minvelmag;
  r = cat(2,r,trx(fly).arena_r(ismoving));
  theta = cat(2,theta,trx(fly).arena_angle(ismoving));
end
  
%% histogram

[edges_theta,centers_theta] = SelectHistEdges(params.nbins_theta,[-pi,pi],'linear');
[edges_r,centers_r] = SelectHistEdges(params.nbins_r,[0,1]*trx.landmark_params.arena_radius_mm,'linear');
% rows = r, cols = theta
frac = hist3([r;theta]','Edges',{edges_r,edges_theta});
frac = [frac(:,1:end-2),frac(:,end-1)+frac(:,end)];
frac = [frac(1:end-2,:);frac(end-1,:)+frac(end,:)];
frac = frac / (sum(frac(:)));

%% convolve over-sampled histogram with Gaussian
% each bin has length arena_radius_mm/nbins_r, 
% 2*pi/nbins_theta
nstd = 3;
binsize = [1/params.nbins_r,2*pi/params.nbins_theta];
hsize = ceil(nstd*[params.sigma_r,params.sigma_theta]./binsize);
sigma = [params.sigma_r,params.sigma_theta];
i = 1;
f = normpdf((-hsize(i):hsize(i))*binsize(i),0,sigma(i));
f = f / sum(f);
fracsmooth = imfilter(frac,f',0,'same','corr');
i = 2;
f = normpdf((-hsize(i):hsize(i))*binsize(i),0,sigma(i));
f = f / sum(f);
fracsmooth = imfilter(fracsmooth,f,0,'same','corr');

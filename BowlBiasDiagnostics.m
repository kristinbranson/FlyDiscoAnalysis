function [bias_diagnostics,hfig] = BowlBiasDiagnostics(expdir,varargin)

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

%% collect data
r = [];
theta = [];
for fly = 1:trx.nflies,
  ismoving = trx(fly).velmag >= params.minvelmag;
  r = cat(2,r,trx(fly).arena_r(ismoving));
  theta = cat(2,theta,trx(fly).arena_angle(ismoving));
end
% all frames
rall = cell2mat(trx.arena_r);
thetaall = cell2mat(trx.arena_angle);
  
%% histogram

[edges_theta,centers_theta] = SelectHistEdges(params.nbins_theta,[-pi,pi],'linear');
[edges_r,centers_r] = SelectHistEdges(params.nbins_r,[0,1]*trx.landmark_params.arena_radius_mm,'linear');
% rows = r, cols = theta
frac = hist3([r;theta]','Edges',{edges_r,edges_theta});
frac = [frac(:,1:end-2),frac(:,end-1)+frac(:,end)];
frac = [frac(1:end-2,:);frac(end-1,:)+frac(end,:)];
% all frames
fracall = hist3([rall;thetaall]','Edges',{edges_r,edges_theta});
fracall = [fracall(:,1:end-2),fracall(:,end-1)+fracall(:,end)];
fracall = [fracall(1:end-2,:);fracall(end-1,:)+fracall(end,:)];


%% normalize by bin size
binareaz = diff(edges_r.^2);
frac = bsxfun(@rdivide,frac,binareaz');
frac = frac / (sum(frac(:)));

fracall = bsxfun(@rdivide,fracall,binareaz');
fracall = fracall / (sum(fracall(:)));

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
fracallsmooth = imfilter(fracall,f',0,'same','corr');
i = 2;
f = normpdf((-hsize(i):hsize(i))*binsize(i),0,sigma(i));
f = f / sum(f);
fracsmooth = imfilter(fracsmooth,f,'circular','same','corr');
fracallsmooth = imfilter(fracallsmooth,f,'circular','same','corr');

bias_diagnostics.frac = frac;
bias_diagnostics.fracsmooth = fracsmooth;
bias_diagnostics.fracall = fracall;
bias_diagnostics.fracallsmooth = fracallsmooth;

%% maximum ratio between frac for different angles
rmin_bin = min(max(1,floor(params.min_r/binsize(1))),params.nbins_r);
sumfracsmooth = sum(fracsmooth(rmin_bin:end,:),1);
sumfracsmooth = sumfracsmooth / sum(sumfracsmooth);
[max_sumfracsmooth,max_angle_bin] = max(sumfracsmooth);
argmaxangle_sumfracsmooth = centers_theta(max_angle_bin);
[min_sumfracsmooth,min_angle_bin] = min(sumfracsmooth);
argminangle_sumfracsmooth = centers_theta(min_angle_bin);
maxratio_sumfracsmooth = max_sumfracsmooth./min_sumfracsmooth;
diffargextremeangle_sumfracsmooth = abs(modrange(argmaxangle_sumfracsmooth - argminangle_sumfracsmooth,-pi,pi));

bias_diagnostics.sumfracsmooth = sumfracsmooth;
bias_diagnostics.max_sumfracsmooth = max_sumfracsmooth;
bias_diagnostics.min_sumfracsmooth = min_sumfracsmooth;
bias_diagnostics.argmaxangle_sumfracsmooth = argmaxangle_sumfracsmooth;
bias_diagnostics.argminangle_sumfracsmooth = argminangle_sumfracsmooth;
bias_diagnostics.maxratio_sumfracsmooth = maxratio_sumfracsmooth;
bias_diagnostics.diffargextremaangle_sumfracsmooth = diffargextremeangle_sumfracsmooth;

% all frames
rmin_bin = min(max(1,floor(params.min_r/binsize(1))),params.nbins_r);
sumfracallsmooth = sum(fracallsmooth(rmin_bin:end,:),1);
sumfracallsmooth = sumfracallsmooth / sum(sumfracallsmooth);
[max_sumfracallsmooth,max_angle_bin] = max(sumfracallsmooth);
argmaxangle_sumfracallsmooth = centers_theta(max_angle_bin);
[min_sumfracallsmooth,min_angle_bin] = min(sumfracallsmooth);
argminangle_sumfracallsmooth = centers_theta(min_angle_bin);
maxratio_sumfracallsmooth = max_sumfracallsmooth./min_sumfracallsmooth;
diffargextremeangle_sumfracallsmooth = abs(modrange(argmaxangle_sumfracallsmooth - argminangle_sumfracallsmooth,-pi,pi));

bias_diagnostics.sumfracallsmooth = sumfracallsmooth;
bias_diagnostics.max_sumfracallsmooth = max_sumfracallsmooth;
bias_diagnostics.min_sumfracallsmooth = min_sumfracallsmooth;
bias_diagnostics.argmaxangle_sumfracallsmooth = argmaxangle_sumfracallsmooth;
bias_diagnostics.argminangle_sumfracallsmooth = argminangle_sumfracallsmooth;
bias_diagnostics.maxratio_sumfracallsmooth = maxratio_sumfracallsmooth;
bias_diagnostics.diffargextremaangle_sumfracallsmooth = diffargextremeangle_sumfracallsmooth;


%% maximum difference in mode in radius
[maxfracsmooth,maxr_bin] = max(fracsmooth,[],1);
argmaxr_fracsmooth = centers_r(maxr_bin);
[max_maxfracsmooth,max_maxr_angle_bin] = max(argmaxr_fracsmooth);
argmaxangle_maxfracsmooth = centers_theta(max_maxr_angle_bin);
[min_maxfracsmooth,min_maxr_angle_bin] = min(argmaxr_fracsmooth);
argminangle_maxfracsmooth = centers_theta(min_maxr_angle_bin);
maxdiff_maxfracsmooth = (max_maxfracsmooth - min_maxfracsmooth) / trx.landmark_params.arena_radius_mm;
diffargextremaangle_maxfracsmooth = abs(modrange(argmaxangle_maxfracsmooth - argminangle_maxfracsmooth,-pi,pi));

bias_diagnostics.maxfracsmooth = maxfracsmooth;
bias_diagnostics.max_maxfracsmooth = max_maxfracsmooth;
bias_diagnostics.min_maxfracsmooth = min_maxfracsmooth;
bias_diagnostics.argmaxangle_maxfracsmooth = argmaxangle_maxfracsmooth;
bias_diagnostics.argminangle_minfracsmooth = argminangle_maxfracsmooth;
bias_diagnostics.maxdiff_maxfracsmooth = maxdiff_maxfracsmooth;
bias_diagnostics.diffargextremaangle_maxfracsmooth = diffargextremaangle_maxfracsmooth;
bias_diagnostics.argmaxr_fracsmooth = argmaxr_fracsmooth;

% all frames
[maxfracallsmooth,maxr_bin] = max(fracallsmooth,[],1);
argmaxr_fracallsmooth = centers_r(maxr_bin);
[max_maxfracallsmooth,max_maxr_angle_bin] = max(argmaxr_fracallsmooth);
argmaxangle_maxfracallsmooth = centers_theta(max_maxr_angle_bin);
[min_maxfracallsmooth,min_maxr_angle_bin] = min(argmaxr_fracallsmooth);
argminangle_maxfracallsmooth = centers_theta(min_maxr_angle_bin);
maxdiff_maxfracallsmooth = (max_maxfracallsmooth - min_maxfracallsmooth) / trx.landmark_params.arena_radius_mm;
diffargextremaangle_maxfracallsmooth = abs(modrange(argmaxangle_maxfracallsmooth - argminangle_maxfracallsmooth,-pi,pi));

bias_diagnostics.maxfracallsmooth = maxfracallsmooth;
bias_diagnostics.max_maxfracallsmooth = max_maxfracallsmooth;
bias_diagnostics.min_maxfracallsmooth = min_maxfracallsmooth;
bias_diagnostics.argmaxangle_maxfracallsmooth = argmaxangle_maxfracallsmooth;
bias_diagnostics.argminangle_minfracallsmooth = argminangle_maxfracallsmooth;
bias_diagnostics.maxdiff_maxfracallsmooth = maxdiff_maxfracallsmooth;
bias_diagnostics.diffargextremaangle_maxfracallsmooth = diffargextremaangle_maxfracallsmooth;
bias_diagnostics.argmaxr_fracallsmooth = argmaxr_fracallsmooth;

%% set up figure

hfig = 1;
figure(hfig);
clf(hfig,'reset');
set(hfig,'Units','Pixels','Position',params.figpos);
nax_r = 2;
nax_c = 4;
hax = createsubplots(nax_r,nax_c,[.025,.05],hfig);
hax = reshape(hax(1:nax_r*nax_c),[nax_r,nax_c]);

%% plot 2-D histogram

gray = [.5,.5,.5];

axi = 1;
polarimagesc(edges_r,edges_theta,[0,0],frac,'Parent',hax(axi));
hold(hax(axi),'on');
plot(hax(axi),[0,cos(argmaxangle_sumfracsmooth)*trx.landmark_params.arena_radius_mm],...
  [0,sin(argmaxangle_sumfracsmooth)*trx.landmark_params.arena_radius_mm],...
  'w-','linewidth',2);
plot(hax(axi),[0,cos(argmaxangle_sumfracsmooth)*trx.landmark_params.arena_radius_mm],...
  [0,sin(argmaxangle_sumfracsmooth)*trx.landmark_params.arena_radius_mm],...
  '--','color',gray,'linewidth',2);
plot(hax(axi),[0,cos(argminangle_sumfracsmooth)*trx.landmark_params.arena_radius_mm],...
  [0,sin(argminangle_sumfracsmooth)*trx.landmark_params.arena_radius_mm],...
  'w-','linewidth',2);
plot(hax(axi),[0,cos(argminangle_sumfracsmooth)*trx.landmark_params.arena_radius_mm],...
  [0,sin(argminangle_sumfracsmooth)*trx.landmark_params.arena_radius_mm],...
  '--','color',gray,'linewidth',2);
axis(hax(axi),'xy','equal',trx.landmark_params.arena_radius_mm*[-1,1,-1,1]);
%set(hax(axi),'Color',black);
title(hax(axi),'Position heatmap, moving flies');
hcb = colorbar('peer',hax(axi),'West');
cbpos = get(hcb,'Position');
cbpos(1) = .01;
set(hcb,'Position',cbpos);

axi = 2;
polarimagesc(edges_r,edges_theta,[0,0],fracsmooth,'Parent',hax(axi));
hold(hax(axi),'on');
plot(hax(axi),[0,cos(argmaxangle_sumfracsmooth)*trx.landmark_params.arena_radius_mm],...
  [0,sin(argmaxangle_sumfracsmooth)*trx.landmark_params.arena_radius_mm],...
  'w-','linewidth',2);
plot(hax(axi),[0,cos(argmaxangle_sumfracsmooth)*trx.landmark_params.arena_radius_mm],...
  [0,sin(argmaxangle_sumfracsmooth)*trx.landmark_params.arena_radius_mm],...
  '--','color',gray,'linewidth',2);
plot(hax(axi),[0,cos(argminangle_sumfracsmooth)*trx.landmark_params.arena_radius_mm],...
  [0,sin(argminangle_sumfracsmooth)*trx.landmark_params.arena_radius_mm],...
  'w-','linewidth',2);
plot(hax(axi),[0,cos(argminangle_sumfracsmooth)*trx.landmark_params.arena_radius_mm],...
  [0,sin(argminangle_sumfracsmooth)*trx.landmark_params.arena_radius_mm],...
  '--','color',gray,'linewidth',2);
axis(hax(axi),'xy','equal',trx.landmark_params.arena_radius_mm*[-1,1,-1,1]);
%set(hax(axi),'Color',black);
title(hax(axi),'Smoothed position heatmap, moving flies');
linkaxes(hax(1:2));

axi = 3;
plot(hax(axi),centers_theta,sumfracsmooth,'k-','linewidth',2);
maxy = max(max_sumfracallsmooth,max_sumfracsmooth);
miny = min(min_sumfracallsmooth,min_sumfracsmooth);
dy = (maxy - miny)*.01;
ax = [edges_theta(1),edges_theta(end),miny-dy,maxy+dy];
hold(hax(axi),'on');
plot(hax(axi),argmaxangle_sumfracsmooth+[0,0],ax(3:4),'k-','linewidth',2);
plot(hax(axi),argmaxangle_sumfracsmooth+[0,0],ax(3:4),'--','color',gray,...
  'linewidth',2);
plot(hax(axi),argminangle_sumfracsmooth+[0,0],ax(3:4),'k-','linewidth',2);
plot(hax(axi),argminangle_sumfracsmooth+[0,0],ax(3:4),'--','color',gray,...
  'linewidth',2);
axis(hax(axi),ax);
%set(hax(axi),'Color',black);
ylabel(hax(axi),'Fraction of time | angle');

axi = 4;
imagesc([centers_theta(1),centers_theta(end)],[centers_r(1),centers_r(end)],fracsmooth,'parent',hax(axi));
hold(hax(axi),'on');
plot(hax(axi),centers_theta,argmaxr_fracsmooth,'w-','linewidth',2);
plot(hax(axi),centers_theta,argmaxr_fracsmooth,'--','color',gray,'linewidth',2);
plot(hax(axi),argmaxangle_sumfracsmooth+[0,0],[edges_r(1),edges_r(end)],'w-','linewidth',2);
plot(hax(axi),argmaxangle_sumfracsmooth+[0,0],[edges_r(1),edges_r(end)],'--','color',gray,'linewidth',2);
plot(hax(axi),argminangle_sumfracsmooth+[0,0],[edges_r(1),edges_r(end)],'w-','linewidth',2);
plot(hax(axi),argminangle_sumfracsmooth+[0,0],[edges_r(1),edges_r(end)],'--','color',gray,'linewidth',2);
plot(hax(axi),[argmaxangle_maxfracsmooth,argminangle_maxfracsmooth],[max_maxfracsmooth,min_maxfracsmooth],'o','color',gray,'markerfacecolor','w');
% text(argmaxangle_sumfracsmooth,centers_r(1),...
%   sprintf('frac(%.2f) = %.3f',argmaxangle_sumfracsmooth,max_sumfracsmooth),...
%   'HorizontalAlignment','right','VerticalAlignment','bottom',...
%   'Color','w','Parent',hax(axi),...
%   'Rotation',90);
% text(argminangle_sumfracsmooth,centers_r(1),...
%   sprintf('frac(%.2f) = %.3f',argminangle_sumfracsmooth,min_sumfracsmooth),...
%   'HorizontalAlignment','right','VerticalAlignment','bottom',...
%   'Color','w','Parent',hax(axi),...
%   'Rotation',90);
axis(hax(axi),[edges_theta(1),edges_theta(end),edges_r(1),edges_r(end)]);
xlabel(hax(axi),'Angle (rad)');
ylabel(hax(axi),'Distance from center (mm)');
title(hax(axi),'Unwrapped, smoothed position heatmap, moving flies');

%% same thing, but all frames

axi = 5;
polarimagesc(edges_r,edges_theta,[0,0],fracall,'Parent',hax(axi));
hold(hax(axi),'on');
plot(hax(axi),[0,cos(argmaxangle_sumfracallsmooth)*trx.landmark_params.arena_radius_mm],...
  [0,sin(argmaxangle_sumfracallsmooth)*trx.landmark_params.arena_radius_mm],...
  'w-','linewidth',2);
plot(hax(axi),[0,cos(argmaxangle_sumfracallsmooth)*trx.landmark_params.arena_radius_mm],...
  [0,sin(argmaxangle_sumfracallsmooth)*trx.landmark_params.arena_radius_mm],...
  '--','color',gray,'linewidth',2);
plot(hax(axi),[0,cos(argminangle_sumfracallsmooth)*trx.landmark_params.arena_radius_mm],...
  [0,sin(argminangle_sumfracallsmooth)*trx.landmark_params.arena_radius_mm],...
  'w-','linewidth',2);
plot(hax(axi),[0,cos(argminangle_sumfracallsmooth)*trx.landmark_params.arena_radius_mm],...
  [0,sin(argminangle_sumfracallsmooth)*trx.landmark_params.arena_radius_mm],...
  '--','color',gray,'linewidth',2);
axis(hax(axi),'xy','equal',trx.landmark_params.arena_radius_mm*[-1,1,-1,1]);
%set(hax(axi),'Color',black);
title(hax(axi),'Position heatmap, all frames');
%hcb = colorbar('peer',hax(axi),'West');
%cbpos = get(hcb,'Position');
%cbpos(1) = .01;
%set(hcb,'Position',cbpos);

axi = 6;
polarimagesc(edges_r,edges_theta,[0,0],fracallsmooth,'Parent',hax(axi));
hold(hax(axi),'on');
plot(hax(axi),[0,cos(argmaxangle_sumfracallsmooth)*trx.landmark_params.arena_radius_mm],...
  [0,sin(argmaxangle_sumfracallsmooth)*trx.landmark_params.arena_radius_mm],...
  'w-','linewidth',2);
plot(hax(axi),[0,cos(argmaxangle_sumfracallsmooth)*trx.landmark_params.arena_radius_mm],...
  [0,sin(argmaxangle_sumfracallsmooth)*trx.landmark_params.arena_radius_mm],...
  '--','color',gray,'linewidth',2);
plot(hax(axi),[0,cos(argminangle_sumfracallsmooth)*trx.landmark_params.arena_radius_mm],...
  [0,sin(argminangle_sumfracallsmooth)*trx.landmark_params.arena_radius_mm],...
  'w-','linewidth',2);
plot(hax(axi),[0,cos(argminangle_sumfracallsmooth)*trx.landmark_params.arena_radius_mm],...
  [0,sin(argminangle_sumfracallsmooth)*trx.landmark_params.arena_radius_mm],...
  '--','color',gray,'linewidth',2);
axis(hax(axi),'xy','equal',trx.landmark_params.arena_radius_mm*[-1,1,-1,1]);
%set(hax(axi),'Color',black);
title(hax(axi),'Smoothed position heatmap, all frames');
linkaxes(hax(1:2));

axi = 7;
plot(hax(axi),centers_theta,sumfracallsmooth,'k-','linewidth',2);
hold(hax(axi),'on');
plot(hax(axi),argmaxangle_sumfracallsmooth+[0,0],ax(3:4),'k-','linewidth',2);
plot(hax(axi),argmaxangle_sumfracallsmooth+[0,0],ax(3:4),'--','color',gray,...
  'linewidth',2);
plot(hax(axi),argminangle_sumfracallsmooth+[0,0],ax(3:4),'k-','linewidth',2);
plot(hax(axi),argminangle_sumfracallsmooth+[0,0],ax(3:4),'--','color',gray,...
  'linewidth',2);
axis(hax(axi),ax);
%set(hax(axi),'Color',black);
ylabel(hax(axi),'Fraction of time | angle');

axi = 8;
imagesc([centers_theta(1),centers_theta(end)],[centers_r(1),centers_r(end)],fracallsmooth,'parent',hax(axi));
hold(hax(axi),'on');
plot(hax(axi),centers_theta,argmaxr_fracallsmooth,'w-','linewidth',2);
plot(hax(axi),centers_theta,argmaxr_fracallsmooth,'--','color',gray,'linewidth',2);
plot(hax(axi),argmaxangle_sumfracallsmooth+[0,0],[edges_r(1),edges_r(end)],'w-','linewidth',2);
plot(hax(axi),argmaxangle_sumfracallsmooth+[0,0],[edges_r(1),edges_r(end)],'--','color',gray,'linewidth',2);
plot(hax(axi),argminangle_sumfracallsmooth+[0,0],[edges_r(1),edges_r(end)],'w-','linewidth',2);
plot(hax(axi),argminangle_sumfracallsmooth+[0,0],[edges_r(1),edges_r(end)],'--','color',gray,'linewidth',2);
plot(hax(axi),[argmaxangle_maxfracallsmooth,argminangle_maxfracallsmooth],[max_maxfracallsmooth,min_maxfracallsmooth],'o','color',gray,'markerfacecolor','w');
% text(argmaxangle_sumfracallsmooth,centers_r(1),...
%   sprintf('fracall(%.2f) = %.3f',argmaxangle_sumfracallsmooth,max_sumfracallsmooth),...
%   'HorizontalAlignment','right','VerticalAlignment','bottom',...
%   'Color','w','Parent',hax(axi),...
%   'Rotation',90);
% text(argminangle_sumfracallsmooth,centers_r(1),...
%   sprintf('fracall(%.2f) = %.3f',argminangle_sumfracallsmooth,min_sumfracallsmooth),...
%   'HorizontalAlignment','right','VerticalAlignment','bottom',...
%   'Color','w','Parent',hax(axi),...
%   'Rotation',90);
axis(hax(axi),[edges_theta(1),edges_theta(end),edges_r(1),edges_r(end)]);
xlabel(hax(axi),'Angle (rad)');
ylabel(hax(axi),'Distance from center (mm)');
title(hax(axi),'Unwrapped, smoothed position heatmap, all frames');

linkaxes(hax([3,7]));
linkaxes(hax([3,4]),'x');
linkaxes(hax([7,8]),'x');

%% save image

savename = fullfile(expdir,trx.dataloc_params.biasdiagnosticsimagefilestr);
if exist(savename,'file'),
  delete(savename);
end
save2png(savename,hfig);
  
%% save to mat file
biasdiagnosticsmatfilename = fullfile(expdir,trx.dataloc_params.biasdiagnosticsmatfilestr);
save(biasdiagnosticsmatfilename,'-struct','bias_diagnostics');

%% write to text file

biasdiagnosticsfilename = fullfile(expdir,trx.dataloc_params.biasdiagnosticsfilestr);
fid = fopen(biasdiagnosticsfilename,'w');
fns = fieldnames(bias_diagnostics);
for i = 1:numel(fns),
  fprintf(fid,'%s',fns{i});
  fprintf(fid,',%f',bias_diagnostics.(fns{i}));
  fprintf(fid,'\n');
end
fclose(fid);
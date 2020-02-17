function [roicenterx,roicentery,roiradii,roiscores] = detectCourtshipBowlROIs(bgmed,params)

nrois = size(params.roimus,1);

isedge = edge(bgmed,'canny',params.cannythresh,params.cannysigma);
% 
% clf;
% imagesc(bgmed);
% axis image;
% colormap gray;
% hold on;
% [r,c] = find(isedge);
% plot(c,r,'r.');

rlim = params.meanroiradius + params.maxdradius*[-1,1];
binedgesr = linspace(rlim(1),rlim(2),params.nbinsradius+1);
bincentersr = (binedgesr(1:end-1)+binedgesr(2:end))/2;
binedgesa = nan(nrois,params.nbinscenter+1);
bincentersb = nan(nrois,params.nbinscenter);
for i = 1:nrois,
  xlim = params.roimus(i,1) + params.maxdcenter*[-1,1];
  ylim = params.roimus(i,2) + params.maxdcenter*[-1,1];
  binedgesa(i,:) = linspace(xlim(1),xlim(2),params.nbinscenter+1);
  binedgesb = linspace(ylim(1),ylim(2),params.nbinscenter+1);
  bincentersb(i,:) = (binedgesb(1:end-1)+binedgesb(2:end))/2;
end

roiradii = nan(1,nrois);
roicenterx = nan(1,nrois);
roicentery = nan(1,nrois);
roiscores = nan(1,nrois);
for i = 1:nrois,
  [roiradii(i),roicenterx(i),roicentery(i),roiscores(i)] = detectcircles(isedge,...
    'binedgesa',binedgesa(i,:),'bincentersb',bincentersb(i,:),...
    'bincentersr',bincentersr,...
    'peaksnhoodsize',[],'peaksthreshold',[],...
    'maxncircles',1,'doedgedetect',false);
end

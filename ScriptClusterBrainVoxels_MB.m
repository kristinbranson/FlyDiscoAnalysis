%% set up paths

addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;

datadir = 'MBClustering';
datafile = fullfile(datadir,'Distribution_of_output_terminals_and_DA_dendrites.txt');

%% parameters

maxdneighbor = 3;
areaopen = 5; % remove connected components with area less than this
allowdistantneighbors = true;
prctcovexpain = 90;
cproplaplace =  0.4553;
kneighbors = 3;

% selected using MBVoxelSymmetrize
rotation = [-0.00420823056252578 -0.00378717698476799 -0.0117908432055422];
% translation = [-1.29790142326049 -0.17598211043431 0.352772397189968];
%rotation = [0,0,0];
%translation = [-1,-1,0];
translation = [-1.00747798915984 0.0310051545304137 -0.00872517032319731];
% this was determined by starting with translation found using
% MBVoxelSymmetrize and then adding the median match translation found
% using hungarian algorithm. 

maxlrmatchcost = 1.5;

%% read data
[data,pos,coltitles,blockids,neuropilnames,neuropilids] = ReadMBVoxelInfo(datafile);
ndata = size(data,1);

blockids0 = blockids;
pos0 = pos;
neuropilnames0 = neuropilnames;
neuropilids0 = neuropilids;
data0 = data;
ndata0 = ndata;

sig = cov(data);
[coeff,eigvals,explained] = pcacov(sig);


% split into left and right 
minpos = min(pos,[],1);
maxpos = max(pos,[],1);
midline = (maxpos+minpos)/2;
posisleft = pos(:,1) < midline(1);
posisright = pos(:,1) > midline(1);

posleft = pos(posisleft,:);
posright = pos(posisright,:);
dataleft = data(posisleft,:);
dataright = data(posisright,:);

posleft0 = posleft;
posright0 = posright;

%% symmetrize left and right sides

posleft1 = posleft0;
posleft1(:,1) = midline(1)-posleft1(:,1);
posright1 = posright0;
posright1(:,1) = posright1(:,1)-midline(1);

Rx = [1,0,0
  0,cos(rotation(1)),-sin(rotation(1))
  0,sin(rotation(1)),cos(rotation(1))];
Ry = [cos(rotation(2)),0,sin(rotation(2))
  0,1,0
  -sin(rotation(2)),0,cos(rotation(2))];
Rz = [cos(rotation(3)),-sin(rotation(3)),0
  sin(rotation(3)),cos(rotation(3)),0
  0,0,1];

posleft1 = bsxfun(@plus,(Rx*Ry*Rz*posleft1')',translation);

% find correspondences
ndataleft = size(posleft1,1);
ndataright = size(posleft1,1);
costmatrix = zeros(ndataleft+ndataright)+maxlrmatchcost;
costmatrix(1:ndataleft,1:ndataright) = dist2(posleft1,posright1);
costmatrix(ndataleft+1:end,ndataright+1:end) = 0;
[assignment,cost] = assignmentoptimal(costmatrix);

l2ridx = assignment(1:ndataleft);
l2ridx(l2ridx > ndataright) = 0;
[~,r2lidx] = ismember((1:ndataright)',assignment);
r2lidx(r2lidx > ndataleft) = 0;

% compare data for visualization
pc1 = coeff(:,1);
projleft1 = dataleft*pc1;
projright1 = dataright*pc1;

figure(1);
clf;
scatter3(posleft1(:,1),posleft1(:,2),posleft1(:,3),[],projleft1,'.');
hold on;
scatter3(posright1(:,1),posright1(:,2),posright1(:,3),[],projright1,'x');
plot3([posleft1(l2ridx>0,1),posright1(l2ridx(l2ridx>0),1)]',...
  [posleft1(l2ridx>0,2),posright1(l2ridx(l2ridx>0),2)]',...
  [posleft1(l2ridx>0,3),posright1(l2ridx(l2ridx>0),3)]','k-');
axis equal;
grid on;
xlabel('x'); ylabel('y'); zlabel('z');

% make one combined data set
poscombined = [(posleft1(l2ridx>0,:)+posright1(l2ridx(l2ridx>0),:))/2
  posleft1(l2ridx==0,:)
  posright1(l2ridx==0,:)];

datacombined = [(dataleft(l2ridx>0,:)+dataright(l2ridx(l2ridx>0),:))/2
  dataleft(l2ridx==0,:)
  dataright(l2ridx==0,:)];
projcombined1 = datacombined*pc1;

figure(2);
clf;
scatter3(posleft0(:,1),posleft0(:,2),posleft0(:,3),[],projleft1,'.')
hold on;
scatter3(posright0(:,1),posright0(:,2),posright0(:,3),[],projright1,'.')
axis equal;
grid on;
xlabel('x'); ylabel('y'); zlabel('z');

figure(3);
clf;
scatter3(posleft1(:,1),posleft1(:,2),posleft1(:,3),[],projleft1,'.')
hold on;
scatter3(posright1(:,1),posright1(:,2),posright1(:,3),[],projright1,'x')
axis equal;
grid on;
xlabel('x'); ylabel('y'); zlabel('z');

figure(4);
clf;
scatter3(poscombined(:,1),poscombined(:,2),poscombined(:,3),[],projcombined1,'.');
axis equal;
grid on;
xlabel('x'); ylabel('y'); zlabel('z');


%% choose neighbors

% [D,idxremove] = MBVoxelSelectNeighbors(pos,maxdneighbor,areaopen,allowdistantneighbors);
% 
% pos = pos0;
% data = data0;
% pos(idxremove,:) = [];
% data(idxremove,:) = [];
% ndata = size(pos,1);

% % left
% [Dleft,idxremove_left] = MBVoxelSelectNeighbors(posleft,maxdneighbor,areaopen,allowdistantneighbors);
% 
% posleft(idxremove_left,:) = [];
% dataleft(idxremove_left,:) = [];
% ndataleft = size(posleft,1);
% 
% % right
% [Dright,idxremove_right] = MBVoxelSelectNeighbors(posright,maxdneighbor,areaopen,allowdistantneighbors);
% 
% posright(idxremove_right,:) = [];
% dataright(idxremove_right,:) = [];
% ndataright = size(posright,1);

[Dcombined,idxremove_combined] = MBVoxelSelectNeighbors(poscombined,maxdneighbor,areaopen,allowdistantneighbors);
poscombined(idxremove_combined,:) = [];
datacombined(idxremove_combined,:) = [];
ndata = size(poscombined,1);

%% visualize neighbors to make sure i didn't screw up

clf;
plot3(poscombined(:,1),poscombined(:,2),poscombined(:,3),'k.');
hold on;
for i = 1:size(poscombined,1),
  idx = ~isinf(Dcombined(i,:));
  z = zeros(nnz(idx),1);
  plot3([poscombined(i,1)+z,poscombined(idx,1)]',[poscombined(i,2)+z,poscombined(idx,2)]',[poscombined(i,3)+z,poscombined(idx,3)]','-');
end

% 
% clf;
% 
% plot3(posleft(:,1),posleft(:,2),posleft(:,3),'k.');
% hold on;
% plot3(posright(:,1),posright(:,2),posright(:,3),'r.');
% for i = 1:ndataleft,
%   idx = ~isinf(Dleft(i,:));
%   z = zeros(nnz(idx),1);
%   plot3([posleft(i,1)+z,posleft(idx,1)]',[posleft(i,2)+z,posleft(idx,2)]',[posleft(i,3)+z,posleft(idx,3)]','-');
% end
% for i = 1:ndataright,
%   idx = ~isinf(Dright(i,:));
%   z = zeros(nnz(idx),1);
%   plot3([posright(i,1)+z,posright(idx,1)]',[posright(i,2)+z,posright(idx,2)]',[posright(i,3)+z,posright(idx,3)]','-');
% end

axis equal
grid on

%% compute affinity matrix

% how many independent dimensions are there, approximately?
% sig = cov(data);
% [coeff,eigvals,explained] = pcacov(sig);
ndimeff = find(cumsum(explained)>=prctcovexpain,1);

% choose sigma
% sigmaleft = MBVoxelChooseSigma(Dleft,dataleft,kneighbors);
% sigmaright = MBVoxelChooseSigma(Dright,dataright,kneighbors);
% sigma = (sigmaleft+sigmaright)/2;
logdatacombined = log(1+datacombined);
[sigma,sigmalocal] = MBVoxelChooseSigma(Dcombined,datacombined,kneighbors);

% [Wleft,wleft] = MBVoxelComputeAffinityMatrix(dataleft,Dleft,ndimeff,cproplaplace,sigma);
% [Wright,wright] = MBVoxelComputeAffinityMatrix(dataright,Dright,ndimeff,cproplaplace,sigma);
[W,w] = MBVoxelComputeAffinityMatrix(datacombined,Dcombined,ndimeff,cproplaplace,sigma);

%% visualize weights

clf;

[r,c] = find(~isinf(Dcombined));
projcombined1 = datacombined*pc1;
scatter3(poscombined(:,1),poscombined(:,2),poscombined(:,3),[],sigmalocal,'.');
hold on;
meanw = mean(w);
[~,order] = sort(w);
for i = order(:)',
  i1 = r(i);
  i2 = c(i);
  wcurr = w(i);
  if wcurr <= meanw,
    colorcurr = .5*wcurr/meanw;
  else
    colorcurr = (wcurr-meanw)*.5 + .5;
  end
  plot3(poscombined([i1,i2],1),poscombined([i1,i2],2),poscombined([i1,i2],3),'-','Color',colorcurr+[0,0,0]);
end

% [r,c] = find(~isinf(Dright));
% plot3(posright(:,1),posright(:,2),posright(:,3),'.','Color',[.7,0,0]);
% hold on;
% [~,order] = sort(wright);
% for i = order(:)',
%   i1 = r(i);
%   i2 = c(i);
%   wcurr = wright(i);
%   plot3(posright([i1,i2],1),posright([i1,i2],2),posright([i1,i2],3),'-','Color',[1,1-wcurr,1-wcurr]);
% end

axis equal
grid on

%% spectral clustering

nclusterstry = 2:10;
% labelsleft = cell(1,numel(nclusterstry));
% for nclustersi = 1:numel(nclusterstry),
%   nclusters = nclusterstry(nclustersi);
%   [NcutDiscrete,NcutEigenvectors,NcutEigenvalues] = ncutW(Wleft,nclusters);
%   [labelsleft{nclustersi},~] = find(NcutDiscrete');
% end
% 
% labelsright = cell(1,numel(nclusterstry));
% for nclustersi = 1:numel(nclusterstry),
%   nclusters = nclusterstry(nclustersi);
%   [NcutDiscrete,NcutEigenvectors,NcutEigenvalues] = ncutW(Wright,nclusters);
%   [labelsright{nclustersi},~] = find(NcutDiscrete');
% end

labelscombined = cell(1,numel(nclusterstry));
for nclustersi = 1:numel(nclusterstry),
  nclusters = nclusterstry(nclustersi);
  [NcutDiscrete,NcutEigenvectors,NcutEigenvalues] = ncutW(W,nclusters);
  [labelscombined{nclustersi},~] = find(NcutDiscrete');
end

% %% plot clustering results
% 
% clf;
% 
% nc = ceil(sqrt(numel(nclusterstry)));
% nr = ceil(numel(nclusterstry)/nc);
% 
% hax = createsubplots(nr,nc,.05);
% 
% for nclustersi = 1:numel(nclusterstry),
%   nclusters = nclusterstry(nclustersi);
%   colors = jet(2*nclusters)*.7;
%   colors = colors(randperm(2*nclusters),:);
%   colorsleft = colors(1:nclusters,:);
%   colorsright = colors(nclusters+1:end,:);
%   markers = '.';
%   markers = markers(randsample(numel(markers),2*nclusters,1));
%   markersleft = markers(1:nclusters);
%   markersright = markers(nclusters+1:end);
%   
%   for i = 1:nclusters,
%     axes(hax(nclustersi));
%     idx = labelsleft{nclustersi} == i;
%     plot3(posleft(idx,1),posleft(idx,2),posleft(idx,3),markersleft(i),'Color',colorsleft(i,:));
%     if i == 1,
%       hold on;
%     end
%     idx = labelsright{nclustersi} == i;
%     plot3(posright(idx,1),posright(idx,2),posright(idx,3),markersright(i),'Color',colorsright(i,:));
%     
%   end
%   axis equal;
%   grid on;
%   xlabel('x'); ylabel('y'); zlabel('z');
%   title(sprintf('k = %d',nclusters*2));
% end
% 
% linkprop(hax(1:numel(nclusterstry)),{'XLim','YLim','ZLim',...
%   'CameraPosition','CameraTarget','CameraUpVector','CameraViewAngle'});

%% plot clustering results

clf;

nc = ceil(sqrt(numel(nclusterstry)));
nr = ceil(numel(nclusterstry)/nc);

hax = createsubplots(nr,nc,.05);

for nclustersi = 1:numel(nclusterstry),
  nclusters = nclusterstry(nclustersi);
  colors = jet(nclusters)*.7;
  colors = colors(randperm(nclusters),:);
  markers = '.';
  markers = markers(randsample(numel(markers),nclusters,1));
  
  for i = 1:nclusters,
    axes(hax(nclustersi));
    idx = labelscombined{nclustersi} == i;
    plot3(poscombined(idx,1),poscombined(idx,2),poscombined(idx,3),markers(i),'Color',colors(i,:));
    if i == 1,
      hold on;
    end
  end
  axis equal;
  grid on;
  xlabel('x'); ylabel('y'); zlabel('z');
  title(sprintf('k = %d',nclusters));
end

linkprop(hax(1:numel(nclusterstry)),{'XLim','YLim','ZLim',...
  'CameraPosition','CameraTarget','CameraUpVector','CameraViewAngle'});

%% output clusterings to tiff files

% timestamp = datestr(now,'yyyymmdd');
% 
% for nclustersi = 1:numel(nclusterstry),
%   
%   im = zeros(maxpos,'uint8');
%   nclusters = nclusterstry(nclustersi);
%   for i = 1:nclusters,
%     
%     idx1 = labelsleft{nclustersi} == i;
%     idx2 = sub2ind(maxpos,posleft(idx1,1),posleft(idx1,2),posleft(idx1,3));
%     im(idx2) = 2*i-1;
%     idx1 = labelsright{nclustersi} == i;
%     idx2 = sub2ind(maxpos,posright(idx1,1),posright(idx1,2),posright(idx1,3));
%     im(idx2) = 2*i;
%     
%   end
% 
%   im = im(minpos(1)+1:end,minpos(2)+1:end,minpos(3)+1:end);
%   im = permute(im,[2,1,3]);
%   outfilename = sprintf('ClusterLabels_k%d_%s.tiff',2*nclusters,timestamp);
%   imwrite(im(:,:,1),outfilename,'tiff');
%   for i = 2:size(im,3),
%     imwrite(im(:,:,i),outfilename,'WriteMode','append');
%   end
%   
% end
  
%% output clusterings to tiff files

timestamp = datestr(now,'yyyymmdd');

maxposcurr = max(round(poscombined),[],1);
minposcurr = min(round(poscombined),[],1);
for nclustersi = 1:numel(nclusterstry),
  
  im = zeros(maxposcurr,'uint8');
  nclusters = nclusterstry(nclustersi);
  for i = 1:nclusters,
    
    idx1 = labelscombined{nclustersi} == i;
    idx2 = sub2ind(maxposcurr,round(poscombined(idx1,1)),round(poscombined(idx1,2)),round(poscombined(idx1,3)));
    im(idx2) = i;
    
  end

  im = im(minposcurr(1)+1:end,minposcurr(2)+1:end,minposcurr(3)+1:end);
  im = permute(im,[2,1,3]);
  outfilename = fullfile('MBClustering',sprintf('ClusterLabels_k%d_%s.tiff',nclusters,timestamp));
  imwrite(im(:,:,1),outfilename,'tiff');
  for i = 2:size(im,3),
    imwrite(im(:,:,i),outfilename,'WriteMode','append');
  end
  
end
  
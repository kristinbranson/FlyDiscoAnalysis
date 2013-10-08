expdir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data/GMR_10B06_AE_01_TrpA_Rig2Plate14BowlB_20110728T094915';
hfig = 8765;
windowradius = 15;
linewidth = 2;

moviefile = fullfile(expdir,'movie.ufmf');
trxfile = fullfile(expdir,'registered_trx.mat');
[readframe,nframes,fid,headerinfo] = get_readframe_fcn(moviefile);
load(trxfile);
sd = load(fullfile(expdir,'scores_pivot_center.mat'));
nflies = numel(trx);

% compute the maximum distance between all pairs of flies for each frame
endframe = max([trx.endframe]);
d = nan(endframe,nflies,nflies);
for fly1 = 1:nflies,
  for fly2 = fly1+1:nflies,
    
    t0 = max(trx(fly1).firstframe,trx(fly2).firstframe);
    t1 = min(trx(fly1).endframe,trx(fly2).endframe);
    if t0 > t1,
      continue;
    end
    
    d(t0:t1,fly1,fly2) = (trx(fly1).x_mm(t0+trx(fly1).off:t1+trx(fly1).off)-...
      trx(fly2).x_mm(t0+trx(fly2).off:t1+trx(fly2).off)).^2 + ...
      (trx(fly1).y_mm(t0+trx(fly1).off:t1+trx(fly1).off)-...
      trx(fly2).y_mm(t0+trx(fly2).off:t1+trx(fly2).off)).^2;
    %maxd(t0:t1) = max(maxd(t0:t1),d);
  end
end

idxremove = all(isnan(d),1);
d = reshape(d,[endframe,nflies*nflies]);
d = d(:,~idxremove(:));

maxd = prctile(d,75,2);
maxd = sqrt(maxd);

npivotswin = nan(nflies,endframe);
for fly = 1:nflies,
  
  npivots = sd.allScores.postprocessed{fly};
  npivots(isnan(npivots)) = 0;
  npivotswin(fly,1:numel(npivots)) = conv(npivots,ones(1,windowradius*2+1),'same');
  
end  

npivotsallflies = nansum(npivotswin,1);

figure(hfig);
clf;
hax = [nan,nan];
hax(1) = subplot(2,1,1);
plot(maxd,'k.-');
axisalmosttight;
title('cluster radius');
hax(2) = subplot(2,1,2);
plot(npivotsallflies,'r.-');
axisalmosttight;
title('npivots');
linkaxes(hax,'x');
input('Zoom in on relevant region, hit enter, then click to select frame.');
[T,~] = ginput(1);
T = round(T);

im = readframe(T);
%%
figure(hfig);
clf;
hax = gca;
% otherwise the colorbar is wrong
imagesc([-windowradius/30,windowradius/30]);
hold on;
him = image(repmat(im,[1,1,3]));
axis image;

T0 = T-windowradius;
T1 = T+windowradius;
colors = jet(2*windowradius+1)*.75;
htri = nan(1,nflies);
hline = nan(2*windowradius+1,nflies);
htrx = nan(1,nflies);

minx = inf;
maxx = 0;
miny = inf;
maxy = 0;

for fly = 1:nflies,
  if T < trx(fly).firstframe || T > trx(fly).endframe,
    continue;
  end
  
  i = T+trx(fly).off;
  htri(fly) = drawflyo(trx(fly).x(i),trx(fly).y(i),trx(fly).theta(i),trx(fly).a(i),trx(fly).b(i),'-','Color','k','LineWidth',linewidth);
  
  i0 = T0+trx(fly).off;
  i1 = T1+trx(fly).off;
  for ii = 1:2*windowradius+1,
    i = i0 + ii - 1;
    if i < 1 || i > trx(fly).nframes,
      continue;
    end
    x0 = trx(fly).x(i)+2*trx(fly).a(i)*cos(trx(fly).theta(i));
    x1 = trx(fly).x(i)-2*trx(fly).a(i)*cos(trx(fly).theta(i));
    y0 = trx(fly).y(i)+2*trx(fly).a(i)*sin(trx(fly).theta(i));
    y1 = trx(fly).y(i)-2*trx(fly).a(i)*sin(trx(fly).theta(i));
    hline(ii,fly) = plot([x0,x1],[y0,y1],...
      '-','LineWidth',linewidth,'Color',colors(ii,:));
    minx = min([minx,x0,x1]);
    maxx = max([maxx,x0,x1]);
    miny = min([miny,y0,y1]);
    maxy = max([maxy,y0,y1]);
  end

  i0 = max(i0,1);
  i1 = min(i1,trx(fly).nframes);
  
  htrx(fly) = plot(trx(fly).x(i0:i1),trx(fly).y(i0:i1),'k.-');
  
end

dx = maxx - minx;
dy = maxy - miny;
ax = [minx-dx*.1,maxx+dx*.1,miny-dy*.1,maxy+dy*.1];
axis(ax);

set(gca,'CLim',[-windowradius/30,windowradius/30]);
colormap(colors);
axis off;
colorbar;

%% save result

[~,name] = fileparts(expdir);
basename = sprintf('clusteringflies/ExampleCluster_%s_%05d',name,T);
if ~exist('clusteringfiles','dir'),
  mkdir clusteringflies;
end
SaveFigLotsOfWays(hfig,basename);

%%
fclose(fid);
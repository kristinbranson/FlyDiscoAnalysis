function [fore2fly,xgrid_isfore,ygrid_isfore,debugdata,istouching] = TrackWings_SegmentFlies(...
  im,isfore_thresh,idxfore_thresh,npxfore_thresh,fore2body,...
  trx,params,XGRID,YGRID,debugdata)

persistent hx hy se_boundary;

%% assign targets to connected components
nflies = numel(trx);

% whether each fly is in the same cc as another
istouching = false(1,nflies);
[L,ncc] = bwlabel(isfore_thresh);
fore2cc = L(isfore_thresh);
fly2cc = nan(1,nflies);
pcurr = [];
[nr,nc] = size(XGRID);
for fly = 1:nflies,
  
  if isempty(trx(fly).x),
    continue;
  end
  x = round(trx(fly).x);
  y = round(trx(fly).y);
  if y >= 1 && y <= nr && x >= 1 && x <= nc && ...
      isfore_thresh(y,x),
    fly2cc(fly) = L(y,x);
  else
    if isempty(pcurr),
      pcurr = [XGRID(isfore_thresh),YGRID(isfore_thresh)];
    end
    mu = [trx(fly).x,trx(fly).y];
    S = axes2cov(trx(fly).a*2,trx(fly).b*2,trx(fly).theta);
    diffs = bsxfun(@minus,pcurr,mu);
    c = chol(S);
    temp = diffs/c;
    dcurr = sum(temp.^2, 2);
    [~,j] = min(dcurr);
    fly2cc(fly) = fore2cc(j);
  end
  
end

%% assign unassigned connected components to flies

fly2nccs = ones(1,nflies);
unassignedcc = setdiff(1:ncc,fly2cc);
if ~isempty(unassignedcc),
  mus = nan(nflies,2);
  for fly = 1:nflies,
    if isempty(trx(fly).x),
      continue;
    end
    mus(fly,:) = [trx(fly).x,trx(fly).y];
  end
  for cci = unassignedcc,
    % remove small connected components
    if nnz(L==cci) < params.min_wingcc_area,
      continue;
    end
    x = XGRID(L==cci);
    y = YGRID(L==cci);
    D = dist2([x(:),y(:)],mus);
    d = min(D,[],1);
    [mind,fly] = min(d);
    mind = mind/(trx(fly).a*4);
    if mind > params.max_wingcc_dist,
      continue;
    end
    L(L==cci) = fly2cc(fly);
    fore2cc(fore2cc==cci) = fly2cc(fly);
    fly2nccs(fly) = fly2nccs(fly) + 1;
  end
end

if debugdata.DEBUG > 1,
  
  DebugPlot_SegmentFlies1();
  
end



%% loop over connected components, assign pixels to flies

unique_cc = unique(fly2cc(~isnan(fly2cc)));
% fore2fly corresponds to isfore_thresh
fore2fly = nan(npxfore_thresh,1);

xgrid_isfore = XGRID(isfore_thresh);
ygrid_isfore = YGRID(isfore_thresh);

% for watershed
if isempty(hy),
  hy = fspecial('sobel');
  hx = hy';
  se_boundary = strel('disk',1);
end

for cci = unique_cc,
  flies = find(fly2cc == cci);
  fore2curr = fore2cc == cci;
  % all pixels belong to one fly
  if numel(flies) == 1,
    fore2fly(fore2curr) = flies;
    continue;
  end
  
  istouching(flies) = true;
  
  % assign body pixels based on mahdist
  
  % indices of foreground pixels in the current cc
  curr2fore = find(fore2curr);
  
  % ncurr x 1, whether current cc pixels are body or not
  curr2body = fore2body(fore2curr);
  
  % nbodycurr x 1, indices of current cc pixels in body
  body2curr = find(curr2body);
  nbodycurr = numel(body2curr);
  
  % compute distance from each current body pixel to fly
  mind = inf(nbodycurr,1);
  body2id = nan(nbodycurr,1);
  body2isback = false(nbodycurr,1);
  pbody = [xgrid_isfore(curr2fore(body2curr)),ygrid_isfore(curr2fore(body2curr))];
  for fly = flies,
    mu = [trx(fly).x,trx(fly).y];
    S = axes2cov(trx(fly).a*2,trx(fly).b*2,trx(fly).theta);
    diffs = bsxfun(@minus,pbody,mu);
    c = chol(S);
    temp = diffs/c;
    dcurr = sum(temp.^2, 2);
    isbackcurr = diffs(:,1)*cos(trx(fly).theta) + diffs(:,2)*sin(trx(fly).theta) <= 0;
    %isbackcurr = abs(modrange(atan2(diffs(:,2),diffs(:,1))-trx(fly).theta,-pi,pi)) >= pi/2;
    
    idx = dcurr < mind;
    body2id(idx) = fly;
    mind(idx) = dcurr(idx);
    body2isback(idx) = isbackcurr(idx);
  end
  
  % compute distance from each non-body pixel to each body pixel
  pcurr = [xgrid_isfore(fore2curr),ygrid_isfore(fore2curr)];
  xlims = [floor(min(pcurr(:,1)))-1,ceil(max(pcurr(:,1)))+1];
  ylims = [floor(min(pcurr(:,2)))-1,ceil(max(pcurr(:,2)))+1];
  xlims = max(1,min(nc,xlims));
  ylims = max(1,min(nr,ylims));
  dx = diff(xlims)+1;
  dy = diff(ylims)+1;
  
  % pixels we can go through
  bw = L(ylims(1):ylims(2),xlims(1):xlims(2)) == cci;
  mind = inf(size(bw));
  flycurr = nan(size(bw));
  for fly = flies,
    % pixels we are finding the distance to
    c = pbody(body2id==fly&body2isback,1)-xlims(1)+1;
    r = pbody(body2id==fly&body2isback,2)-ylims(1)+1;
    
    % find distance within foreground pixels
    dcurr = bwdistgeodesic(bw,c,r,'quasi-euclidean');
    
    % assign if smaller
    tmpidx = dcurr < mind;
    mind(tmpidx) = dcurr(tmpidx);
    flycurr(tmpidx) = fly;
  end
  
  % pixels that have infinite distance
  infidx = bw & isinf(mind);
  if any(infidx(:)),
    [infy,infx] = find(infidx);
    infmind = inf(1,numel(infy));
    infflycurr = nan(1,numel(infy));
    for fly = flies,
      % pixels we are finding the distance to
      c = pbody(body2id==fly&body2isback,1)-xlims(1)+1;
      r = pbody(body2id==fly&body2isback,2)-ylims(1)+1;
      if isempty(c),
        continue;
      end

      % find distance to a foreground pixel through any pixel
      dcurr = min(dist2([c(:),r(:)],[infx(:),infy(:)]),[],1);
    
      % assign if smaller
      tmpidx = dcurr < infmind;
      infmind(tmpidx) = dcurr(tmpidx);
      infflycurr(tmpidx) = fly;
    end
    mind(infidx) = infmind;
    flycurr(infidx) = infflycurr;
  end
  
  % superpixel segmentation
  imbb = im(ylims(1):ylims(2),xlims(1):xlims(2));
  imbb(~isfore_thresh(ylims(1):ylims(2),xlims(1):xlims(2))) = ...
    mean(imbb(~isfore_thresh(ylims(1):ylims(2),xlims(1):xlims(2))));
  Iy = imfilter(imbb, hy, 'replicate');
  Ix = imfilter(imbb, hx, 'replicate');
  gradmag = sqrt(Ix.^2 + Iy.^2);
  isboundary = ~bw & imdilate(bw,se_boundary);
  gradmag(isboundary) = inf;
  
  Lseg = watershed(gradmag);
  
  % remove bg segments
  nsegcurr = max(Lseg(:));
  for segi = 1:nsegcurr,
    segcurr = Lseg==segi;
    tmp = bw(segcurr);
    tmp(isnan(tmp)) = [];
    if isempty(tmp),
      continue;
    end
    if nnz(tmp) < numel(tmp)/2,
      Lseg(segcurr) = nan;
    end
  end
  
  imbb(~bw|isnan(Lseg)) = nan;
  
  % assign watershed pixels based on intensity difference
  [rw,cw] = find(Lseg==0);
  iw = sub2ind(size(Lseg),rw,cw);
  nwatershed = numel(rw);
  mind = inf(nwatershed,1);
  Lw = nan(nwatershed,1);
  for dc = -1:1,
    for dr = -1:1,
      if dc == 0 && dr == 0,
        continue;
      end
      % adjacent pixel
      r1 = rw + dr;
      c1 = cw + dc;
      % make sure it is in bounds
      idxw = find(r1 >= 1 & r1 <= dy & c1 >= 1 & c1 <= dx);
      iw1 = sub2ind(size(Lseg),r1(idxw),c1(idxw));
      % also make sure it does not have Lseg == 0, and is foreground
      idxw(Lseg(iw1)==0 | ~bw(iw1)) = [];
      iw1 = sub2ind(size(Lseg),r1(idxw),c1(idxw));
      % compute intensity difference
      dcurr = nan(nwatershed,1);
      lcurr = nan(nwatershed,1);
      dcurr(idxw) = abs(imbb(iw(idxw))-imbb(iw1));
      lcurr(idxw) = Lseg(iw1);
      idxw1 = dcurr < mind;
      if any(isnan(lcurr(idxw1))),
        keyboard;
      end
      %fprintf('dr = %d, dc = %d, assigning for: %s\n',dr,dc,mat2str(find(idxw1)));
      mind(idxw1) = dcurr(idxw1);
      Lw(idxw1) = lcurr(idxw1);
    end
  end
  
  Lseg(iw) = Lw;
  Lseg(~bw) = nan;
  
  % for each superpixel, take the majority vote
  nsegcurr = max(Lseg(:));
  flycurr1 = nan(size(flycurr));
  for segi = 1:nsegcurr,
    segcurr = Lseg==segi;
    tmp = flycurr(segcurr);
    tmp(isnan(tmp)) = [];
    if isempty(tmp),
      continue;
    end
    flyseg = mode(tmp);
    flycurr1(segcurr) = flyseg;
  end
  
  flycurr1(~bw|isnan(flycurr1)) = 0;
  
  % store
  fore2fly(fore2curr) = flycurr1(bw);
  
end

if debugdata.DEBUG > 1,
  DebugPlot_SegmentFlies2();
end

%% debug plot stuff


  function DebugPlot_SegmentFlies1()
    
    axcurr = 2;
    if isnan(debugdata.hims(axcurr)),
      cla(debugdata.hax(axcurr));
    end
    
    imtmp = repmat(debugdata.im(:),[1,3]);
    for dfly = 1:nflies,
      idx1 = L == fly2cc(dfly);
      imtmp(idx1,:) = min(bsxfun(@plus,imtmp(idx1,:)*3,255*debugdata.colors(dfly,:))/4,255);
    end
    imtmp = uint8(reshape(imtmp,[nr,nc,3]));
    
    if isnan(debugdata.hims(axcurr)),
      debugdata.hims(axcurr) = image(imtmp,'Parent',debugdata.hax(axcurr)); axis(debugdata.hax(axcurr),'image','off');
    else
      set(debugdata.hims(axcurr),'CData',imtmp);
      delete(debugdata.htext{axcurr}(ishandle(debugdata.htext{axcurr})));
    end
    debugdata.htext{axcurr} = [];
    
    for dcci = unique(fly2cc),
      flies = find(fly2cc == dcci);
      s = sprintf('%d,',flies);
      s = s(1:end-1);
      x = mean(XGRID(L==dcci));
      y = mean(YGRID(L==dcci));
      debugdata.htext{axcurr} = [debugdata.htext{axcurr},text(x,y,s,'HorizontalAlignment','center','VerticalAlignment','middle','Clipping','on','Parent',debugdata.hax(axcurr))];
    end
    
    title(debugdata.hax(axcurr),sprintf('Connected components for frame %d',debugdata.t));
    drawnow;
    
  end

  function DebugPlot_SegmentFlies2()
    
    
    axcurr = 3;
    if isnan(debugdata.hims(axcurr)),
      cla(debugdata.hax(axcurr));
    end
    
    imtmp = repmat(debugdata.im(:),[1,3]);
    for dfly = 1:nflies,
      idx1 = idxfore_thresh(fore2fly==dfly);
      imtmp(idx1,:) = min(bsxfun(@plus,imtmp(idx1,:)*3,255*debugdata.colors(dfly,:))/4,255);
    end
    
    imtmp = uint8(reshape(imtmp,[nr,nc,3]));
    if isnan(debugdata.hims(axcurr)),
      debugdata.hims(axcurr) = image(imtmp,'Parent',debugdata.hax(axcurr)); axis(debugdata.hax(axcurr),'image','off');
    else
      set(debugdata.hims(axcurr),'CData',imtmp);
      delete(debugdata.htext{axcurr}(ishandle(debugdata.htext{axcurr})));
    end
    debugdata.htext{axcurr} = [];
    
    for dfly = 1:nflies,
      s = sprintf('%d',dfly);
      idx1 = idxfore_thresh(fore2fly==dfly);
      x = mean(XGRID(idx1));
      y = mean(YGRID(idx1));
      debugdata.htext{axcurr} = [debugdata.htext{axcurr},text(x,y,s,'HorizontalAlignment','center','VerticalAlignment','middle','Clipping','on','Parent',debugdata.hax(axcurr))];
    end
    
    title(debugdata.hax(axcurr),sprintf('Segmentation for frame %d',debugdata.t));
    drawnow;
    
  end

end

%% plot segmentation of one connected component into flies
%
% tmpplot = repmat(imbb(:),[1,3]);
% debugdata.colors = lines(nflies);
% debugdata.colors = debugdata.colors(randperm(size(debugdata.colors,1)),:);
% for i = 1:nflies,
%   idx1 = flycurr1==flies(i);
%   tmpplot(idx1,:) = min(bsxfun(@plus,tmpplot(idx1,:),255*debugdata.colors(i,:))/2,255);
% end
%
% image(uint8(reshape(tmpplot,[size(imbb),3]))); axis(debugdata.hax(axcurr),'image','off')

%% plot segmentation of one connected component into superpixels
%
% tmpplot = repmat(imbb(:),[1,3]);
% nseg = double(max(Lseg(:)))+1;
% debugdata.colors = hsv(nseg);
% debugdata.colors = debugdata.colors(randperm(size(debugdata.colors,1)),:);
% for i = 1:nseg+1,
%   idx1 = Lseg==i-1;
%   if ~any(idx1(:)),
%     continue;
%   end
%   tmpplot(idx1,:) = min(bsxfun(@plus,tmpplot(idx1,:)*3,255*debugdata.colors(i,:))/4,255);
% end
%
% clf;
% image(uint8(reshape(tmpplot,[size(imbb),3]))); axis(debugdata.hax(axcurr),'image','off')
% hold on;
% for i = 1:nseg+1,
%   [rtmp,ctmp] = find(Lseg==i-1);
%   if ~isempty(rtmp),
%     text(mean(ctmp),mean(rtmp),num2str(i-1),'horizontalalignment','center');
%   end
% end

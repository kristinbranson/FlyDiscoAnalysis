function [wingtrx,debugdata] = TrackWings_FitWings(fore2fly,xgrid_isfore,ygrid_isfore,isfore_thresh,iswing,fore2dbkgd,idxfore_thresh,trx,wingtrxprev,params,debugdata)

% initialize output

nflies = numel(trx);
wingtrx = struct('wing_anglel',cell(1,nflies),...
  'wing_angler',cell(1,nflies),...
  'nwingsdetected',cell(1,nflies),...
  'wing_areal',cell(1,nflies),...
  'wing_arear',cell(1,nflies),...
  'wing_trough_angle',cell(1,nflies));

fore2wing = iswing(isfore_thresh);
fore2flywing = zeros(size(fore2fly));
for fly = 1:nflies,
  
  % find pixels that belong to each fly's wings
  if isempty(trx(fly).x),
    continue;
  end
  x = trx(fly).x;
  y = trx(fly).y;
  theta = trx(fly).theta;
  idxcurr = fore2wing&fore2fly==fly;
  xwing = xgrid_isfore(idxcurr);
  ywing = ygrid_isfore(idxcurr);
  dthetawing = modrange(atan2(ywing-y,xwing-x)-(theta+pi),-pi,pi);
  isallowed = abs(dthetawing) <= params.max_wingpx_angle;
  idxcurr = find(idxcurr);
  idxcurr(~isallowed) = [];
  fore2flywing(idxcurr) = fly;
  dthetawing(~isallowed) = [];
  
  % fit wings
  
  switch params.wing_fit_method,
    
    case 'peaks',
      
      wingtrx(fly) = TrackWings_FitWings_Peak(dthetawing,params);
      
    case 'gmm',
      
      wingtrx(fly) = TrackWings_FitWings_GMM(dthetawing,fore2dbkgd,wingtrxprev,params);
      
  end
  
end

% plot

if debugdata.DEBUG,
  
  DebugPlot_FitWings();
  
end

  function DebugPlot_FitWings()
    
    debugdata.hfig = 2;
    axcurr = 4;

    imtmp = repmat(debugdata.im(:),[1,3]);
    for dfly = 1:nflies,
      idx1 = idxfore_thresh(fore2flywing==dfly);
      imtmp(idx1,:) = min(bsxfun(@plus,imtmp(idx1,:)*3,255*debugdata.colors(dfly,:))/4,255);
    end
    [nr,nc,~] = size(debugdata.im);
    
    if isnan(debugdata.hims(axcurr)) || ~ishandle(debugdata.hims(axcurr)),
      fprintf('clearing axis\n');
      figure(debugdata.hfig);
      clf;
      debugdata.hax1 = gca;
      cla(debugdata.hax1);
      debugdata.hims(axcurr) = image(uint8(reshape(imtmp,[nr,nc,3])),'Parent',debugdata.hax1);
      axis(debugdata.hax1,'image','off');
      hold(debugdata.hax1,'on');
      if debugdata.DEBUG > 1,
        linkaxes([debugdata.hax,debugdata.hax1]);
      end
    else
      set(debugdata.hims(axcurr),'CData',uint8(reshape(imtmp,[nr,nc,3])));
    end
    
    if isfield(debugdata,'hwing'),
      delete(debugdata.hwing(ishandle(debugdata.hwing)));
    end
    %if isfield(debugdata,'debugdata.htext2'),
    %  delete(debugdata.htext2(ishandle(debugdata.htext2)));
    %end
    if isfield(debugdata,'htrough'),
      delete(debugdata.htrough(ishandle(debugdata.htrough)));
    end
    debugdata.hwing = [];
    debugdata.htrough = [];
    %debugdata.htext2 = [];
    for dfly = 1:nflies,
      if isempty(trx(dfly).x),
        continue;
      end
      xwing = [nan,trx(dfly).x,nan];
      ywing = [nan,trx(dfly).y,nan];
      wing_angles = [wingtrx(dfly).wing_anglel,wingtrx(dfly).wing_angler];
      xwing([1,3]) = trx(dfly).x + 4*trx(dfly).a*cos(trx(dfly).theta+pi+wing_angles);
      ywing([1,3]) = trx(dfly).y + 4*trx(dfly).a*sin(trx(dfly).theta+pi+wing_angles);
      xtrough = trx(dfly).x+2*trx(dfly).a*cos(trx(dfly).theta+pi+wingtrx(dfly).wing_trough_angle);
      ytrough = trx(dfly).y+2*trx(dfly).a*sin(trx(dfly).theta+pi+wingtrx(dfly).wing_trough_angle);
      debugdata.hwing(end+1) = plot(debugdata.hax1,xwing,ywing,'.-','color',debugdata.colors(dfly,:));
      debugdata.htrough(end+1) = plot(debugdata.hax1,xtrough,ytrough,'x','color',debugdata.colors(dfly,:));
      %debugdata.htext2(end+1) = text(xwing(1),ywing(1),sprintf('%.1f',wingtrx(dfly).wing_areal));
      %debugdata.htext2(end+1) = text(xwing(3),ywing(3),sprintf('%.1f',wingtrx(dfly).wing_arear));
    end
    
    if isfield(debugdata,'htext') && numel(debugdata.htext) >= axcurr,
      delete(debugdata.htext{axcurr}(ishandle(debugdata.htext{axcurr})));
    end
    
    debugdata.htext{axcurr} = [];
    
    for dfly = 1:nflies,
      if isempty(trx(dfly).x),
        continue;
      end
      s = sprintf('%d',dfly);
      x = trx(dfly).x;
      y = trx(dfly).y;
      debugdata.htext{axcurr} = [debugdata.htext{axcurr},text(x,y,s,'HorizontalAlignment','center','VerticalAlignment','middle','Clipping','on','Parent',debugdata.hax1,'Color','w')];
    end
    
    if debugdata.DEBUG > 1,
      input('');
      %drawnow;
      %     tmpoutfilename = sprintf('WingTrackingExamples_Segmentation_%05d.pdf',t);
      %     savefig(tmpoutfilename,1,'pdf');
      %     tmpoutfilename = sprintf('WingTrackingExamples_Results_%05d.pdf',t);
      %     savefig(tmpoutfilename,2,'pdf');
    else
      drawnow;
    end
  end
end

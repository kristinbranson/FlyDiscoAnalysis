function saveRegistrationImage(imsavename, hfig, figpos, bkgdImage, circleCenterX, circleCenterY, circleRadius, ...
                               nBowlMarkers, bowlMarkerPoints, circleRadius_mm, originX, originY, offTheta, A, ...
                               bowlMarkerPairTheta_true, scale) 

if isempty(imsavename) ,
  return
end

if isempty(hfig) ,
  hfig = figure() ;
else
  figure(hfig);
end
clf(hfig);
set(hfig,'Units','pixels','Position',figpos);

set(hfig,'Units','pixels','Position',figpos);
nimsplot = 2;
hax = createsubplots(1,nimsplot,0.025,hfig);

% plot background image in jet colormap
imagesc(bkgdImage,'parent',hax(1),[0,255]);
axis(hax(1),'image');
axis(hax(1),'xy');  % This means image will be upside-sown from how Matlab normally displays images
hold(hax(1),'on');

% plot detected/labeled circles
chamber_count = numel(circleRadius) ;
theta_line = linspace(0,2*pi,ceil(2*pi*max(circleRadius)));
for k = 1 : chamber_count ,
  r = circleRadius(k) ;
  xc  = circleCenterX(k) ;
  yc  = circleCenterY(k) ;  
  plot(...
    hax(1), ...
    xc + r*cos(theta_line),...
    yc + r*sin(theta_line),...
    'w-','linewidth',2);
  plot(...
    hax(1), ...
    xc + r*cos(theta_line),...
    yc + r*sin(theta_line),...
    'k--','linewidth',2);
  plot(hax(1),xc,yc,'ks','markerfacecolor','k');
  plot(hax(1),xc,yc,'wd');
end

% plot bowl marker(s)
if nBowlMarkers > 0,
  plot(hax(1),bowlMarkerPoints(1,:),bowlMarkerPoints(2,:),'ks','markersize',12);
  plot(hax(1),bowlMarkerPoints(1,:),bowlMarkerPoints(2,:),'kx','markersize',12);
end

% plot axes
xangle = 0;
yangle = pi/2;
l = circleRadius_mm/2 / scale;
xc_mean = mean(originX) ;
yc_mean = mean(originY) ;
% plot x axis and label
vxx = cos(xangle-offTheta)*l ;
vxy = sin(xangle-offTheta)*l ;
quiver(hax(1),xc_mean,yc_mean,vxx,vxy,0,'k-');
quiver(hax(1),xc_mean,yc_mean,vxx,vxy,0,'w--');
text(xc_mean+vxx,yc_mean+vxy,'x','parent',hax(1));
% plot y axis and label
vyx = cos(yangle-offTheta)*l ;
vyy = sin(yangle-offTheta)*l ;
quiver(hax(1),xc_mean,yc_mean,vyx,vyy,0,'w-');
quiver(hax(1),xc_mean,yc_mean,vyx,vyy,0,'k--');
text(xc_mean+vyx,yc_mean+vyy,'y','parent',hax(1));
% Make title
title(hax(1),'Input background image');

% transform image
A_mean = affineTransformMatrixFromOffsetsAndScale(-xc_mean,-yc_mean,offTheta,scale) ;
T = maketform('affine',A_mean);  %#ok<MTFA1>
[imreg,xdata,ydata] = imtransform(bkgdImage,T,'XYScale',scale,'FillValues',nan);  %#ok<DIMTRNS>

% draw registered image
imagesc(xdata,ydata,imreg,'parent',hax(2),[0,255]);
hold(hax(2),'on');

% draw angle to bowl marker
l = circleRadius_mm;
plot(hax(2),[0,cos(bowlMarkerPairTheta_true)]*l,...
  [0,sin(bowlMarkerPairTheta_true)]*l,'k-o','markerfacecolor','k');
plot(hax(2),[0,cos(bowlMarkerPairTheta_true)]*l,...
  [0,sin(bowlMarkerPairTheta_true)]*l,'w:');
text(cos(bowlMarkerPairTheta_true)*l/2,...
  sin(bowlMarkerPairTheta_true)*l/2,sprintf('%.2fmm',l));
text(cos(bowlMarkerPairTheta_true)*l/5,...
  sin(bowlMarkerPairTheta_true)*l/5,sprintf('%.1fdeg',bowlMarkerPairTheta_true*180/pi));

% draw axes
l = circleRadius_mm/2;
quiver(hax(2),0,0,cos(xangle)*l,sin(xangle)*l,0,'k-');
quiver(hax(2),0,0,cos(xangle)*l,sin(xangle)*l,0,'w--');
text(cos(xangle)*l,sin(xangle)*l,'x','parent',hax(2));
quiver(hax(2),0,0,cos(yangle)*l,sin(yangle)*l,0,'w-');
quiver(hax(2),0,0,cos(yangle)*l,sin(yangle)*l,0,'k--');
text(cos(yangle)*l,sin(yangle)*l,'y','parent',hax(2));

axis(hax(2),'image');
axis(hax(2),'xy');  % This means image will be upside-sown from how Matlab normally displays images
title(hax(2),'Registered bkgd image');
colormap(hfig,'jet');

try
  if exist(imsavename,'file'),
    try
      delete(imsavename);
    catch ME,
      warning('Could not delete file %s:\n %s',imsavename,getReport(ME));
    end
  end
  set(hfig,'Units','pixels','Position',figpos);
  save2png(imsavename,hfig);
catch ME,
  warning('Could not save to %s:\n%s',imsavename,getReport(ME));
end

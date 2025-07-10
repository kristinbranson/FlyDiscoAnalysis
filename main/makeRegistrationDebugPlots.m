function makeRegistrationDebugPlots(isInDebugMode, iscircle, bkgdImage, circleRadius_mm, pairDist_mm, originX, originY, offTheta, ...
                                    circleCenterX, circleCenterY, circleRadius, registrationPoints, bowlMarkerPoints, ...
                                    bowlMarkerPairTheta_true, markerPairAngle_true, doesYAxisPointUp)

if isInDebugMode,
  hfig = figure() ;
  set(hfig,'Position',[20,20,1500,600]);
  clf;
  if iscircle,
    nsubplots = 3;
  else
    nsubplots = 4;
  end
  hax = createsubplots(1,nsubplots,0.05);
  %hax(1) = subplot(1,nsubplots,1);
  axes(hax(1)); 
  imagesc(bkgdImage);
  hold on;
  if iscircle,
    l = circleRadius_mm/2 / scale;
  else
    l = pairDist_mm/4 / scale;
  end
  xangle = 0;
  yangle = pi/2;
  quiver(originX,originY,cos(xangle-offTheta)*l,sin(xangle-offTheta)*l,0,'k--');
  text(originX+cos(xangle-offTheta)*l,originY+sin(xangle-offTheta)*l,'x');
  quiver(originX,originY,cos(yangle-offTheta)*l,sin(yangle-offTheta)*l,0,'k-');
  text(originX+cos(yangle-offTheta)*l,originY+sin(yangle-offTheta)*l,'y');
  if iscircle,
    tmp = linspace(0,2*pi,50);
    plot(circleCenterX + circleRadius*cos(tmp),circleCenterY + circleRadius*sin(tmp),'k-');
  else
    plot(registrationPoints(1,:),registrationPoints(2,:),'ks');
  end
  if nBowlMarkers > 0,
    plot(bowlMarkerPoints(1,:),bowlMarkerPoints(2,:),'mo');
  end
  axis image ; 
  if doesYAxisPointUp ,
    axis xy ;  % This means image will be upside-down from how Matlab normally displays images
  else
    axis ij ;
  end
  axes(hax(2)); 
  %hax(2) = subplot(1,nsubplots,2);
  imagesc(bkgdImage);
  hold on;
  quiver(originX,originY,cos(xangle-offTheta)*l,sin(xangle-offTheta)*l,0,'k--');
  text(originX+cos(xangle-offTheta)*l,originY+sin(xangle-offTheta)*l,'x');
  quiver(originX,originY,cos(yangle-offTheta)*l,sin(yangle-offTheta)*l,0,'k-');
  text(originX+cos(yangle-offTheta)*l,originY+sin(yangle-offTheta)*l,'y');
  if iscircle,
    tmp = linspace(0,2*pi,50);
    plot(circleCenterX + circleRadius*cos(tmp),circleCenterY + circleRadius*sin(tmp),'k-');
  else
    scatter(registrationPoints(1,:),registrationPoints(2,:),50,1:size(registrationPoints,2),'s');
  end
  if nBowlMarkers > 0,
    plot(bowlMarkerPoints(1,:),bowlMarkerPoints(2,:),'mo');
  end
  axis image ;
  if doesYAxisPointUp ,
    axis xy ;  % This means image will be upside-down from how Matlab normally displays images
  else
    axis ij ;
  end
  
  axes(hax(3)); 
  %hax(3) = subplot(1,nsubplots,3);
  imagesc(circleim);
  hold on;
  quiver(originX,originY,cos(xangle-offTheta)*l,sin(xangle-offTheta)*l,0,'w--');
  text(originX+cos(xangle-offTheta)*l,originY+sin(xangle-offTheta)*l,'x','color','w');
  quiver(originX,originY,cos(yangle-offTheta)*l,sin(yangle-offTheta)*l,0,'w-');
  text(originX+cos(yangle-offTheta)*l,originY+sin(yangle-offTheta)*l,'y','color','w');
  if iscircle,
    tmp = linspace(0,2*pi,50);
    plot(circleCenterX + circleRadius*cos(tmp),circleCenterY + circleRadius*sin(tmp),'k-');
  else
    scatter(registrationPoints(1,:),registrationPoints(2,:),50,1:size(registrationPoints,2),'s');
  end
  if nBowlMarkers > 0,
    plot(bowlMarkerPoints(1,:),bowlMarkerPoints(2,:),'mo');
  end
  axis image ;
  if doesYAxisPointUp ,
    axis xy ;  % This means image will be upside-down from how Matlab normally displays images
  else
    axis ij ;
  end  
  linkaxes(hax(1:3));
  if ~iscircle
    axes(hax(4)); 
    %subplot(1,nsubplots,4);
    tmp = linspace(0,2*pi,100);
    plot(cos(tmp)*pairDist_mm/2,sin(tmp)*pairDist_mm/2,'b');
    hold on;
    tmp = pi/2:pi/2:2*pi;
    plot([cos(bowlMarkerPairTheta_true+markerPairAngle_true/2+tmp)*pairDist_mm/2;zeros(size(tmp))],...
      [sin(bowlMarkerPairTheta_true+markerPairAngle_true/2+tmp)*pairDist_mm/2;zeros(size(tmp))],'g-');
    plot([cos(bowlMarkerPairTheta_true-markerPairAngle_true/2+tmp)*pairDist_mm/2;zeros(size(tmp))],...
      [sin(bowlMarkerPairTheta_true-markerPairAngle_true/2+tmp)*pairDist_mm/2;zeros(size(tmp))],'g-');
    plot([0,cos(bowlMarkerPairTheta_true)*pairDist_mm/2],[0,sin(bowlMarkerPairTheta_true)*pairDist_mm/2],'m-');
    [x,y] = registerfn(registrationPoints(1,:),registrationPoints(2,:));
    plot(x,y,'ks');
    [x,y] = registerfn(bowlMarkerPoint(1),bowlMarkerPoint(2));
    plot(x,y,'mo');
    quiver(0,0,cos(xangle)*pairDist_mm/4,sin(xangle)*pairDist_mm/4,0,'k');
    text(cos(xangle)*pairDist_mm/4,sin(xangle)*pairDist_mm/4,'x');
    quiver(0,0,cos(yangle)*pairDist_mm/4,sin(yangle)*pairDist_mm/4,0,'k');
    text(cos(yangle)*pairDist_mm/4,sin(yangle)*pairDist_mm/4,'y');
    quiver(0,0,0,pairDist_mm/4,0,'k');
    axis equal;
    axisalmosttight();
  end  
end

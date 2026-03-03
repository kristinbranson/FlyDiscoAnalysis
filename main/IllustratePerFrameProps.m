function IllustratePerFrameProps(trx,fn,varargin)

[hax,hfig,fly,t] = myparse(varargin,'hax',[],'hfig',[],'fly',1,'t',1);

if isempty(hax) && isempty(hfig),
  hfig = figure;
  hax = gca;
elseif isempty(hax),
  figure(hfig);
  hax = gca;
elseif isempty(hfig),
  hfig = get(hax,'Parent');
end
cla(hax,'reset');

switch fn
  
  case {'corfrac_min','corfrac_maj','du_cor','dv_cor','absdv_cor'},
    corfrac_min = trx(fly).corfrac_min(t);
    corfrac_maj = trx(fly).corfrac_maj(t);
    z = zeros(1,trx(fly).nframes-1);
    o = ones(1,trx(fly).nframes-1);
        
    [x1,y1,x2,y2] = rfrac2center(trx,fly,[corfrac_maj;corfrac_min]);
    [xhead1,yhead1,xhead2,yhead2] = rfrac2center(trx,fly,[o;z]);
    [xtail1,ytail1,xtail2,ytail2] = rfrac2center(trx,fly,[-o;z]);
    [xleft1,yleft1,xleft2,yleft2] = rfrac2center(trx,fly,[z;-o]);
    [xright1,yright1,xright2,yright2] = rfrac2center(trx,fly,[z;o]);
    theta_mm = trx(fly).theta_mm(t);
    dt = trx(fly).dt(t);
    du = trx(fly).du_cor(t)*dt;
    dv = trx(fly).dv_cor(t)*dt;
    
    hfly = nan(1,2);
    hfly(1) = trx.drawfly(fly,t,'registered',true,'shape','ellipse','color','k','parent',hax);
    hold(hax,'on');
    hfly(2) = trx.drawfly(fly,t+1,'registered',true,'shape','ellipse','color',[.7,.7,.7],'parent',hax);
    hhead = plot(hax,[xhead1(t),xhead2(t)],[yhead1(t),yhead2(t)],'.-','color','b');
    htail = plot(hax,[xtail1(t),xtail2(t)],[ytail1(t),ytail2(t)],'.-','color','k');
    hleft = plot(hax,[xleft1(t),xleft2(t)],[yleft1(t),yleft2(t)],'.-','color','k');
    hright = plot(hax,[xright1(t),xright2(t)],[yright1(t),yright2(t)],'.-','color','k');
    %hcor = plot(hax,[x1(t),x2(t)],[y1(t),y2(t)],'-','color',[.7,0,0],'linewidth',3);
    hcor2 = plot(hax,x2(t),y2(t),'o','color',[.7,.7,.7],'markerfacecolor',[.7,.7,.7]);
    hcor1 = plot(hax,x1(t),y1(t),'ko','markerfacecolor','k');
    htexthead = text(mean([xhead1(t),xhead2(t)]),mean([yhead1(t),yhead2(t)]),'head','horizontalalignment','center','parent',hax);
    hcor = text(mean([x1(t),x2(t)]),mean([y1(t),y2(t)]),'center of rotation','horizontalalignment','center','parent',hax);
    hdu = quiver(x1(t),y1(t),du*cos(theta_mm),du*sin(theta_mm),0,'parent',hax);
    hdv = quiver(x1(t),y1(t),-dv*sin(theta_mm),dv*cos(theta_mm),0,'parent',hax);
    dx = du*cos(theta_mm)-dv*sin(theta_mm);
    dy = du*sin(theta_mm)+dv*cos(theta_mm);
    hdudv = quiver(x1(t),y1(t),dx,dy,0,'color',[.7,0,0],'linewidth',2,'parent',hax);
    axis(hax,'equal');
    axisalmosttight(.1,hax);
    title(hax,sprintf('cor, fly %d, t = %d',fly,t-trx(fly).off));
    
  case {'angle2wall','dist2wall'},
    
    x = trx(fly).x_mm(t);
    y = trx(fly).y_mm(t);
    theta_mm = trx(fly).theta_mm(t);
    angle2wall = trx(fly).angle2wall(t);
    dist2wall = trx(fly).dist2wall(t);
    theta = linspace(0,2*pi,100);
    hfly = trx.drawfly(fly,t,'registered',true,'shape','triangle','color','k','parent',hax);
    hold(hax,'on');
    quiver(x,y,cos(theta_mm),sin(theta_mm),'parent',hax);
    quiver(x,y,dist2wall*cos(theta_mm+angle2wall),dist2wall*sin(theta_mm+angle2wall),0,'parent',hax);
    axisalmosttight(.2,hax);
    axis(hax,'equal');
    ax = axis(hax);
    harena = plot(hax,trx.landmark_params.arena_center_mm(1)+trx.landmark_params.arena_radius_mm*cos(theta),...
      trx.landmark_params.arena_center_mm(2)+trx.landmark_params.arena_radius_mm*sin(theta),'k--');
    axis(hax,ax);
    title(hax,sprintf('angle2wall, fly %d, t = %d',fly,t-trx(fly).off));
    
  case {'dcenter','anglefrom1to2_center'},
    
    fly2 = trx(fly).closestfly_center(t);
    off = trx(fly2).off - trx(fly).off;
    x = trx(fly).x_mm(t);
    y = trx(fly).y_mm(t);
    dcenter = trx(fly).dcenter(t);
    theta_mm = trx(fly).theta_mm(t);
    anglefrom1to2 = trx(fly).anglefrom1to2_center(t);
    hfly1 = trx.drawfly(fly,t,'registered',true,'shape','triangle','color','r','parent',hax);
    hold(hax,'on');
    hfly2 = trx.drawfly(fly2,t+off,'registered',true,'shape','ellipse','color','k','parent',hax);
    quiver(x,y,cos(theta_mm),sin(theta_mm),'parent',hax);
    quiver(x,y,dcenter*cos(theta_mm+anglefrom1to2),dcenter*sin(theta_mm+anglefrom1to2),0,'parent',hax);
    axis(hax,'equal');
    axisalmosttight(1,hax);
    flies = trx.exp2flies{trx.fly2exp(fly)};
    for fly3 = flies,
      if fly3 == fly || fly3 == fly2, continue; end
      off = trx(fly3).off - trx(fly).off;
      if t+off < 1 || t+off > trx(fly3).nframes,
        continue;
      end
      trx.drawfly(fly3,t+off,'registered',true,'shape','ellipse','color',[.5,.5,.5],'parent',hax);
      hold(hax,'on');
    end

    title(hax,sprintf('dcenter, fly %d, t = %d',fly,t-trx(fly).off));
    
  case 'dnose2ell',
    
    fly2 = trx(fly).closestfly_nose2ell(t);
    off = trx(fly2).off - trx(fly).off;
    x = trx(fly).x_mm(t);
    y = trx(fly).y_mm(t);
    xnose = trx(fly).xnose_mm(t);
    ynose = trx(fly).ynose_mm(t);
    dnose2ell = trx(fly).dnose2ell(t);
    theta_mm = trx(fly).theta_mm(t);
    anglefrom1to2 = trx(fly).anglefrom1to2_nose2ell(t);
    hfly1 = trx.drawfly(fly,t,'registered',true,'shape','triangle','color','r','parent',hax);
    hold(hax,'on');
    hfly2 = trx.drawfly(fly2,t+off,'registered',true,'shape','ellipse','color','k','parent',hax);
    quiver(x,y,cos(theta_mm),sin(theta_mm),'parent',hax);
    quiver(xnose,ynose,dnose2ell*cos(theta_mm+anglefrom1to2),dnose2ell*sin(theta_mm+anglefrom1to2),0,'parent',hax);
    thetas = linspace(theta_mm+anglefrom1to2-pi/6,theta_mm+anglefrom1to2+pi/6,100);
    plot(hax,xnose+cos(thetas)*dnose2ell,ynose+sin(thetas)*dnose2ell);
    axis(hax,'equal');
    axisalmosttight(1,hax);
    flies = trx.exp2flies{trx.fly2exp(fly)};
    for fly3 = flies,
      if fly3 == fly || fly3 == fly2, continue; end
      off = trx(fly3).off - trx(fly).off;
      if t+off < 1 || t+off > trx(fly3).nframes,
        continue;
      end
      trx.drawfly(fly3,t+off,'registered',true,'shape','ellipse','color',[.5,.5,.5],'parent',hax);
      hold(hax,'on');
    end

    title(hax,sprintf('dnose2ell, fly %d, t = %d',fly,t-trx(fly).off));
    
    case 'dell2nose',
    
      fly2 = trx(fly).closestfly_nose2ell(t);
      off = trx(fly2).off - trx(fly).off;
      x = trx(fly).x_mm(t);
      y = trx(fly).y_mm(t);
      xnose = trx(fly2).xnose_mm(t+off);
      ynose = trx(fly2).ynose_mm(t+off);
      dell2nose = trx(fly).dell2nose(t);
      theta_mm = trx(fly).theta_mm(t);
      anglefrom1to2 = trx(fly).anglefrom1to2_nose2ell(t);
      hfly1 = trx.drawfly(fly,t,'registered',true,'shape','ellipse','color','r','parent',hax);
      hold(hax,'on');
      hfly2 = trx.drawfly(fly2,t+off,'registered',true,'shape','ellipse','color','k','parent',hax);
      quiver(x,y,cos(theta_mm),sin(theta_mm),'parent',hax);
      quiver(xnose,ynose,dell2nose*cos(theta_mm+anglefrom1to2+pi),dell2nose*sin(theta_mm+anglefrom1to2+pi),0,'parent',hax);
      thetas = linspace(theta_mm+pi+anglefrom1to2-pi/6,theta_mm+pi+anglefrom1to2+pi/6,100);
      plot(hax,xnose+cos(thetas)*dell2nose,ynose+sin(thetas)*dell2nose);
      axis(hax,'equal');
      axisalmosttight(1,hax);
      flies = trx.exp2flies{trx.fly2exp(fly)};
      for fly3 = flies,
        if fly3 == fly || fly3 == fly2, continue; end
        off = trx(fly3).off - trx(fly).off;
        if t+off < 1 || t+off > trx(fly3).nframes,
          continue;
        end
        trx.drawfly(fly3,t+off,'registered',true,'shape','ellipse','color',[.5,.5,.5],'parent',hax);
        hold(hax,'on');
      end
      
      title(hax,sprintf('dell2nose, fly %d, t = %d',fly,t-trx(fly).off));
      
  case 'anglesub',
    
    fly2 = trx(fly).closestfly_anglesub(t);
    off = trx(fly2).off - trx(fly).off;
    x = trx(fly).x_mm(t);
    y = trx(fly).y_mm(t);
    x2 = trx(fly2).x_mm(t+off);
    y2 = trx(fly2).y_mm(t+off);
    xnose = trx(fly).xnose_mm(t);
    ynose = trx(fly).ynose_mm(t);
    
    dcenter = sqrt( (x-x2)^2 + (y-y2)^2 );
    theta_mm = trx(fly).theta_mm(t);
    anglefrom1to2 = trx(fly).anglefrom1to2_anglesub(t);
    anglesub = trx(fly).anglesub(t);
    
    hfly1 = trx.drawfly(fly,t,'registered',true,'shape','triangle','color','r','parent',hax);
    hold(hax,'on');
    hfly2 = trx.drawfly(fly2,t+off,'registered',true,'shape','ellipse','color','k','parent',hax);
    quiver(x,y,cos(theta_mm),sin(theta_mm),'parent',hax);
    quiver(xnose,ynose,1.25*dcenter*cos(theta_mm+anglefrom1to2),1.25*dcenter*sin(theta_mm+anglefrom1to2),0,'parent',hax);
    if modrange(anglesub/2+anglefrom1to2,-pi,pi) > trx.perframe_params.fov/2,
      angle2 = anglefrom1to2+trx.perframe_params.fov/2;
      angle1 = angle2 - anglesub;
    elseif modrange(anglefrom1to2-anglesub/2,-pi,pi) < -trx.perframe_params.fov/2,
      angle1 = anglefrom1to2-trx.perframe_params.fov/2;
      angle2 = angle1 + anglesub;
    else
      angle1 = anglefrom1to2-anglesub/2;
      angle2 = anglefrom1to2+anglesub/2;
    end
    quiver(xnose,ynose,1.25*dcenter*cos(theta_mm+angle1),1.25*dcenter*sin(theta_mm+angle1),0,'parent',hax);
    quiver(xnose,ynose,1.25*dcenter*cos(theta_mm+angle2),1.25*dcenter*sin(theta_mm+angle2),0,'parent',hax);
    thetas = linspace(theta_mm+angle1,theta_mm+angle2,100);
    plot(hax,xnose+cos(thetas)*dcenter,ynose+sin(thetas)*dcenter);
    axis(hax,'equal');
    axisalmosttight(1,hax);
    flies = trx.exp2flies{trx.fly2exp(fly)};
    for fly3 = flies,
      if fly3 == fly || fly3 == fly2, continue; end
      off = trx(fly3).off - trx(fly).off;
      if t+off < 1 || t+off > trx(fly3).nframes,
        continue;
      end
      trx.drawfly(fly3,t+off,'registered',true,'shape','ellipse','color',[.5,.5,.5],'parent',hax);
      hold(hax,'on');
    end
    
    title(hax,sprintf('anglesub, fly %d, t = %d',fly,t-trx(fly).off));
    
  case {'phi','yaw'},

    x = trx(fly).x_mm(t);
    y = trx(fly).y_mm(t);
    a = trx(fly).a_mm(t);
    phi = trx(fly).phi(t);
    theta_mm = trx(fly).theta_mm(t);
    yaw = trx(fly).yaw(t);
    hfly1 = trx.drawfly(fly,t,'registered',true,'shape','triangle','color','k','parent',hax);
    hold(hax,'on');
    hfly2 = trx.drawfly(fly,t+1,'registered',true,'shape','triangle','color',[.5,.5,.5],'parent',hax);
    
    quiver(x,y,cos(phi),sin(phi),4*a,'parent',hax);
    quiver(x,y,cos(theta_mm),sin(theta_mm),4*a,'parent',hax);
    angle1 = theta_mm;
    angle2 = theta_mm+yaw;
    thetas = linspace(angle1,angle2,100);
    plot(hax,x+cos(thetas)*4*a,y+sin(thetas)*4*a);

    axis(hax,'equal');
    axisalmosttight(.1,hax);
    title(hax,sprintf('phi, fly %d, t = %d',fly,t-trx(fly).off));

  case {'magveldiff','magveldiff_center','magveldiff_nose2ell','magveldiff_ell2nose','magveldiff_anglesub'},
    
    m = regexp(fn,'magveldiff_(.*)','tokens','once');
    if isempty(m),
      type = 'center';
    else
      type = m{1};
    end
    
    fly2 = trx(fly).(['closestfly_',type])(t);
    magveldiff = trx(fly).(['magveldiff_',type])(t);
    off = trx(fly2).off - trx(fly).off;
    x = trx(fly).x_mm(t);
    y = trx(fly).y_mm(t);
    x2 = trx(fly2).x_mm(t+off);
    y2 = trx(fly2).y_mm(t+off);
    theta_mm = trx(fly).theta_mm(t);
    dx1 = diff(trx(fly).x_mm(t:t+1));
    dy1 = diff(trx(fly).y_mm(t:t+1));
    dx2 = diff(trx(fly2).x_mm(t+off:t+off+1));
    dy2 = diff(trx(fly2).y_mm(t+off:t+off+1));
    angle = atan2(dy2-dy1,dx2-dx1);
    
    hfly1 = trx.drawfly(fly,t,'registered',true,'shape','triangle','color','r','parent',hax);
    hold(hax,'on');
    hfly2 = trx.drawfly(fly2,t+off,'registered',true,'shape','triangle','color','k','parent',hax);
    quiver(x,y,dx1,dy1,0,'parent',hax);
    quiver(x,y,dx2,dy2,0,'parent',hax);
    quiver(x2,y2,dx2,dy2,0,'parent',hax);
    quiver(x+dx1,y+dy1,cos(angle)*magveldiff,sin(angle)*magveldiff,0,'parent',hax);
    axis(hax,'equal');
    axisalmosttight(.1,hax);
    
    title(hax,sprintf('magveldiff, fly %d, t = %d',fly,t-trx(fly).off));
    
  case {'veltoward','veltoward_center','veltoward_nose2ell','veltoward_ell2nose','veltoward_anglesub'},
    
    m = regexp(fn,'veltoward_(.*)','tokens','once');
    if isempty(m),
      type = 'center';
    else
      type = m{1};
    end
    
    fly2 = trx(fly).(['closestfly_',type])(t);
    veltoward = trx(fly).(['veltoward_',type])(t);
    off = trx(fly2).off - trx(fly).off;
    x = trx(fly).x_mm(t);
    y = trx(fly).y_mm(t);
    x2 = trx(fly2).x_mm(t+off);
    y2 = trx(fly2).y_mm(t+off);
    angle = atan2(y2-y,x2-x);
    theta_mm = trx(fly).theta_mm(t);
    dx1 = diff(trx(fly).x_mm(t:t+1));
    dy1 = diff(trx(fly).y_mm(t:t+1));
    
    hfly1 = trx.drawfly(fly,t,'registered',true,'shape','triangle','color','r','parent',hax);
    hold(hax,'on');
    hfly2 = trx.drawfly(fly2,t+off,'registered',true,'shape','triangle','color','k','parent',hax);
    quiver(x,y,dx1,dy1,0,'parent',hax);
    quiver(x,y,cos(angle)*veltoward,sin(angle)*veltoward,0,'parent',hax);
    axis(hax,'equal');
    axisalmosttight(.1,hax);
    
    title(hax,sprintf('veltoward, fly %d, t = %d',fly,t-trx(fly).off));
    
end
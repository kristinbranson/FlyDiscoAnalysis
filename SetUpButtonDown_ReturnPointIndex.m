function hdata = SetUpButtonDown_ReturnPointIndex(hax,x,y,callback,maxdfrac)

hdata = struct;
hdata.hax = hax;
hdata.hchil = setdiff(findall(hax,'HitTest','on'),hax);
if nargin < 5,
  maxdfrac = 1/100;
end
hdata.maxdfrac = maxdfrac;
set(hax,'HitTest','on');
set(hdata.hchil,'HitTest','off');

if iscell(callback),
  args = callback(2:end);
  callback = callback{1};
else
  args = {};
end

set(hax,'ButtonDownFcn',@ButtonDownFcn_ReturnPointIndex);

  function ButtonDownFcn_ReturnPointIndex(hObject,eventdata,varargin)
    
    try
    
    ax = axis(hax);
    dx = diff(ax(1:2));
    dy = diff(ax(3:4));
    
    pt = get(hax,'CurrentPoint');
    xclick = pt(1,1);
    yclick = pt(1,2);
    [d,j] = min((x-xclick).^2/dx^2 + (y-yclick).^2/dy^2);
    d = sqrt(d);
    if d < maxdfrac,
      fprintf('%d: %f, at (%f, %f)\n',j,d,xclick,yclick);
      eventdata.pointindex = j;
      eventdata.distance = d;
      eventdata.xclick = xclick;
      eventdata.yclick = yclick;
      
      callback(hObject,eventdata,args{:});
    end

    catch ME,
      
      warning('Error in ButtonDownFcn, removing Callback: %s\n',getReport(ME));
      ClearButtonDown_ReturnPointIndex(hdata);
      
    end
    
  end
end
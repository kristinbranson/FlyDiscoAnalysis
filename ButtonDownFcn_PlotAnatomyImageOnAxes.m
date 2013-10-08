function ButtonDownFcn_PlotAnatomyImageOnAxes(hax,eventdata,line_names,xs,ys,imwidth,hdata,varargin)

line_name = line_names{eventdata.pointindex};
x = xs(eventdata.pointindex);
y = ys(eventdata.pointindex);
him = PlotAnatomyImageOnAxes(line_name,x,y,imwidth,hax,hdata,varargin{:});
ud = get(hax,'UserData');
if ~isstruct(ud),
  ud = struct;
end
htmp = findall(him,'HitTest','on');
set(htmp,'HitTest','off');
if isfield(ud,'hims_anatomy');
  ud.hims_anatomy(end+1) = him;
else
  ud.hims_anatomy = him;
end
set(hax,'UserData',ud);

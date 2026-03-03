function ButtonDownFcn_PlotCompartmentHistogramOnAxes(hax,eventdata,line_names,xs,ys,imwidth,imheight,hdata,int_manual,varargin)

line_name = line_names{eventdata.pointindex};
x = xs(eventdata.pointindex);
y = ys(eventdata.pointindex);
him = PlotCompartmentHistogramOnAxes(line_name,x,y,imwidth,imheight,hax,hdata,int_manual,varargin{:});
ud = get(hax,'UserData');
htmp = findall(him,'HitTest','on');
set(htmp,'HitTest','off');
if ~isstruct(ud),
  ud = struct;
end
if isfield(ud,'hhist_anatomy');
  ud.hhist_anatomy(end+1:end+numel(him)) = him;
else
  ud.hhist_anatomy = him;
end
set(hax,'UserData',ud);

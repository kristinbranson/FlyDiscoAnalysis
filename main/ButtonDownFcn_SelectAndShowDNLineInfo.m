function ButtonDownFcn_SelectAndShowDNLineInfo(hObject,eventdata,line_names,hpts,...
  isdnexpr,dnnames,varargin)

userdata = get(hObject,'UserData');
if ~isstruct(userdata),
  userdata = struct('isselected',false(1,numel(line_names)));
elseif ~isfield(userdata,'isselected') || numel(userdata.isselected) ~= numel(line_names),
  userdata.isselected = false(1,numel(line_names));
end

if ~isfield(userdata,'htext') || numel(userdata.htext) ~= numel(line_names),
  userdata.htext = nan(1,numel(line_names));
end

if strcmpi(eventdata.selectiontype,'alt'),
  if userdata.isselected(eventdata.pointindex),
    line_names_curr = line_names(userdata.isselected);
    ShowDNLineInfo(line_names_curr,isdnexpr(userdata.isselected,:),dnnames,varargin{:});
  end
  return;
end

if strcmpi(eventdata.selectiontype,'extend'),
  if userdata.isselected(eventdata.pointindex),
    line_names_curr = line_names(eventdata.pointindex);
    ShowDNLineInfo(line_names_curr,isdnexpr(userdata.isselected,:),dnnames,varargin{:});
  end
  return;
end


if ~strcmpi(eventdata.selectiontype,'normal'),
  return;
end

if ~userdata.isselected(eventdata.pointindex),
  color = get(hpts(eventdata.pointindex),'Color');
  userdata.oldmarkers{eventdata.pointindex} = get(hpts(eventdata.pointindex),'Marker');
  set(hpts(eventdata.pointindex),'Marker','o','MarkerFaceColor',color);

  shortlinename = line_names{eventdata.pointindex};
  shortlinename = regexprep(shortlinename,'GMR_','R');
  shortlinename = regexprep(shortlinename,'_AE_01','');
  shortlinename = regexprep(shortlinename,'_AD_01','D');
  if ishandle(userdata.htext(eventdata.pointindex)),
    set(userdata.htext(eventdata.pointindex),'String',[' ',shortlinename],'Color',color,'Visible','on',...
      'Position',[eventdata.xpoint,eventdata.ypoint],'HitTest','off',...
      'HorizontalAlignment','left','VerticalAlignment','middle');
  else
    userdata.htext(eventdata.pointindex) = text(eventdata.xpoint,eventdata.ypoint,[' ',shortlinename],...
      'Color',color,'Parent',hObject,'HitTest','off',...
      'HorizontalAlignment','left','VerticalAlignment','middle');
  end
  drawnow;
else
  set(hpts(eventdata.pointindex),'Marker',userdata.oldmarkers{eventdata.pointindex},'MarkerFaceColor','none');
  if ishandle(userdata.htext(eventdata.pointindex)),
    set(userdata.htext(eventdata.pointindex),'Visible','off');
  end
end
userdata.isselected(eventdata.pointindex) = ~userdata.isselected(eventdata.pointindex);
set(hObject,'UserData',userdata);

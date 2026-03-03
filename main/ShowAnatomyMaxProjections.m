function iselect = ShowAnatomyMaxProjections(handles,is)

persistent imdata;
if isempty(imdata),
  
end

hfig = 1000;
figure(hfig);
clf;
n = numel(is);
nc = ceil(sqrt(n));
nr = ceil(n/nc);
hax = createsubplots(nr,nc,[.01,.1]);
if numel(hax) > n,
  delete(hax(n+1:end));
end
set(hfig,'UserData',0);

hbutton = nan(1,n);
for ii = 1:n,
  
  i = is(ii);
  
  try
    im = imread(handles.imdata(i).maxproj_file_system_path);
    image(im,'Parent',hax(ii));
    axis(hax(ii),'image','off');
  catch ME,
    warning('Could not read %s: %s',handles.imdata(i).maxproj_file_system_path,getReport(ME));
  end
  
  capture_date = handles.imdata(i).capture_date;
  if numel(capture_date) > 10,
    capture_date = capture_date(1:10);
  end
  
  s = sprintf('ann=%d, view=%d, qi=%.2f, %s',handles.isannotated(i),handles.imdata(i).isviewed,handles.imdata(i).qi,capture_date);
%   hti = title(hax(ii),s);
%   set(hti,'Interpreter','none','Units','normalized');
%   pos = get(hti,'Extent');
%   delete(hti);
  
  axpos = get(hax(ii),'Position');
  
  pos = [axpos(1),axpos(2)-.07,axpos(3),.06];
  title(hax(ii),handles.imdata(i).line,'Interpreter','none');
  
  hbutton(ii) = uicontrol(hfig,'style','pushbutton','Units','normalized',...
    'Position',pos,'String',s,'Callback',{@ChooseImage_ButtonCallback,i});
    
end

set(hfig,'ToolBar','figure','Visible','on');

waitfor(hfig,'UserData',1);

if ishandle(hfig),
  set(hfig,'Visible','off');
end

  function ChooseImage_ButtonCallback(hObject,eventdata,i) %#ok<INUSL>
    
    iselect = i;
    set(hfig,'UserData',1);
    
  end

end
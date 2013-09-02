function ClearButtonDown_ReturnPointIndex(hdata)

hchil = hdata.hchil(ishandle(hdata.hchil));
set(hchil,'HitTest','on');
set(hdata.hax,'ButtonDownFcn','');

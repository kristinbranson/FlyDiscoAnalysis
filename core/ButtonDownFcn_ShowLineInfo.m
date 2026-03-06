function ButtonDownFcn_ShowLineInfo(hObject,eventdata,line_names,varargin)

line_names_curr = line_names(eventdata.pointindex);
ShowLineInfo(line_names_curr,varargin{:});
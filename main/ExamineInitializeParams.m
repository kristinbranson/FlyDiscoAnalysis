function varargout = ExamineInitializeParams(varargin)
% EXAMINEINITIALIZEPARAMS MATLAB code for ExamineInitializeParams.fig
%      EXAMINEINITIALIZEPARAMS, by itself, creates a new EXAMINEINITIALIZEPARAMS or raises the existing
%      singleton*.
%
%      H = EXAMINEINITIALIZEPARAMS returns the handle to a new EXAMINEINITIALIZEPARAMS or the handle to
%      the existing singleton*.
%
%      EXAMINEINITIALIZEPARAMS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EXAMINEINITIALIZEPARAMS.M with the given input arguments.
%
%      EXAMINEINITIALIZEPARAMS('Property','Value',...) creates a new EXAMINEINITIALIZEPARAMS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ExamineInitializeParams_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ExamineInitializeParams_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ExamineInitializeParams

% Last Modified by GUIDE v2.5 16-May-2011 13:51:51

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ExamineInitializeParams_OpeningFcn, ...
                   'gui_OutputFcn',  @ExamineInitializeParams_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before ExamineInitializeParams is made visible.
function ExamineInitializeParams_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ExamineInitializeParams (see VARARGIN)

% Choose default command line output for ExamineInitializeParams
handles.output = hObject;
handles.DATASOURCE_SAGE = 1;
handles.DATASOURCE_LOAD = 2;

dateranges = varargin{1};
daterange_idx = varargin{2};
username = varargin{3};
loadcacheddata = varargin{4};
datafilename = varargin{5};

% set default values

% username
set(handles.edit_username,'String',username);

% data source
if isempty(loadcacheddata) || ~loadcacheddata,
  set(handles.popupmenu_datasource,'Value',handles.DATASOURCE_SAGE);
else
  set(handles.popupmenu_datasource,'Value',handles.DATASOURCE_LOAD);
end

UpdateControlVisibility(handles);

% date range
set(handles.popupmenu_daterange,'String',dateranges,'Value',daterange_idx);

% data file name
set(handles.edit_datafilename,'String',datafilename,'Enable','inactive');

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ExamineInitializeParams wait for user response (see UIRESUME)
uiwait(handles.figure1);

function UpdateControlVisibility(handles)

if get(handles.popupmenu_datasource,'Value') == handles.DATASOURCE_SAGE,
  set(handles.text_datafilename,'Visible','off');
  set(handles.edit_datafilename,'Visible','off');
  set(handles.text_daterange,'Visible','on');
  set(handles.popupmenu_daterange,'Visible','on');
else
  set(handles.popupmenu_datasource,'Value',handles.DATASOURCE_LOAD);
  set(handles.text_datafilename,'Visible','on');
  set(handles.edit_datafilename,'Visible','on');
  set(handles.text_daterange,'Visible','off');
  set(handles.popupmenu_daterange,'Visible','off');
end

% --- Outputs from this function are returned to the command line.
function varargout = ExamineInitializeParams_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.success;
varargout{2} = get(handles.popupmenu_daterange,'Value');
varargout{3} = get(handles.edit_username,'String');
if get(handles.popupmenu_datasource,'Value') == handles.DATASOURCE_LOAD,
  varargout{4} = get(handles.edit_datafilename,'String');
else
  varargout{4} = '';
end
delete(handles.figure1);

function edit_username_Callback(hObject, eventdata, handles)
% hObject    handle to edit_username (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_username as text
%        str2double(get(hObject,'String')) returns contents of edit_username as a double

% --- Executes during object creation, after setting all properties.
function edit_username_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_username (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_daterange.
function popupmenu_daterange_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_daterange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_daterange contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_daterange

% --- Executes during object creation, after setting all properties.
function popupmenu_daterange_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_daterange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
s = get(handles.edit_datafilename,'String');
if get(handles.popupmenu_datasource,'Value') == handles.DATASOURCE_LOAD && ...
    ( isempty(s) || ~exist(s,'file') ),
  uiwait(warndlg(sprintf('Invalid data file name >%s<',s)));
  return;
end
handles.success = true;
guidata(hObject,handles);
uiresume(handles.figure1);

% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.success = false;
guidata(hObject,handles);
uiresume(handles.figure1);


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
handles.success = false;
guidata(hObject,handles);
uiresume(handles.figure1);


% --- Executes on selection change in popupmenu_datasource.
function popupmenu_datasource_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_datasource (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_datasource contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_datasource
UpdateControlVisibility(handles);

% --- Executes during object creation, after setting all properties.
function popupmenu_datasource_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_datasource (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_datafilename_Callback(hObject, eventdata, handles)
% hObject    handle to edit_datafilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_datafilename as text
%        str2double(get(hObject,'String')) returns contents of edit_datafilename as a double


% --- Executes during object creation, after setting all properties.
function edit_datafilename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_datafilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over edit_datafilename.
function edit_datafilename_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to edit_datafilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename,pathname] = uigetfile('*.mat','MAT file containing cached data',get(hObject,'String'));
if ~ischar(filename),
  return;
end
set(hObject,'String',fullfile(pathname,filename));

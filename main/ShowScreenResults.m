function varargout = ShowScreenResults(varargin)
% SHOWSCREENRESULTS MATLAB code for ShowScreenResults.fig
%      SHOWSCREENRESULTS, by itself, creates a new SHOWSCREENRESULTS or raises the existing
%      singleton*.
%
%      H = SHOWSCREENRESULTS returns the handle to a new SHOWSCREENRESULTS or the handle to
%      the existing singleton*.
%
%      SHOWSCREENRESULTS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SHOWSCREENRESULTS.M with the given input arguments.
%
%      SHOWSCREENRESULTS('Property','Value',...) creates a new SHOWSCREENRESULTS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ShowScreenResults_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ShowScreenResults_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ShowScreenResults

% Last Modified by GUIDE v2.5 02-Oct-2012 15:43:09

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ShowScreenResults_OpeningFcn, ...
                   'gui_OutputFcn',  @ShowScreenResults_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1}) && exist(varargin{1}),
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before ShowScreenResults is made visible.
function ShowScreenResults_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ShowScreenResults (see VARARGIN)


% parse input

collecteddatamatfile = '/groups/branson/home/bransonk/behavioranalysis/code/Jdetect/Jdetect/experiments/gal4screen/ExperimentsFracTimePerSet20120906T160000.mat';
normalize_datatype = 'control';
normalize_unit = 'experiment';

[collecteddatamatfile,handles.normalize_datatype,...
  handles.normalize_unit] = myparse(varargin,...
  'collecteddatamatfile',collecteddatamatfile,...
  'normalize_datatype',normalize_datatype,...
  'normalize_unit',normalize_unit);

% load in data
handles = InitializeData(handles,collecteddatamatfile);

% normalization data
handles = UpdateNormalizationData(handles);

% initialize gui
handles = InitializeGUI(handles);

% normalize
handles = UpdatePlotData(handles);

% update gui
handles = UpdateGUI(handles);

% Choose default command line output for ShowScreenResults
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ShowScreenResults wait for user response (see UIRESUME)
% uiwait(handles.figure_ScreenResults);

function handles = InitializeData(handles,collecteddatamatfile)

if ischar(collecteddatamatfile),
  tmp = load(collecteddatamatfile);
else
  tmp = collecteddatamatfile;
end
fns = fieldnames(tmp);
for i = 1:numel(fns),
  handles.(fns{i}) = tmp.(fns{i});
end
handles.nbehaviors = size(handles.fractime,2);
handles.behaviors = handles.scoresfilestr(:,1);

[handles.linenames,~,handles.lineidx_perset] = unique(handles.linename_perset);
[~,handles.lineidx_perexp] = ismember({handles.metadata.line_name},handles.linenames);
handles.nlines = numel(handles.linenames);
handles.mean_fractime_perline = nan(handles.nlines,handles.nbehaviors);
for linei = 1:handles.nlines,
  handles.mean_fractime_perline(linei,:) = nanmean(handles.mean_fractime_perset(handles.lineidx_perset==linei,:),1);
end

function handles = UpdateGUI(handles)

set(handles.hexp,'XData',handles.rank(handles.lineidx_perexp),...
  'YData',handles.expdata);
for linei = 1:handles.nlines,
  idx = handles.lineidx_perset == linei;
  set(handles.hset(linei),'XData',handles.rank(linei)+zeros(1,nnz(idx)),...
    'YData',handles.setdata(idx));
  set(handles.hline(linei),'XData',handles.rank(linei),...
    'YData',handles.linedata(linei));
end
  


function handles = InitializeGUI(handles)

handles.showstati = 1;
set(handles.popupmenu_ShowStat,'String',handles.behaviors,'Value',handles.showstati);

handles.normalizeby = [];
tmp = handles.behaviors;
tmp{end+1} = 'None';
set(handles.listbox_NormalizeBy,'String',tmp,'Value',handles.normalizeby);

handles.colors = jet(handles.nlines);
handles.hexp = plot(nan,nan,'k.');
hold(handles.axes_ScreenResults,'on');

handles.hset = nan(1,handles.nlines);
for linei = 1:handles.nlines,
  handles.hset(linei) = plot(nan,nan,'*','Color',handles.colors(linei,:)*.7+.3);
end

handles.hline = nan(1,handles.nlines);
for linei = 1:handles.nlines,
  handles.hline(linei) = plot(nan,nan,'o','Color',handles.colors(linei,:)*.7,...
    'MarkerFaceColor',handles.colors(linei,:)*.7);
  set(handles.hline(linei),'ButtonDownFcn',...
    @(h,evt) LineButtonDownFcn(h,evt,linei,guidata(h)));
end

set(handles.axes_ScreenResults,'XLim',[-5,handles.nlines+4],'YLim',[0,max(handles.fractime(:))]);
set(handles.pushbutton_UpdatePlot,'Enable','off');

handles.selectedlinei = [];

if strcmpi(handles.normalize_datatype,'Control'),
  set(handles.menu_Norm_Control,'Checked','on');
  set(handles.menu_Norm_GAL4,'Checked','off');
  set(handles.menu_Norm_Line,'Enable','off');
else
  set(handles.menu_Norm_Control,'Checked','off');
  set(handles.menu_Norm_GAL4,'Checked','on');
  set(handles.menu_Norm_Line,'Enable','on');
end

set([handles.menu_Norm_Experiment,handles.menu_Norm_Set,handles.menu_Norm_Line],'Checked','off');
if strcmpi(handles.normalize_unit,'Experiment'),
  set(handles.menu_Norm_Experiment,'Checked','on');
elseif strcmpi(handles.normalize_unit,'Set'),
  set(handles.menu_Norm_Set,'Checked','on');
else
  set(handles.menu_Norm_Line,'Checked','on');
end

function handles = UpdateNormalizationData(handles)

% normalize by control data
if strcmpi(handles.normalize_datatype,'control'),
  
  if strcmpi(handles.normalize_unit,'set'),
    
    idx = strcmp(handles.linename_perset,'pBDPGAL4U');
    handles.X = handles.mean_fractime_perset(idx,:);
    
  else
    
    idx = strcmp({handles.metadata.line_name},'pBDPGAL4U');
    handles.X = handles.fractime(idx,:);
    
  end
  
  % normalize by all GMRs
elseif strcmpi(handles.normalize_datatype,'gmr'),
  
  if strcmpi(handles.normalize_unit,'line'),
    
    handles.X = handles.mean_fractime_perline;
    
  elseif strcmpi(handles.normalize_unit,'set'),
    
    handles.X = handles.mean_fractime_perset;
    
  elseif strcmpi(handles.normalize_unit,'experiment'),
    
    handles.X = handles.fractime;
    
  end
  
end

function handles = UpdatePlotData(handles)

if isempty(handles.normalizeby),
  handles.expdata = handles.fractime(:,handles.showstati);
  handles.setdata = handles.mean_fractime_perset(:,handles.showstati);
  handles.linedata = handles.mean_fractime_perline(:,handles.showstati);
  s = 'Data not normalized';
else
  [handles.normfcn,handles.normcoeffs] = GetBehaviorNormalizationFunction(handles.X(:,handles.normalizeby),handles.X(:,handles.showstati));
  handles.expdata = handles.normfcn(handles.fractime(:,handles.showstati),handles.fractime(:,handles.normalizeby));
  handles.setdata = handles.normfcn(handles.mean_fractime_perset(:,handles.showstati),handles.mean_fractime_perset(:,handles.normalizeby));
  handles.linedata = handles.normfcn(handles.mean_fractime_perline(:,handles.showstati),handles.mean_fractime_perline(:,handles.normalizeby));
  s = {};
  s{end+1} = sprintf('%s = %f',handles.behaviors{handles.showstati},handles.normcoeffs(1));
  for i = 1:numel(handles.normalizeby),
    if handles.normcoeffs(i+1) == 0,
      continue;
    elseif handles.normcoeffs(i+1) < 0,
      c = '-';
    else
      c = '+';
    end
    s{end+1} = sprintf(' %s %f %s',c,abs(handles.normcoeffs(i+1)),handles.behaviors{handles.normalizeby(i)});
  end
end

[~,handles.order] = sort(handles.linedata);
[~,handles.rank] = sort(handles.order);
handles.maxyvalue = max(handles.expdata);
handles.minyvalue = min(handles.expdata);
set(handles.axes_ScreenResults,'YLim',[handles.minyvalue,handles.maxyvalue]);
set(handles.text_NormalizationInfo,'String',s);

UpdateSelectionInfo(handles);

% --- Outputs from this function are returned to the command line.
function varargout = ShowScreenResults_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in listbox_NormalizeBy.
function listbox_NormalizeBy_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_NormalizeBy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_NormalizeBy contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_NormalizeBy

v = get(hObject,'Value');
if any(v > handles.nbehaviors),
  handles.normalizeby = [];
  set(hObject,'Value',[]);
else
  handles.normalizeby = v;
end
guidata(hObject,handles);
set(handles.pushbutton_UpdatePlot,'Enable','on');

% --- Executes during object creation, after setting all properties.
function listbox_NormalizeBy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_NormalizeBy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_ShowStat.
function popupmenu_ShowStat_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_ShowStat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_ShowStat contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_ShowStat

handles.showstati = get(hObject,'Value');
guidata(hObject,handles);
set(handles.pushbutton_UpdatePlot,'Enable','on');


% --- Executes during object creation, after setting all properties.
function popupmenu_ShowStat_CreateFcn(hObject, eventdata, linei, handles)
% hObject    handle to popupmenu_ShowStat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_UpdatePlot.
function pushbutton_UpdatePlot_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_UpdatePlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = UpdatePlotData(handles);
handles = UpdateGUI(handles);
guidata(hObject,handles);
set(hObject,'Enable','off');

% --- Executes on mouse press over axes background.
function axes_ScreenResults_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes_ScreenResults (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function UpdateSelectionInfo(handles)

linei = handles.selectedlinei;
s = handles.linenames(linei);
s{end+1} = sprintf('Rank=%d',numel(handles.linenames)-handles.rank(linei));
s{end+1} = sprintf('%s = %f',handles.behaviors{handles.showstati},handles.mean_fractime_perline(linei,handles.showstati));
if ~isempty(handles.normalizeby),
  snormalizeby = sprintf('%s,',handles.behaviors{handles.normalizeby});
  snormalizeby = snormalizeby(1:end-1);
  s{end+1} = sprintf('%s | %s = %+f',handles.behaviors{handles.showstati},snormalizeby,handles.linedata(linei));
end
set(handles.text_SelectionInfo,'String',s);

function LineButtonDownFcn(hObject,eventdata,linei,handles)

if ~isempty(handles.selectedlinei),
  set(handles.hline(handles.selectedlinei),'MarkerSize',6,...
    'Color',handles.colors(handles.selectedlinei,:),...
    'MarkerFaceColor',handles.colors(handles.selectedlinei,:));
end
set(handles.hline(linei),'MarkerSize',12,'MarkerFaceColor','k','Color','w');
handles.selectedlinei = linei;

UpdateSelectionInfo(handles);

selectiontype = get(handles.figure_ScreenResults,'SelectionType');
if strcmp(selectiontype,'open'),
  ShowLine(handles);

end
guidata(hObject,handles);

function ShowLine(handles)

linei = handles.selectedlinei;
if isempty(linei),
  warndlg('No line selected');
  return;
end

line_name = handles.linenames{linei};
% create an html page with links to all experiments
filenamecurr = fullfile(tempdir,[line_name,'.html']);
fid = fopen(filenamecurr,'w');
if fid >= 0,
  fprintf(fid,'<html>\n<title>%s</title>\n<body>\n',line_name);
  fprintf(fid,'<h1>%s</h1>\n',line_name);
  fprintf(fid,'<ul>\n');
  
  if ~isempty(handles.normalizeby),
    snormalizeby = sprintf('%s,',handles.behaviors{handles.normalizeby});
    snormalizeby = snormalizeby(1:end-1);
  end
  
  idx = find(handles.lineidx_perexp == linei);
  [~,order] = sort(-handles.expdata(idx));
  idx = idx(order);
  for i = idx,
    
    s = sprintf('%s = %f',handles.behaviors{handles.showstati},handles.fractime(i,handles.showstati));
    if ~isempty(handles.normalizeby),
      s = [s,sprintf(', %s | %s = %+f',handles.behaviors{handles.showstati},snormalizeby,handles.expdata(i))];
    end
    
    expdir = handles.metadata(i).file_system_path;
    [~,name] = fileparts(expdir);
    fprintf(fid,'  <li><a href="file://%s">%s</a>: %s</li>\n',expdir,name,s);
  end
  fprintf(fid,'</ul>\n');
  fprintf(fid,'</body>\n</html>\n');
  fclose(fid);
end
if ~exist(filenamecurr,'file'),
  warning('Could not open temporary file %s',filenamecurr);
  return;
end
web(filenamecurr,'-browser');

% --------------------------------------------------------------------
function menu_Callback(hObject, eventdata, handles)
% hObject    handle to menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_Open_Callback(hObject, eventdata, handles)
% hObject    handle to menu_Open (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ShowLine(handles);

% --------------------------------------------------------------------
function menu_NormData_Callback(hObject, eventdata, handles)
% hObject    handle to menu_NormData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function menu_Norm_Control_Callback(hObject, eventdata, handles)
% hObject    handle to menu_Norm_Control (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.normalize_datatype = 'control';
if strcmpi(handles.normalize_unit,'line'),
  handles.normalize_unit = 'set';
end
set(hObject,'Checked','on');
set(handles.menu_Norm_GAL4,'Checked','off');

handles = UpdateNormalizationData(handles);
handles = UpdatePlotData(handles);
handles = UpdateGUI(handles);
guidata(hObject,handles);

% --------------------------------------------------------------------
function menu_Norm_GAL4_Callback(hObject, eventdata, handles)
% hObject    handle to menu_Norm_GAL4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(hObject,'Checked','on');
set(handles.menu_Norm_Control,'Checked','off');

handles.normalize_datatype = 'gmr';
handles = UpdateNormalizationData(handles);
handles = UpdatePlotData(handles);
handles = UpdateGUI(handles);
guidata(hObject,handles);


% --------------------------------------------------------------------
function menu_NormUnit_Callback(hObject, eventdata, handles)
% hObject    handle to menu_NormUnit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_Norm_Experiment_Callback(hObject, eventdata, handles)
% hObject    handle to menu_Norm_Experiment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(setdiff([handles.menu_Norm_Experiment,handles.menu_Norm_Set,handles.menu_Norm_Line],hObject),...
  'Checked','off');
set(hObject,'Checked','on');

handles.normalize_unit = 'experiment';
handles = UpdateNormalizationData(handles);
handles = UpdatePlotData(handles);
handles = UpdateGUI(handles);
guidata(hObject,handles);

% --------------------------------------------------------------------
function menu_Norm_Set_Callback(hObject, eventdata, handles)
% hObject    handle to menu_Norm_Set (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(setdiff([handles.menu_Norm_Experiment,handles.menu_Norm_Set,handles.menu_Norm_Line],hObject),...
  'Checked','off');
set(hObject,'Checked','on');

handles.normalize_unit = 'set';
handles = UpdateNormalizationData(handles);
handles = UpdatePlotData(handles);
handles = UpdateGUI(handles);
guidata(hObject,handles);


% --------------------------------------------------------------------
function menu_Norm_Line_Callback(hObject, eventdata, handles)
% hObject    handle to menu_Norm_Line (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(setdiff([handles.menu_Norm_Experiment,handles.menu_Norm_Set,handles.menu_Norm_Line],hObject),...
  'Checked','off');
set(hObject,'Checked','on');

handles.normalize_unit = 'line';
handles = UpdateNormalizationData(handles);
handles = UpdatePlotData(handles);
handles = UpdateGUI(handles);
guidata(hObject,handles);

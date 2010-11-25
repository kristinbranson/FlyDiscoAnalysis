function varargout = BgModelReview(varargin)
% BGMODELREVIEW MATLAB code for BgModelReview.fig
%      BGMODELREVIEW, by itself, creates a new BGMODELREVIEW or raises the existing
%      singleton*.
%
%      H = BGMODELREVIEW returns the handle to a new BGMODELREVIEW or the handle to
%      the existing singleton*.
%
%      BGMODELREVIEW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BGMODELREVIEW.M with the given input arguments.
%
%      BGMODELREVIEW('Property','Value',...) creates a new BGMODELREVIEW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before BgModelReview_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to BgModelReview_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help BgModelReview

% Last Modified by GUIDE v2.5 13-Nov-2010 14:39:30

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @BgModelReview_OpeningFcn, ...
                   'gui_OutputFcn',  @BgModelReview_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
warning off MATLAB:str2func:invalidFunctionName;
if nargin && ischar(varargin{1}) && exist(varargin{1},'file'),
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before BgModelReview is made visible.
function BgModelReview_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to BgModelReview (see VARARGIN)

% Choose default command line output for BgModelReview
handles.output = hObject;

%% parameters

% for now, the ann and mat files are in 
handles.DEBUG = true;
handles.resultsdir = '/groups/branson/home/bransonk/tracking/data/olympiad/FlyBowl';

% name of movie within the experiment directory
handles.moviefilestr = 'movie.ufmf';
% name of annotation file
handles.annfilestr = 'movie.ufmf.ann';
% name of trajectory mat file
handles.trxfilestr = 'movie.mat';
% name of file to save results to
handles.savefilestr = 'CheckTrackingResults.txt';

% currently show first movie
handles.expdiri = 1;

% views of video to show

% variable in annotation to name in dropdown list
handles.var2ShowLabel = {'background_center','Background Center'
  'background_dev','Background Deviation'
  'background_median','Background Median'
  'background_mad','Background Median Absolute Deviation'
  'fracframesisback','Fraction of Frames in Bg Model'};
handles.annMenuShowLabels = handles.var2ShowLabel(:,2)';
handles.annVars = handles.var2ShowLabel(:,1)';

% views that depend on the current frame
handles.frameMenuShowLabels = {'Current Frame','Difference Image','Is Foreground','Connected Components'};

% all show labels
handles.menuShowLabels = ...
  [handles.annMenuShowLabels,...
  handles.frameMenuShowLabels];

% if in this list, then we scale the image from 0 to 255
handles.showRangeIs0to255 = {'Background Center','Background Median'};

% initialize to show background center
handles.showType = 'Background Center';

%% parse inputs
if isempty(varargin),
  error('Usage: BgModelReview(expdirs,[optional params]');
end
if ischar(varargin{1}),
  handles.expdirs = {varargin{1}}; %#ok<CCAT1>
else
  handles.expdirs = varargin{1};
end
[handles.rootdir,handles.ExpBGFGModelFile] = myparse(varargin(2:end),'rootdir','',...
  'ExpBGFGModelFile','');

%% initialize structures

% grab flags from GUI
handles.popupmenu_Flag_String = get(handles.popupmenu_Flag,'String');

% open ExpBGFGModel
if ~isempty(handles.ExpBGFGModelFile) && exist(handles.ExpBGFGModelFile,'file'),
  handles.ExpBGFGModel = readExpBGFGModel(handles.ExpBGFGModelFile);
end

handles.IsFirstExpDir = true;
handles = openExpDir(handles);

% initialize display
handles = InitGUI(handles);

% set to have data in it
handles = UpdateGUI(handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes BgModelReview wait for user response (see UIRESUME)
% uiwait(handles.figure1);

function handles = openExpDir(handles)

handles = closeExpDir(handles);

handles.expdir = fullfile(handles.rootdir,handles.expdirs{handles.expdiri});

% prepare to read movie
handles.moviefile = fullfile(handles.expdir,handles.moviefilestr);
if ~exist(handles.moviefile,'file'),
  error('Movie %s does not exist',handles.moviefile);
end
[handles.readframe,handles.nframes,handles.fid,handles.headerinfo] = get_readframe_fcn(handles.moviefile);
im = handles.readframe(1);
handles.nr = size(im,1);
handles.nc = size(im,2);
handles.ncolors = size(im,3);

% read annotation
if handles.DEBUG,
  tmp = regexp(handles.expdir,'/?([^/]+)/?$','tokens','once');
  handles.expdir_base = tmp{1};
  handles.annfile = fullfile(handles.resultsdir,[handles.expdir_base,'.ann']);
else
  handles.annfile = fullfile(handles.expdir,handles.annfilestr);
end
if ~exist(handles.annfile,'file');
  error('Annotation file %s does not exist',handles.annfile);
end
handles.ann = read_ann(handles.annfile);
if isempty(fieldnames(handles.ann)),
  warning('Could not parse annotation file %s',handles.annfile);
end

% resize images read from annotation
handles.annVarIsRead = false(1,size(handles.var2ShowLabel,1));
fns = handles.annVars;
for i = 1:length(fns),
  fn = fns{i};
  if isfield(handles.ann,fn),
    handles.ann.(fn) = reshape(handles.ann.(fn),[handles.nr,handles.nc,numel(handles.ann.(fn))/(handles.nr*handles.nc)]);
    handles.annVarIsRead(i) = true;
  end
end
if ~isfield('background_center',handles.ann) && isfield(handles.ann,'background_median'),
  handles.ann.background_center = handles.ann.background_median;
end
if ~isfield('background_dev',handles.ann) && isfield(handles.ann,'background_mad'),
  handles.ann.background_dev = min(max(handles.ann.background_mad,handles.ann.bg_std_min),handles.ann.bg_std_max);
end

% compute ExpBGFGModel log-lik of bg center
handles = UpdateExpBGFGModelLikelihood(handles);

% read trx
if handles.DEBUG,
  handles.trxfile = fullfile(handles.resultsdir,[handles.expdir_base,'.mat']);
else
  handles.trxfile = fullfile(handles.expdir,handles.trxfilestr);
end
if ~exist(handles.trxfile,'file');
  %error('Mat file %s does not exist',handles.trxfile);
end
%handles.trx = load_tracks(handles.trxfile);

% read the checking results
if handles.DEBUG,
  handles.savefile = fullfile(handles.resultsdir,[handles.expdir_base,'_CheckTrackingResults.txt']);
else
  handles.savefile = fullfile(handles.expdir,handles.savefilestr);
end
[handles.isOK,handles.isRedo,handles.flag,handles.notes,didread] = readCheckTrackingResults(handles.savefile);

set(handles.listbox_ExpDirs,'Value',handles.expdiri);
set(handles.popupmenu_Flag,'Value',find(strcmp(handles.popupmenu_Flag_String,handles.flag),1));
set(handles.edit_Notes,'String',handles.notes);

% current frame
handles.f = 1;

% make sure that the current view is available
i = find(strcmp(handles.showType,handles.annMenuShowLabels),1);
if isempty(i) || ~isfield(handles.ann,handles.annVars{i}),
  handles.showType = 'Current Frame';
end

handles = InitMenuShow(handles);

handles = InitSliderFrame(handles);

%handles = InitEditNotes(handles);

function handles = UpdateExpBGFGModelLikelihood(handles)

if ~isfield(handles.ann,'background_center') || ~isfield(handles,'ExpBGFGModel'),
  return;
end

[handles.ExpBGFGModel_Bg_LogLikRatio,...
  handles.ExpBGFGModel_Bg_LogLikGivenFg,...
  handles.ExpBGFGModel_Bg_LogLikGivenBg,...
  handles.ExpBGFGModel_Bg_LogLikAppearanceGivenFg,...
  handles.ExpBGFGModel_Bg_LogLikAppearanceGivenBg,...
  handles.ExpBGFGModel_Bg_LogLikDistObjGivenFg,...
  handles.ExpBGFGModel_Bg_LogLikDistObjGivenBg] = ...
  ExpBGFGModel_compute_log_lik_ratio(handles.background_center,handles.ExpBGFGModel);

function handles = closeExpDir(handles)

if handles.IsFirstExpDir,
  handles.IsFirstExpDir = false;
  return;
end

if isfield(handles,'fid') && handles.fid > 0,
  fclose(handles.fid);
end

%handles.flag{handles.expdiri} = handles.popupmenu_Flag_String{get(handles.popupmenu_Flag,'Value')};
%handles.notes{handles.expdiri} = get(handles.edit_Notes,'String');

function handles = InitGUI(handles)

hold(handles.axes_Video,'off');
if handles.ncolors == 1,
  handles.himage = imagesc(zeros(handles.nr,handles.nc),'Parent',handles.axes_Video,[0,255]);
  axis(handles.axes_Video,'image','off','xy');
else
  handles.himage = image(zeros([handles.nr,handles.nc,handles.ncolors],'uint8'),'Parent',handles.axes_Video);
end

% colorbar
axespos = get(handles.axes_Video,'Position');
listboxpos = get(handles.listbox_ExpDirs,'Position');
border = .01;
handles.colorbarPos = [axespos(1)+axespos(3)+border,axespos(2),listboxpos(1)-(axespos(1)+axespos(3))-border*2,axespos(4)];
handles.hColorbar = colorbar('peer',handles.axes_Video,'East');
set(handles.hColorbar,'Position',handles.colorbarPos,'XColor','w','YColor','w');

%handles = InitMenuShow(handles);

handles = InitListboxExpDirs(handles);

%handles = InitSliderFrame(handles);

function handles = InitListboxExpDirs(handles)

handles.listbox_ExpDirs_String = handles.expdirs;

set(handles.listbox_ExpDirs,'String',handles.listbox_ExpDirs_String,'Value',handles.expdiri);

function handles = InitMenuShow(handles)

chil = get(handles.menu_Show,'children');
if ~isempty(chil),
  delete(chil);
end
handles.menu_Show_Children = nan(1,length(handles.menuShowLabels));
for i = 1:length(handles.menuShowLabels),
  label = handles.menuShowLabels{i};
  j = find(strcmp(label,handles.annMenuShowLabels),1);
  if ~isempty(j) && ~handles.annVarIsRead(j),
    continue;
  end
  handles.menu_Show_Children(i) = uimenu(handles.menu_Show,'Label',...
    handles.menuShowLabels{i},'Checked','off',...
    'Callback',@(hObject,eventdata)BgModelReview('menu_Show_Children_Callback',hObject,eventdata,guidata(hObject)));
  if strcmp(handles.showType,label),
    set(handles.menu_Show_Children(i),'Checked','on');
  end
end

function handles = InitSliderFrame(handles)

% set slider steps
handles.sliderstep = [1/(handles.nframes-1),min(1,100/(handles.nframes-1))];
set(handles.slider_Frame,'Value',0,'SliderStep',handles.sliderstep);

%function handles = InitEditNotes(handles)
%
%set(handles.edit_Notes,'String',handles.notes{handles.expdiri});

function handles = UpdateGUI(handles)

switch handles.showType,
  case handles.annMenuShowLabels,
    i = find(strcmp(handles.showType,handles.annMenuShowLabels),1);
    annfn = handles.var2ShowLabel{i,1};
    if isfield(handles.ann,annfn),
      set(handles.himage,'CData',handles.ann.(annfn));
      if ismember(handles.showType,handles.showRangeIs0to255),
        set(handles.axes_Video,'CLim',[0,255]);
      else
        clim = double([min(handles.ann.(annfn)(:)),max(handles.ann.(annfn)(:))]);
        if clim(1) == clim(2),
          clim(1) = clim(1) - .5;
          clim(2) = clim(2) + .5;
        end
        set(handles.axes_Video,'CLim',clim);
      end
    else
      warning('BgModelReview:MissingAnnData','%s not read from annotation file',annfn);
    end
    colormap jet;
  case 'Current Frame',
    handles.im = handles.readframe(handles.f);
    set(handles.himage,'CData',handles.im);
    set(handles.axes_Video,'CLim',[0,255]);
    colormap jet;
  case 'Difference Image',
    handles = computeDiffIm(handles);
    set(handles.himage,'CData',handles.diffim);
    set(handles.axes_Video,'Clim',[0,max(handles.diffim(:))]);
    colormap jet;
  case 'Is Foreground',
    handles = computeIsFore(handles);
    set(handles.himage,'CData',handles.isfore);
    set(handles.axes_Video,'Clim',[0,1]);
    colormap gray;
  case 'Connected Components',
    handles = computeCC(handles);
    set(handles.himage,'CData',handles.ccim);
    set(handles.axes_Video,'Clim',[0,handles.cc.NumObjects]);
    colormap kjet;
end

set(handles.slider_Frame,'Value',(handles.f-1)/(handles.nframes-1));
set(handles.edit_Frame,'String',num2str(handles.f));

function handles = computeIsFore(handles)

handles = computeDiffIm(handles);

isfore0 = handles.diffim >= handles.ann.n_bg_std_thresh_low;
handles.isfore = false(size(isfore0));
handles.cc = bwconncomp(isfore0);
removecc = false(1,handles.cc.NumObjects);
if handles.ann.n_bg_std_thresh > handles.ann.n_bg_std_thresh_low;
  for i = 1:handles.cc.NumObjects,
    if max(handles.diffim(handles.cc.PixelIdxList{i})) > handles.ann.n_bg_std_thresh,
      handles.isfore(handles.cc.PixelIdxList{i}) = true;
    else
      removecc(i) = true;
    end
  end
end
handles.cc.NumObjects = nnz(~removecc);
handles.cc.PixelIdxList(removecc) = [];

function handles = computeCC(handles)

handles = computeIsFore(handles);
handles.ccim = zeros(handles.nr,handles.nc);
colors = randperm(handles.cc.NumObjects);
for i = 1:handles.cc.NumObjects,
  handles.ccim(handles.cc.PixelIdxList{i}) = colors(i);
end

function handles = computeDiffIm(handles)

handles.im = double(handles.readframe(handles.f));
if handles.ann.bg_type == 0,
  handles.diffim = max(0,handles.im - handles.ann.background_center);
elseif handles.ann.bg_type == 1,
  handles.diffim = max(0,handles.ann.background_center-handles.im);
else
  handles.diffim = max(handles.im - handles.ann.background_center,...
    handles.ann.background_center-handles.im);
end
handles.diffim = handles.diffim ./ handles.ann.background_dev;


% --- Outputs from this function are returned to the command line.
function varargout = BgModelReview_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --------------------------------------------------------------------
function menu_File_Callback(hObject, eventdata, handles)
% hObject    handle to menu_File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_File_Save_Callback(hObject, eventdata, handles)
% hObject    handle to menu_File_Save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.notes = get(handles.edit_Notes,'String');
writeCheckTrackingResults(handles.savefile,handles.isOK,handles.isRedo,handles.flag,handles.notes);

% --------------------------------------------------------------------
function menu_Show_Callback(hObject, eventdata, handles)
% hObject    handle to menu_Show (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function menu_Show_Children_Callback(hObject, eventdata, handles)
% hObject    handle to menu_Show (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.showType = get(hObject,'Label');
for i = 1:length(handles.menu_Show_Children),
  if ~ishandle(handles.menu_Show_Children(i)), continue; end
  if handles.menu_Show_Children(i) == hObject,
    set(handles.menu_Show_Children(i),'Checked','on');
  else
    set(handles.menu_Show_Children(i),'Checked','off');
  end
end
handles = UpdateGUI(handles);
guidata(hObject,handles);

% --- Executes when figure1 is resized.
function figure1_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on slider movement.
function slider_Frame_Callback(hObject, eventdata, handles)
% hObject    handle to slider_Frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
v = get(hObject,'Value');
handles.f = round(1 + v * (handles.nframes - 1));
if ismember(handles.showType,handles.frameMenuShowLabels),
  handles = UpdateGUI(handles);
end
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function slider_Frame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_Frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit_Frame_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Frame as text
%        str2double(get(hObject,'String')) returns contents of edit_Frame as a double
v = str2double(get(hObject,'String'));
if isnan(v),
  set(hObject,'String',num2str(handles.f));
  return;
end
handles.f = min(handles.nframes,max(1,v));
if ismember(handles.showType,handles.frameMenuShowLabels),
  handles = UpdateGUI(handles);
end
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_Frame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox_ExpDirs.
function listbox_ExpDirs_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_ExpDirs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_ExpDirs contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_ExpDirs
oldexpdiri = handles.expdiri;
handles.expdiri = get(hObject,'Value');
if oldexpdiri == handles.expdiri,
  return;
end
handles = openExpDir(handles);
handles = UpdateGUI(handles);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function listbox_ExpDirs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_ExpDirs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_OK.
function pushbutton_OK_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_OK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

s = ['<HTML><FONT color="green">',handles.expdirs{handles.expdiri},'</FONT></html>'];
handles.listbox_ExpDirs_String{handles.expdiri} = s;
set(handles.listbox_ExpDirs,'String',handles.listbox_ExpDirs_String);
handles.isOK = true;
handles.isRedo = false;
handles.notes = get(handles.edit_Notes,'String');
writeCheckTrackingResults(handles.savefile,handles.isOK,handles.isRedo,handles.flag,handles.notes);

if handles.expdiri < length(handles.expdirs),
  handles.expdiri = handles.expdiri+1;
end
handles = openExpDir(handles);
handles = UpdateGUI(handles);
guidata(hObject,handles);

guidata(hObject,handles);

% --- Executes on button press in pushbutton_Redo.
function pushbutton_Redo_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Redo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

s = ['<HTML><FONT color="red">',handles.expdirs{handles.expdiri},'</FONT></html>'];
handles.listbox_ExpDirs_String{handles.expdiri} = s;
set(handles.listbox_ExpDirs,'String',handles.listbox_ExpDirs_String);
handles.isOK = false;
handles.isRedo = true;
handles.notes = get(handles.edit_Notes,'String');
writeCheckTrackingResults(handles.savefile,handles.isOK,handles.isRedo,handles.flag,handles.notes);
if handles.expdiri < length(handles.expdirs),
  handles.expdiri = handles.expdiri+1;
end
handles = openExpDir(handles);
handles = UpdateGUI(handles);
guidata(hObject,handles);

guidata(hObject,handles);

guidata(hObject,handles);

function edit_Notes_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Notes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Notes as text
%        str2double(get(hObject,'String')) returns contents of edit_Notes as a double
handles.notes = get(hObject,'String');
writeCheckTrackingResults(handles.savefile,handles.isOK,handles.isRedo,handles.flag,handles.notes);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_Notes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Notes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_Flag.
function popupmenu_Flag_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_Flag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_Flag contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_Flag
handles.flag = handles.popupmenu_Flag_String{get(hObject,'Value')};
handles.notes = get(handles.edit_Notes,'String');
writeCheckTrackingResults(handles.savefile,handles.isOK,handles.isRedo,handles.flag,handles.notes);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function popupmenu_Flag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_Flag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'fid') && handles.fid > 0,
  fclose(handles.fid);
  handles.fid = -1;
end

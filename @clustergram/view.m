function view(obj, varargin)
%VIEW View a clustergram object in a MATLAB figure window.
%
%   VIEW(CGOBJ) shows a clustergram object. The hierarchical clustering is
%   shown with a heat map and dendrograms for the rows and columns.
%   Analysis tools for the dendrograms are accessible through the mouse
%   left/right buttons.
% 
%   Examples:
%       
%       view(cgobj)
%
%   See also CLUSTERGRAM, CLUSTERGRAM/PLOT.

%   Copyright 2007-2010 The MathWorks, Inc.
%   $Revision: 1.1.6.14 $  $Date: 2010/03/22 03:40:31 $

if isempty(obj)
    return;
end

if nargin == 2 && ishandle(obj.FigureHandle)
    appdata = getappdata(obj.FigureHandle, 'DendrogramData');
    if strcmpi(obj.ShowDendrogram, 'off')
        set(appdata.rowLines, 'Visible', 'off')
        set(appdata.colLines, 'Visible', 'off')
        disableDendrogramUIControls(obj.FigureHandle);
        set(obj.FigureHandle,'WindowButtonMotionFcn',[]);
    else
        set(appdata.rowLines, 'Visible', 'on')
        set(appdata.colLines, 'Visible', 'on')
        createDendrogramUIControls(obj);
        set(obj.FigureHandle, 'WindowButtonMotionFcn', {@localWindowButtonMotion, obj});
    end
    
    return;
elseif nargin == 4 && ishandle(obj.FigureHandle)
    % For clusterGroup calls
    groupIdx = varargin{1};
    dim = varargin{2};
    color = varargin{3};
    if isempty(color)
        highlightSelectGroup(obj, groupIdx, dim);
    else
        colorSelectGroup(obj, groupIdx, dim, color)
    end
    return;
end

if isempty(obj.FigureHandle) || ~ishandle(obj.FigureHandle)
    figureName = bioinfoprivate.indexedFigureName('Clustergram', 'Clustergram');
    obj.FigureHandle = figure('Renderer',' zbuffer',...
        'Name',figureName,...
        'NumberTitle','off',...
        'Tag', 'Clustergram',...
        'IntegerHandle','off',...
        'Handlevisibility', 'callback',...
        'Visible', 'on');

   updateToolbar(obj);
   updateUIMenus(obj);
else
    disableUIControls(obj.FigureHandle);
    rmappdata(obj.FigureHandle, 'DendrogramData');
    rmappdata(obj.FigureHandle, 'HeatMapListeners')
end

%== Plot the clustergram;
initDendrogramAppData(obj.FigureHandle)
delete(findall(obj.FigureHandle, 'Type', 'axes'))
plot(obj, obj.FigureHandle);

%%== Add marker dots to groups
createDendrogramUIControls(obj);
%== Add color marker Ui controls
createColorMarkerUIControls(obj)

%== Window Button Motion
set(obj.FigureHandle, 'WindowButtonMotionFcn', {@localWindowButtonMotion, obj});
set(obj.FigureHandle, 'WindowButtonUpFcn', {@localWindowButtonUp});

%== Handle data cursor on heatmap
updateFigureModes(obj, @heatmapDataCursorUpdatefcn)
end

function updateToolbar(obj)
% Helper function to update the toolbar
% Remove save, open edit plot buttons

oldSH = get(0,'ShowHiddenHandles');
set(0,'ShowHiddenHandles','on')

set(obj.FigureHandle,'toolbar','figure')  % needs to update because uicontrols turn it off

% Fix toolbar options, we keep: ZoomIn, ZoomOut, DataCursor, ColorBar,
% Print. Legend need more work.
hw = findall(obj.FigureHandle,'type','uitoolbar');
hf = get(hw,'Children');
h1 = findall(hf,'Tag','Exploration.ZoomOut');
h2 = findall(hf,'Tag','Exploration.ZoomIn');
% h3 = findall(hf,'Tag','Annotation.InsertColorbar');
% % h4 = findall(hf,'Tag','Annotation.InsertLegend');
h3 = findall(hf,'Tag','Exploration.DataCursor');
h4 = findall(hf,'Tag','Standard.PrintFigure');
delete(setxor(hf,[h1,h2,h3,h4]))

%InsertColorbar
iconfile = fullfile(toolboxdir('matlab'),'icons','tool_colorbar.gif');
hColorbar = bioinfoprivate.createToolbarButton(hw, 2, iconfile,...
    'Insert Colorbar', 'HMInsertColorbar', 'State', obj.Colorbar, 'separator', 'on');
set(hColorbar, 'ClickedCallback', {@insertColorbarCB, obj})
%Annote
iconfile = fullfile(toolboxdir('bioinfo'),'proteins','icons','reset_view.gif');
hAnnot = bioinfoprivate.createToolbarButton(hw, 2, iconfile,...
    'Annotate', 'HMAnnotateText', 'State', obj.Annotate);
set(hAnnot, 'ClickedCallback', {@annotateTextCB, obj})
%ShowDendromgram
iconfile = fullfile(toolboxdir('bioinfo'),'proteins','icons','hm_dendro.gif');
hShowDendro = bioinfoprivate.createToolbarButton(hw, 2, iconfile,...
    'Show Dendrogram', 'ShowDendrogramFlag', 'State', obj.ShowDendrogram);
set(hShowDendro, 'ClickedCallback', {@showDendrogramCB, obj})


set(0,'ShowHiddenHandles',oldSH)
end
%---------------------------------------------------------
function updateUIMenus(obj)
% Helper function to set UI menus
% Remove File->Save, Save as etc. menu items

oldSH = get(0,'ShowHiddenHandles');
set(0,'ShowHiddenHandles','on')

%== Delete figure menus not used
h1 = findall(obj.FigureHandle,'Type','uimenu', 'Label','&Edit');
h2 = findall(obj.FigureHandle,'Type','uimenu', 'Label','&View');
h3 = findall(obj.FigureHandle,'Type','uimenu', 'Label','&Insert');
delete([h1,h2, h3])

%== Repair "File" menu
hw = findall(obj.FigureHandle,'Type','uimenu', 'Label','&File');
hf = get(hw,'children');
h1 = findall(hw,'Label','Expo&rt Setup...');
h3 = findall(hw,'Label','Print Pre&view...');
h4 = findall(hw,'Label','&Print...');
delete(setxor(hf,[h1,h3,h4]))

uimenu(hw ,'Label','Print to Figure',...
                   'Position', 1,...
                   'Callback',{@dotPrintToFigureCB, obj, false});
uimenu(hw,'Label', 'Copy to New Clustergram',...
                   'Position',2,...
                   'Callback',{@dotExportCB, obj, false, false});
uimenu(hw,'Label', 'Export to Workspace',...
                   'Position',3,...
                   'Callback',{@dotExportCB, obj, true, false});
uimenu(hw,'Label', 'Exit','Separator','on',...
                   'Position',7,'Callback',...
                   {@closeClustergramFigure, obj})
set(h1,'Separator','on')
    
%==Repair "Tools" menu
hw = findall(obj.FigureHandle,'Type','uimenu','Label','&Tools');
hf = get(hw,'children');
h1 = findall(hw,'Tag','figMenuZoomIn');   
h2 = findall(hw,'Tag','figMenuZoomOut'); 
% % h3 = findall(hw,'Tag','figMenuPan');     
h3 = findall(hw,'Tag','figMenuResetView'); 
h4 = findall(hw,'Tag','figMenuOptions');
h5 = findall(hw,'Tag','figMenuDatatip');
set([h1,h3], 'separator','off')
delete(setxor(hf,[h1,h2,h3,h4,h5]))

% Add Colorbar switch
h6 = bioinfoprivate.createMenuItem(obj.FigureHandle, hw, 'Insert Colorbar',...
                                   'HMInsertColorbar',...
                                   'Position', 6,...
                                   'Checked', obj.Colorbar);
set(h6,'Separator','on')
% Add Annotation switch
h7 = bioinfoprivate.createMenuItem(obj.FigureHandle, hw, 'Annotate',...
                                   'HMAnnotateText',...
                                   'Position', 7,...
                                   'Checked', obj.Annotate);
% Add Dendrogram switch
h8 = bioinfoprivate.createMenuItem(obj.FigureHandle, hw, 'Show Dendrogram',...
                                   'ShowDendrogramFlag',...
                                   'Position', 8,...
                                   'Checked', obj.ShowDendrogram);                               
setappdata(obj.FigureHandle, 'ClustergramMenuItems', [h6, h7, h8]);

% Repair "Help" menu
bioinfoprivate.bioFigureHelpMenu(obj.FigureHandle, 'Clustergram', 'clustergram_refpage')

set(0,'ShowHiddenHandles',oldSH)
end
%------------------

function createDendrogramUIControls(obj)
% % function createDendrogramUIControls(hFig, obj)
% Create other UI components for dendrograms
appdata = getappdata(obj.FigureHandle, 'DendrogramData');

%== Add group marker dots
if ~isempty(appdata.rowLines)
    [appdata.hRowGroupDot,...
     appdata.xRowGroupDot,...
     appdata.yRowGroupDot]= getGroupDots(appdata.rowLines, 1);
end

if ~isempty(appdata.colLines)
    [appdata.hColGroupDot,...
     appdata.xColGroupDot,...
     appdata.yColGroupDot] = getGroupDots(appdata.colLines, 2);
end

% Add context menu to the dots
appdata.dotContextMenu = uicontextmenu('Parent',   obj.FigureHandle,...
                                       'Tag',      'GroupDotContextmenu',...
                                       'Callback', @hideDatatipCB);
uimenu(appdata.dotContextMenu, 'Label',    'Set Group Color',...
                               'Callback', {@dotChangeColorCB, obj});
uimenu(appdata.dotContextMenu, 'Label',     'Print Group to Figure',...
                               'Separator', 'on',...
                               'Callback',  {@dotPrintToFigureCB, obj, true});
uimenu(appdata.dotContextMenu, 'Label',    'Copy Group to New Clustergram',...
                               'Callback', {@dotExportCB, obj, false, true});
uimenu(appdata.dotContextMenu, 'Label',     'Export Group to Workspace',...
                               'Callback', {@dotExportCB, obj, true, true});
uimenu(appdata.dotContextMenu, 'Label',    'Export Group Info to Workspace',...
                               'Callback', {@dotExportGroupInfoCB, obj});
uimenu(appdata.dotContextMenu, 'Label',    'Show details',...
                               'Callback', {@dotShowGroupDetails, obj});

%== Add Buttondown function to dots
if ~isempty(appdata.hColGroupDot)
    set(appdata.hColGroupDot,'ButtonDownFcn',{@dotButtonDownCB, obj},...
                             'UIContextMenu', appdata.dotContextMenu)
end

if ~isempty(appdata.hRowGroupDot)
    set(appdata.hRowGroupDot,'ButtonDownFcn',{@dotButtonDownCB, obj},...
                              'UIContextMenu', appdata.dotContextMenu)
end

%== Create dendrogram datatip
appdata.datatip = getDatatipText('DendrogramDatatip', obj.HMAxesHandle);
appdata.highlight = line(1,1, 'Parent',    obj.HMAxesHandle,...    
                              'Color',     appdata.highlightColor,...
                              'LineWidth', 1.5,...
                              'Visible',   'off');
bioma.util.disableHGBehaviors(appdata.highlight)
setappdata(obj.FigureHandle, 'DendrogramData', appdata);
end

function closeClustergramFigure(h, evt, obj) %#ok
    close(obj.FigureHandle)
end
%--------------------------

function createColorMarkerUIControls(obj)
% Create other UI components for color markers
appdata = getappdata(obj.FigureHandle, 'DendrogramData');

% There are not markers
if isempty( appdata.colMarkers) &&  isempty(appdata.rowMarkers)
    return;
end

% Context menu for the markers
appdata.markerContextMenu = uicontextmenu('Parent',    obj.FigureHandle,...
                                           'Tag',      'ColorMarkerContextmenu',...
                                           'Callback', @hideDatatipCB);
uimenu(appdata.markerContextMenu, 'Label',     'Change Color',...
                                  'Callback', {@cmarkerChangeColorCB, obj});
uimenu(appdata.markerContextMenu, 'Label',    'Update Annotation',...
                                  'Callback', {@cmarkerUpdateAnnoCB, obj});

%== Add Buttondown function and context menu to color markers
if ~isempty(appdata.colMarkers)
    set(appdata.colMarkers, 'ButtonDownFcn', {@cmarkerButtonDownCB},...
                            'UIContextMenu', appdata.markerContextMenu)
end

if ~isempty(appdata.rowMarkers)
    set(appdata.rowMarkers, 'ButtonDownFcn', {@cmarkerButtonDownCB},...
                            'UIContextMenu', appdata.markerContextMenu)
end

%== Create dendrogram datatip
appdata.markerDatatip = getDatatipText('ColorMarkerDatatip', obj.HMAxesHandle);
setappdata(obj.FigureHandle, 'DendrogramData', appdata);
end
%------------------
function [hdots, xdots, ydots] = getGroupDots(lineH, dir)
% Return handles to the dots for each dendrogram group. The dot is to be
% drawn at the center of the group lines -lineH
% Dir - direction of the axes, 1-along column, 2-along row.

nGroup = numel(lineH); 
hdots = zeros(nGroup,1);
xdots = zeros(nGroup,1);
ydots = zeros(nGroup,1);
axes = get(lineH(1), 'Parent');

for i=1:nGroup
    if dir == 1
        xdots(i) = max(get(lineH(i), 'XData'));
        ydots(i) =(max(get(lineH(i), 'YData')) + min(get(lineH(i), 'YData')))/2;
    elseif dir == 2
        xdots(i) =(max(get(lineH(i), 'XData')) + min(get(lineH(i), 'XData')))/2;
        ydots(i) = max(get(lineH(i), 'YData'));
    end
    
    hdots(i) = line(xdots(i), ydots(i), 'Parent',         axes,...
                                        'Marker',         'o',...
                                        'MarkerSize',     5,...
                                        'MarkerEdgeColor','r',...
                                        'MarkerFaceColor','r',...
                                        'Tag',             num2str(i),...
                                        'Visible', 'off');
    bioma.util.disableHGBehaviors(hdots(i));
end
end

%== callbacks
function hideDatatipCB(~,varargin)
appdata = getappdata(gcbf, 'DendrogramData');
set(appdata.datatip, 'Visible', 'off')
end
%---------------------------
function cmarkerChangeColorCB(~, ~, obj)
% Change color marker color
rect = gco;
if rect
    c = uisetcolor;
    if c == 0
        return;
    end
    
    haxes = get(rect, 'Parent');    
    idx = str2double(get(rect, 'Tag'));
    switch get(haxes, 'Tag')
        case 'CMarkerRowAxes'
            obj.RowGroupMarker(idx).Color = c;
        case 'CMarkerColAxes'
            obj.ColumnGroupMarker(idx).Color = c;
    end
    set(rect, 'FaceColor', c);
end
end

%---------------------------
function cmarkerUpdateAnnoCB(~, ~, obj)
% Change color marker color
% % appdata = getappdata(gcbf, 'DendrogramData');
appdata = getappdata(obj.FigureHandle, 'DendrogramData');

rect = gco;
if rect
    s = inputdlg('Enter Annotation','Update Color Marker Annotation',1);
    if isempty(s)
        return;
    end
    s = s{:};
    haxes = get(rect, 'Parent');
    dim = 1;
    if strcmpi(get(haxes, 'Tag'), 'CMarkerColAxes')
        dim = 2;
    end

   idx = str2double(get(rect, 'Tag'));
    if dim == 1
        %-- Update the text string
        old_s = get(rect, 'DisplayName');
        tidx = strcmpi(get(appdata.rowMarkerAnnotation, 'String'), old_s);
        set(appdata.rowMarkerAnnotation(tidx), 'String', s);
        
        obj.RowGroupMarker(idx).Annotation = s;
    elseif dim == 2
        obj.ColumnGroupMarker(idx).Annotation = s;
    end
    set(rect, 'DisplayName', s);
    
end
end
%------------------
function dotButtonDownCB(h, ~, obj)
% Callback to mouse button click on a branch dot. 
% When clicked it should the group data tip and highlight the nodes under
% the group
hfig = gcbf;

groupidx = str2double(get(h, 'Tag'));
dim = 1;
if strcmpi(get(get(h, 'Parent'), 'Tag'), 'DendroColAxes')
    dim = 2;
end

% Show datatip 
showDatatip(hfig, obj, groupidx, dim)
end
%---------------------------
function dotChangeColorCB(~, ~, obj)
% Change color group color
appdata = getappdata(obj.FigureHandle, 'DendrogramData');
if appdata.selectDot
    c = uisetcolor;
    if c==0
        return;
    end
    updateGroupColors(obj, appdata.selectDot, c, 2);
end
end
%---------------------------------------------
function dotPrintToFigureCB(~, ~, obj, groupflag)
% Print selected group to a new Figure
hFig = gcbf;
if groupflag
    newcg_o = getGroupObject(obj, hFig);
else
    newcg_o = createNewClustergramObj(obj);
    %= Also need to get group selection colors
    appdata = getappdata(hFig, 'DendrogramData');
    newcg_o.DendroRowLineColor = appdata.rowLineColor;
    newcg_o.DendroColLineColor = appdata.colLineColor;
end
plot(newcg_o)
end
%-------------------------------------------------------
function newcg_o = getGroupObject(obj, hfig)
% Return a new clustergram object from the selected groups.
appdata = getappdata(hfig, 'DendrogramData');
newcg_o = [];
if appdata.selectDot
    groupIdx = str2double(get(appdata.selectDot, 'Tag'));
    
    switch get(get(appdata.selectDot, 'Parent'), 'Tag')
        case 'DendroRowAxes'
             fDim = 1;
        case 'DendroColAxes'
            fDim = 2;
    end
    
    newcg_o = getDendroGroupObject(obj, groupIdx, fDim);
end
end

%-------------------------------------------------------
function dotExportCB(~, ~, obj, towsflag, groupflag)
% Export selected groups or clustergram to workspace

hFig = gcbf;
if groupflag
    newcg_o = getGroupObject(obj, hFig);
else
    newcg_o = createNewClustergramObj(obj);
    %= Also need to get group selection colors
    appdata = getappdata(hFig, 'DendrogramData');
    newcg_o.DendroRowLineColor = appdata.rowLineColor;
    newcg_o.DendroColLineColor = appdata.colLineColor;
    
    if ~obj.CopyOnly
        newcg_o.CopyOnly = false;
    end
end
if towsflag % export to workspace
    s = inputdlg('Workspace variable name ?','Export to Workspace',1);
    while ~(isempty(s) || isvarname(s{1}) || isempty(s{1}))
        s = inputdlg('Not a valid variable name, type a MATLAB variable name ?',...
            'Export to Workspace',1);
    end
    if ~(isempty(s) || isempty(s{1}))
        assignin('base',s{1}, newcg_o)
    end
else % no, then export to other viewer
    view(newcg_o);
end

end


%-------------------------------------------------------
function dotExportGroupInfoCB(~, ~, obj)
% Export selected group information to workspace
% Information structure
appdata = getappdata(gcbf, 'DendrogramData');
s = inputdlg('Workspace variable name ?','Export Group Information to Workspace',1);
while ~(isempty(s) || isvarname(s{1}) || isempty(s{1}))
    s = inputdlg('Not a valid variable name, type a MATLAB variable name ?',...
        'Export to Workspace',1);
end

if ~(isempty(s) || isempty(s{1}))    
    dim = 1;
    if strcmpi(get(get(appdata.selectDot, 'Parent'), 'Tag'), 'DendroColAxes')
        dim = 2;
    end
    
    groupIdx = str2double(get(appdata.selectDot, 'Tag'));
   
    istruct = getGroupInfo(obj, groupIdx, dim);
    assignin('base',s{1}, istruct)
end
end

function dotShowGroupDetails(h,~,obj)

persistent metadata;

if isempty(metadata),
  [f,p] = uigetfile('*.mat','Select metadata mat file');
  if ~ischar(f),
    return;
  end
  metadatafile = fullfile(p,f);
  metadata = load(metadatafile,'metadata');
  metadata = metadata.metadata;
end
  

hfigout = 12345;
hfigout2 = 12346;
anatomydir = '/nobackup/branson/AverageAnatomyData20130618';

hFig = gcbf;
appdata = getappdata(hFig, 'DendrogramData');

% if this is a column
if strcmpi(get(get(appdata.selectDot, 'Parent'), 'Tag'), 'DendroRowAxes'),
  fprintf('Not implemented.\n');
  return;
end

newcg_o = getGroupObject(obj, hFig);
short_line_names = get(newcg_o,'ColumnLabels');

isad = ~cellfun(@isempty,regexp(short_line_names,'D$','once'));
line_names = cell(size(short_line_names));
if any(~isad),
  line_names(~isad) = cellfun(@(x) ['GMR_',x(2:end),'_AE_01'],short_line_names(~isad),'UniformOutput',false);
end
if any(isad),
  line_names(isad) = cellfun(@(x) ['GMR_',x(2:end-1),'_AD_01'],short_line_names(isad),'UniformOutput',false);
end
nlines = numel(line_names);

% create an html page with links to all experiments
timestamp = datestr(now,'yyyymmddTHHMMSS');
filenamecurr = fullfile(tempdir,sprintf('selectedgroup_%s.html',timestamp));
fid = fopen(filenamecurr,'w');

fprintf(fid,'<html>\n<title>Group selected %s</title>\n<body>\n',timestamp);
for i = 1:nlines,
  fprintf(fid,'<h1>%s</h1>\n',line_names{i});
  idx = find(strcmp({metadata.line_name},line_names{i}));
  if isempty(idx),
    fprintf(fid,'<p>No experiments found.</p>\n');
  else
    fprintf(fid,'<ul>\n');
    for j = idx(:)',
      [~,name] = fileparts(metadata(j).file_system_path);
      filename = fullfile(metadata(j).file_system_path,...
        sprintf('ctrax_results_movie_%s.avi',name));
      if ~exist(filename,'file'),
        continue;
      end
      fprintf(fid,'  <li><a href="file://%s">%s</a></li>\n',filename,name);
    end
    fprintf(fid,'</ul>\n');
  end
end
fprintf(fid,'</body>\n</html>\n');
fclose(fid);

if ~exist(filenamecurr,'file'),
  warning('Could not open temporary file %s',filenamecurr);
  return;
end
% open this page
web(filenamecurr,'-browser');

figure(hfigout);
clf;
nc = ceil(sqrt(nlines));
nr = ceil((nlines)/nc);
hax = createsubplots(nr,nc,[.01,.01],hfigout);
if numel(hax) > nlines,
  delete(hax(nlines+1:end));
end

mu = single(0);
nread = 0;
hwait = mywaitbar(0,'Reading images');
for i = 1:nlines,
  filename = fullfile(anatomydir,sprintf('meanim_%s.png',line_names{i}));
  hwait = mywaitbar((i-1)/(nlines+1),hwait,sprintf('Reading %s...\n',line_names{i}));
  if exist(filename,'file'),
    im = imread(filename);
    image(im,'Parent',hax(i));
    axis(hax(i),'image');
  end
  title(hax(i),line_names{i},'Interpreter','none');
  set(hax(i),'XTick',[],'YTick',[]);
  drawnow;
  filename = fullfile(anatomydir,sprintf('meanim_%s.mat',line_names{i}));
  if exist(filename,'file'),
    tmp = load(filename,'meanim');
    mu = mu + tmp.meanim;
    nread = nread+1;
  end
end

mu = mu / nread;
figure(hfigout2);
clf;
hax(nlines+1) = gca;
imagesc(max(mu,[],3)','Parent',hax(nlines+1));
axis(hax(nlines+1),'image');
title(hax(nlines+1),'Mean image','Interpreter','none');
colormap(hfigout,kjetsmooth(256));
set(hax(nlines+1),'XTick',[],'YTick',[]);
colorbar('peer',hax(nlines+1));
impixelinfo(hfigout);
impixelinfo(hfigout2);
linkaxes(hax(1:nlines+1));

if ishandle(hwait),
  delete(hwait);
end


end

%---------------------------------
function cmarkerButtonDownCB(src,~)
% src - the object that is the source of the event
% evnt - empty for this property
appdata = getappdata(gcbf, 'DendrogramData');
rect = src;
haxes = get(rect, 'Parent');
dim = 1;
if strcmpi(get(haxes, 'Tag'), 'CMarkerColAxes')
    dim = 2;
end
if rect    
    set(appdata.markerDatatip, 'Parent', haxes,...
                                'String', get(rect, 'DisplayName'))
                            
    %== Adjust for position
    cphm = get(haxes,'CurrentPoint'); % Data unit
    xl = get(haxes,'Xlim'); % Data unit
    yl = get(haxes,'Ylim');
    xDelta=diff(xl)/100;
    yDelta=diff(yl)/100;
    halign = 'Left';
    valign = 'Top';
    ext = get(appdata.markerDatatip, 'Extent');

    if dim == 1 % Row Cmarkers
        topEdge = cphm(1,2) + yDelta;
        leftEdge = cphm(1,1)-1+ 3*xDelta;
    else % Col CMarkers
        topEdge = cphm(1,2) + yDelta;
        if ext(3)+cphm(1,1) > xl(2)
            leftEdge = cphm(1,1)- 3*xDelta;
            halign = 'Right';
        else
            leftEdge = cphm(1,1)+ 3*xDelta;
        end
        valign = 'Bottom';
    end

    set(appdata.markerDatatip, 'Position', [leftEdge, topEdge, 1],...
        'VerticalAlignment', valign,...
        'Horizontalalignment',halign,...
        'Visible','on')
end
end
%-----------------------------------------------------------
function localWindowButtonUp(~, ~)
% Window mouse released 
hFig = gcbf;
appdata = getappdata(hFig, 'DendrogramData');
if ~isempty(appdata.markerDatatip)
    set(appdata.markerDatatip, 'Visible', 'off');
end
set(appdata.datatip, 'Visible', 'off');
setappdata(hFig, 'DendrogramData', appdata);
end

%----------------------------------------------
function localWindowButtonMotion(~, ~, obj)
% Callback function activated when moving over the axes, checks location of
% the mouse and puts datatip if over an active node.

% set a virtual grid to get the point
hFig = gcbf;
appdata = getappdata(hFig, 'DendrogramData');
set(appdata.selectDot, 'visible', 'off')
set(appdata.highlight, 'visible', 'off')

% Get the drendrogram axes the mouse is on
onAxes = [];
groupDots = [];
updateGroupColors(obj, appdata.selectDot, [], 3)

if ~isempty(appdata.rowDendroAxes)
    xl = get(appdata.rowDendroAxes,'Xlim');
    yl = get(appdata.rowDendroAxes,'Ylim');
    cp = get(appdata.rowDendroAxes,'CurrentPoint');
    xPos = cp(1,1); yPos = cp(1,2);
    
    if xPos >= xl(1) && xPos <= xl(2) &&  yPos >= yl(1) && yPos <= yl(2)
        onAxes = appdata.rowDendroAxes; 
        xdendro = appdata.xRowGroupDot;
        ydendro = appdata.yRowGroupDot;
        groupDots = appdata.hRowGroupDot;
    end
    set(groupDots, 'visible', 'off')
end

if isempty(onAxes) && ~isempty(appdata.colDendroAxes)
    xl = get(appdata.colDendroAxes,'Xlim');
    yl = get(appdata.colDendroAxes,'Ylim');
    cp = get(appdata.colDendroAxes,'CurrentPoint');
    xPos = cp(1,1); yPos = cp(1,2);
    
    if xPos >= xl(1) && xPos <= xl(2) &&  yPos >= yl(1) && yPos <= yl(2)
        onAxes = appdata.colDendroAxes;
        xdendro = appdata.xColGroupDot;
        ydendro = appdata.yColGroupDot;
        groupDots = appdata.hColGroupDot;
    end
    set(groupDots, 'visible', 'off')
end

if ~isempty(onAxes)
    xThres=diff(xl)/70; %100;
    yThres=diff(yl)/70; %100;

    hp = xdendro < (xPos+xThres) & xdendro > (xPos-xThres) & ...
         ydendro < (yPos+yThres) & ydendro > (yPos-yThres);
    hp = find (hp);

    if isempty(hp)
        set(appdata.datatip, 'visible','off')
        set(groupDots, 'visible', 'off')
    elseif numel(hp)
% %         %== Turn on the dot for the selected group
% %         % Turn off activated modes 
% %         deactivateMode(hFig);

        %== If any UITools mode is active, do nothing.
        hMode = deactivateMode(hFig, false);
        if isempty(hMode)
            set(groupDots(hp(1)), 'visible', 'on')
            appdata.selectDot = groupDots(hp(1));
            updateGroupColors(obj, appdata.selectDot, appdata.highlightColor, 1);
        end
    end
    setappdata(hFig, 'DendrogramData', appdata);
end
end

%------------
function highlightSelectGroup(obj, groupIdx, dim)
%== Turn on the dot for the selected group with specified groupIdx

% Turn off activated modes
appdata = getappdata(obj.FigureHandle, 'DendrogramData');
switch dim
    case 1
        groupDots = appdata.hRowGroupDot;
    case 2
        groupDots = appdata.hColGroupDot;
end
%== Don't deactivate any UITools mode.
% % deactivateMode(obj.FigureHandle);
set(groupDots(groupIdx), 'visible', 'on')
appdata.selectDot = groupDots(groupIdx);
updateGroupColors(obj, appdata.selectDot, appdata.highlightColor, 1);
setappdata(obj.FigureHandle, 'DendrogramData', appdata);
end

function hMode = deactivateMode(hFig, deactiveFlag) 
    %==Turn off all activated modes except datatcursor if deactive. Return
    %the handle to current UITools.mode: zoom or datacursor etc.
    hManager = uigetmodemanager(hFig);
    hMode = get(hManager,'CurrentMode');
    
    if deactiveFlag
        if strcmpi(get(hMode, 'Name'), 'Exploration.Datacursor')
            return;
        end
        
        %== Set turn off current mode
        set(hManager,'CurrentMode',[]);
    end
end
%== Helper functions---------------------------------------------

function xedge = getGroupEdge(xd)
% Helper function for find the edge of a group
xedge = zeros(1, 2);
if size(xd, 1) > 1
    xd=cell2mat(xd);
    xedge(1) = fix(min(min(xd)));
    xedge(2) = fix(max(max(xd)));
else
    xedge(1) = fix(min(xd));
    xedge(2) = fix(max(xd));
end
end
%--------------
function updateGroupColors(obj, hdot, grpColor, type)
% Set the colors of all under selected group
% Type: highlight, update and recover
hFig = obj.FigureHandle;
appdata = getappdata(hFig, 'DendrogramData');

if hdot
    idx = str2double(get(hdot, 'Tag'));
    
    switch get(get(hdot, 'Parent'), 'Tag')
        case 'DendroRowAxes'
            dim = 1;
            dendroLH = appdata.rowLines;
            tickText = getappdata(obj.HMAxesHandle, 'YTickLabelTextHandles');
        case 'DendroColAxes'
            dim = 2;
            dendroLH = appdata.colLines;
            tickText = getappdata(obj.HMAxesHandle, 'XTickLabelTextHandles');
    end

    % Propagate to get all the children under the group
    sel_idx = clusterPropagation(obj, idx, dim);
    
    switch type
        case 1 % highlight
            % Update color
            set(dendroLH(sel_idx), 'Color', grpColor)
            if ~isempty(tickText)
                idx = getTickTextIndex(dendroLH(sel_idx), numel(tickText),...
                                       obj.HMAxesHandle, dim);
                set(tickText(idx), 'Color', grpColor)
            end
            % Set highligh box
            xd = get(dendroLH(sel_idx), 'XData');
            yd = get(dendroLH(sel_idx), 'YData');
            if dim == 1 %row
                xl = get(obj.HMAxesHandle, 'Xlim');
                yl = getGroupEdge(yd) + [-0.5 0.5];
                hl_xd = [xl(1) xl(2) xl(2) xl(1) xl(1)];
                hl_yd = [yl(1) yl(1) yl(2) yl(2) yl(1)];
            else % column
                xl = getGroupEdge(xd)+ [-0.5 0.5];
                yl = get(obj.HMAxesHandle, 'Ylim');               
                hl_xd = [xl(1) xl(2) xl(2) xl(1) xl(1)];
                hl_yd = [yl(2) yl(2) yl(1) yl(1) yl(2)];
                
            end
            set(appdata.highlight, 'XData', hl_xd,...
                                   'YData', hl_yd,...
                                   'Visible', 'on')
        case 2 % update
            set(dendroLH(sel_idx), 'Color', grpColor)
            if dim == 1
                appdata.rowLineColor(sel_idx,:) = repmat(grpColor, numel(sel_idx), 1);
            elseif dim == 2
                appdata.colLineColor(sel_idx,:) = repmat(grpColor, numel(sel_idx), 1);
            end
             setappdata(hFig, 'DendrogramData', appdata);
        case 3 % recover
            if dim == 1
                cd = appdata.rowLineColor;
            else
                cd = appdata.colLineColor;
            end
            for i = 1:size(cd, 1)
                set(dendroLH(i), 'Color', cd(i, :))
            end
            if ~isempty(tickText)
                set(tickText, 'Color', 'k')
            end
    end
end
end

function colorSelectGroup(obj, groupIdx, dim, grpColor)
% Update dendrogram line color to specified group color
appdata = getappdata(obj.FigureHandle, 'DendrogramData');
switch dim
    case 1
        hDendroLines = appdata.rowLines;
    case 2
        hDendroLines = appdata.colLines;
end

% Propagate to get all the children under the group
selIdx = clusterPropagation(obj, groupIdx, dim);
set(hDendroLines(selIdx), 'Color', grpColor)

if dim == 1
    appdata.rowLineColor(selIdx,:) = repmat(grpColor, numel(selIdx), 1);
elseif dim == 2
    appdata.colLineColor(selIdx,:) = repmat(grpColor, numel(selIdx), 1);
end
setappdata(obj.FigureHandle, 'DendrogramData', appdata);
end


function idx = getTickTextIndex(selLines, numTixkText, hAxes, dim)
switch dim
    case 1
        xd = get(selLines, 'Ydata');
        alim = ceil(get(hAxes, 'Ylim'));
    case 2
        xd = get(selLines,     'Xdata');
        alim = ceil(get(hAxes, 'Xlim'));
end
xl = getGroupEdge(xd);
idx = find(ismember(alim(1):alim(2), (xl(1):xl(2))));
idx(idx > numTixkText) = [];
end

%--------------------------------
function showDatatip(hFig, obj, groupidx, dim)
% Show dendrogram datatips
appdata = getappdata(hFig, 'DendrogramData');

switch dim
    case 1
        labels = obj.RowLabels;
        objGroupNames = appdata.rowGroupNames;
    case 2
        labels = obj.ColumnLabels;
        objGroupNames = appdata.colGroupNames;
end

[~, sel_nodes] = clusterPropagation(obj, groupidx, dim);

%== Place text
name = objGroupNames{groupidx};
numChil = numel(sel_nodes);
if dim ==1
    childrenNames = char(flipud(labels(sel_nodes))); 
else
    childrenNames = char(labels(sel_nodes));
end
childrenNames=[repmat('   ',size(childrenNames,1),1) childrenNames];

if numChil
    set(appdata.datatip,'String',char([
        {[name '  (' num2str(numChil) ' nodes)']};...
        mat2cell(childrenNames,ones(size(childrenNames,1),1),...
        size(childrenNames,2))]))
else
    set(appdata.datatip,'string', name)
end
extraLines = 1;

%compute some values before adjusting data tip
% Use point units
figunit_orig = get(hFig, 'unit');
set(hFig, 'Units', 'points')

fp = get(hFig,'Position'); % fig position in points
fh = fp(4);      % fig size (height & width) in points
reqPt  = (numChil+extraLines)*14+2;  % required datatip height in pts

%Adjust string of datatip if it will not fit
if reqPt > fh
    str = get(appdata.datatip,'String');
    set(appdata.datatip,'String',str(1:extraLines,:));
    reqPt  = (extraLines)*14+2;
end

cphm = get(obj.HMAxesHandle,'CurrentPoint'); % Data unit
xl = get(obj.HMAxesHandle,'Xlim'); % Data unit
yl = get(obj.HMAxesHandle,'Ylim');
xDelta=diff(xl)/100;
yDelta=diff(yl)/100;
halign = 'Left';
ext = get(appdata.datatip, 'Extent');
if dim == 1
    yPosPt = ((cphm(1,2)+3*yDelta))*fh/diff(yl); % Datatip top edge in points
    % datatip position in data unit
    if reqPt <= yPosPt
        topEdge = cphm(1,2)+3*yDelta;
    else
        topEdge = cphm(1,2)+3*yDelta+ diff([yPosPt, reqPt]) * (diff(yl)/fh);
    end
    
    del = topEdge - ext(4);
    if del < yl(1)
       topEdge = topEdge + (yl(1)-del)+3*yDelta;
    end
    
    leftEdge = cphm(1,1) + 3*xDelta;
else
    topEdge = cphm(1,2) + yDelta;
    if ext(3)+cphm(1,1) > xl(2)
        leftEdge = cphm(1,1)- 3*xDelta;
        halign = 'Right';
    else
        leftEdge = cphm(1,1)+ 3*xDelta;
    end
end

set(appdata.datatip, 'Position',[leftEdge, topEdge, 1],...
                     'Horizontalalignment',halign,...
                     'Visible','on')

set(hFig, 'Units', figunit_orig)

setappdata(hFig, 'DendrogramData', appdata);
end
%---------------------------------------
%== Data cursor for heatmap
function datacursorLabel = heatmapDataCursorUpdatefcn(~, event_obj, dc_obj, obj)
% Display data cursor in heatmap image.
datacursorLabel = [];
if strcmpi(get(dc_obj, 'SnapToDataVertex'), 'off')
    return;
end

tg = get(event_obj, 'Target');
if strcmpi(get(tg, 'Type'), 'image')
    type = get(tg, 'Tag');
    if strcmpi(type, 'HeatmapImage')        
        pos = get(event_obj, 'Position');
        cdata = get(tg, 'CData');
        val = cdata(pos(2), pos(1));
        
        origval = obj.OriginalData(pos(2), pos(1));
        scale = getScale(obj, pos(2), pos(1));

        if isempty(scale)
            valstr=sprintf('%0.2f(Value:%0.2f)', val, origval);
        else
            valstr=sprintf('%0.2f(Value:%0.2f,Std:%0.2f)', val, origval,scale);
        end

        datacursorLabel = {valstr,...
                           obj.RowLabels{pos(2)},...
                           obj.ColumnLabels{pos(1)}};
    else
        datacursorLabel = [];
        return;
    end
else
    datacursorLabel = [];
    return;
end
end

function scale = getScale(obj, row, col)
scale = obj.Scales;
switch obj.Standardize
    case 'COLUMN'
        scale = scale(obj.DendroColPerm);
        if ~isempty(obj.ColNodes)
            scale = scale(obj.ColNodes);
        end
        scale = scale(col);
    case 'ROW'
        scale = scale(obj.DendroRowPerm);
        if ~isempty(obj.RowNodes)
            scale = scale(obj.RowNodes);
        end
        scale = scale(row);
end
end
%== Axes helper functions
%------------------------------------------------
function datatip = getDatatipText(tag, hax)
% Return a text object for datatip with specified tag.
datatip = text(0,1,1,'k', 'Parent', hax,...
                          'Tag', tag,...
                          'BackgroundColor',[1 1 .93],...
                          'Color', [0 0 0],...
                          'EdgeColor', [0.8 0.8 0.8],...
                          'VerticalAlignment','Top',...
                          'Clipping','off',...
                          'Visible','off',...
                          'Fontsize',8,...
                          'Interpreter','none');
end

function initDendrogramAppData(hFig)
    appdata = guihandles(hFig);
    appdata.highlightColor = [0 0 1];
    appdata.selectDot = [];
    appdata.hColGroupDot = [];
    appdata.hRowGroupDot = [];
    appdata.colDendroAxes = [];
    appdata.rowDendroAxes = [];
    appdata.colLines = [];
    appdata.rowLines = [];
   
    appdata.markerDatatip = [];
    setappdata(hFig, 'DendrogramData', appdata)
end

%----------------------------------------------------------
function disableUIControls(hFig)
% Disable all the UI controls
disableDendrogramUIControls(hFig)
disableColorMarkerUIControls(hFig)

%== Turn off interactive mode
activateuimode(hFig, '')

%== Disable window button function
set(hFig,'WindowButtonDownFcn',[]);
set(hFig,'WindowButtonUpFcn',[])
set(hFig,'WindowButtonMotionFcn',[]);
end

function disableDendrogramUIControls(hFig)
% Disable all the UI controls
appdata = getappdata(hFig, 'DendrogramData');

%== Disable group dot UI controls
if ~isempty(appdata.hColGroupDot)
    set(appdata.hColGroupDot, 'ButtonDownFcn',[],...
                              'UIContextMenu',[])
end
if ~isempty(appdata.hRowGroupDot)
    set(appdata.hRowGroupDot, 'ButtonDownFcn',[],...
                              'UIContextMenu', [])
end
appdata.dotContextMenu = [];
appdata.datatip = [];
appdata.highlight = [];
setappdata(hFig, 'DendrogramData', appdata)
end

function disableColorMarkerUIControls(hFig)
% Disable all the UI controls
appdata = getappdata(hFig, 'DendrogramData');

%== Color marker ui control
if ~isempty(appdata.colMarkers)
    set(appdata.colMarkers, 'ButtonDownFcn',[],...
                            'UIContextMenu', [])
end
if ~isempty(appdata.rowMarkers)
    set(appdata.rowMarkers, 'ButtonDownFcn',[],...
                            'UIContextMenu', [])
end
appdata.markerContextMenu = [];
appdata.markerDatatip = [];
setappdata(hFig, 'DendrogramData', appdata)
end


%-----------------------
function annotateTextCB(hSrc, ~, obj)
% % hFig= gcbf;
hFig = obj.FigureHandle;
state = bioinfoprivate.toggleState(hFig, hSrc);

switch state
    case 'on'
        obj.Annotate = true;
    case 'off'
        obj.Annotate = false;
end
end % end of function

function showDendrogramCB(hSrc, ~, obj)
hFig = obj.FigureHandle;
state = bioinfoprivate.toggleState(hFig, hSrc);
obj.ShowDendrogram = state;
end % end of function

function insertColorbarCB(hSrc, hEvt, obj) %#ok
hFig= gcbf;
state = bioinfoprivate.toggleState(hFig, hSrc);

switch state
    case 'on'
        obj.Colorbar = true;
    case 'off'
        obj.Colorbar = false;
end
end % end of function





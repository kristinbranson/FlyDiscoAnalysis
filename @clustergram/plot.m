function varargout = plot(obj, hFig, varargin)
%PLOT Render a clustergram heat map and dendrograms.
%
%   PLOT(CGOBJ) renders a heat map and dendrograms of a clustergram object
%   CGOBJ, in a MATLAB figure window.
%
%   PLOT(CGOBJ,HFIG) renders a heat map and dendrograms of a clustergram
%   object CGOBJ, in a MATLAB figure window with the handle HFIG.
%
%   H = PLOT(...) returns the handle to the figure. The graphic properties
%   are stored as application data in the figure handle.
%
%   Examples:
%
%       plot(cgobj)
%
%   See also CLUSTERGRAM, CLUSTERGRAM/VIEW.

%   Copyright 2007-2010 The MathWorks, Inc.
%   $Revision: 1.1.6.13 $  $Date: 2010/11/08 02:14:08 $

if isempty(obj)
    return;
end

%== Get figure handle
if nargin < 2
    hFig = [];
end
%== Remove all axes first
if nargin == 3 && ~isempty(varargin{1})
    imAxes = plot@HeatMap(obj, hFig, varargin{:});
else    
    if isempty(hFig) % Print to figure case
        hFig = figure('Renderer',     'ZBuffer',...
                      'Visible',      'on');
    end
    
    appdata = getappdata(hFig, 'DendrogramData');
    appdata.colDendroAxes = [];
    appdata.rowDendroAxes = [];
    appdata.colMarkerAxes = [];
    appdata.rowMarkerAxes = [];
    hRowLines = [];
    hColLines = [];
    colMarkers=[];
    appdata.colMarkerAnnotation = [];
    rowMarkers = [];
    appdata.rowlMarkerAnnotation = [];
    
    axesPos = obj.AxesPositions;
    %== Before calling super plot
    colMarkerFlag = false;
    rowMarkerFlag = false;
    
    if ~isempty(obj.ColumnGroupMarker)
        colMarkerFlag = true;
    end
    if ~isempty(obj.RowGroupMarker)
        rowMarkerFlag = true;
    end
    
    dendroAxesProps = {'Visible',         'off',...
                       'YTick',            [],...
                       'XTick',            [],...
                       'HandleVisibility', 'callback',...
                       'HitTest', 'off'};
    revAxesProps = {'YAxisLocation', 'right',...
                    'Xdir',          'reverse'};
    setappdata(hFig, 'DendrogramData', appdata);
        appdata = getappdata(hFig, 'DendrogramData');
        switch obj.Cluster
            case 'COLUMN' % 1
                appdata.rowDendroAxes = axes('Parent', hFig, dendroAxesProps{:},...
                    'Xlim', getLimit(obj.DendroRowLineX(obj.RowGroups, :)),...
                    'Ylim', getLimit2(numel(obj.RowLabels)),...
                    'Position', axesPos(2,:),...
                    revAxesProps{:});
                
                hRowLines = plotDendrogram(obj.DendroRowLineX(obj.RowGroups, :),...
                    obj.DendroRowLineY(obj.RowGroups, :),...
                    obj.DendroRowLineColor(obj.RowGroups, :),...
                    appdata.rowDendroAxes, 2,...
                    obj.ShowDendrogram);
                
                if rowMarkerFlag
                    appdata.rowMarkerAxes = axes('Parent', hFig,...
                        'Ylim', getLimit2(numel(obj.RowLabels)),...
                        'Position', axesPos(4,:),...
                        'HitTest', 'off',...
                        revAxesProps{:});
                    axis(appdata.rowMarkerAxes, 'off')
                end
            case 'ROW'
                appdata.colDendroAxes = axes('Parent', hFig, dendroAxesProps{:},...
                    'Xlim', getLimit2(numel(obj.ColumnLabels)),...
                    'Ylim', getLimit(obj.DendroColLineY(obj.ColGroups, :)),...
                    'Position', axesPos(3,:));
                
                hColLines = plotDendrogram(obj.DendroColLineX(obj.ColGroups, :),...
                    obj.DendroColLineY(obj.ColGroups, :),...
                    obj.DendroColLineColor(obj.ColGroups, :),...
                    appdata.colDendroAxes, 1,...
                    obj.ShowDendrogram);
                
                if colMarkerFlag
                    appdata.colMarkerAxes = axes('Parent', hFig,...
                                                 'Position', axesPos(5,:),...
                                                 'HitTest', 'off');
                    axis(appdata.colMarkerAxes, 'off')
                end
            case 'ALL'
                appdata.rowDendroAxes = axes('Parent', hFig, dendroAxesProps{:},...
                    'Xlim', getLimit(obj.DendroRowLineX(obj.RowGroups, :)),...
                    'Ylim', getLimit2(numel(obj.RowLabels)),...
                    'Position',axesPos(2,:),...
                    revAxesProps{:});
                
                appdata.colDendroAxes = axes('Parent', hFig, dendroAxesProps{:},...
                    'Xlim', getLimit2(numel(obj.ColumnLabels)),...
                    'Ylim', getLimit(obj.DendroColLineY(obj.ColGroups, :)),...
                    'Position',axesPos(3,:)); 
                
                hRowLines = plotDendrogram(obj.DendroRowLineX(obj.RowGroups, :),...
                    obj.DendroRowLineY(obj.RowGroups, :),...
                    obj.DendroRowLineColor(obj.RowGroups, :),...
                    appdata.rowDendroAxes, 2,...
                    obj.ShowDendrogram);
                
                hColLines = plotDendrogram(obj.DendroColLineX(obj.ColGroups, :),...
                    obj.DendroColLineY(obj.ColGroups, :),...
                    obj.DendroColLineColor(obj.ColGroups, :),...
                    appdata.colDendroAxes, 1,...
                    obj.ShowDendrogram);
                
                if rowMarkerFlag
                    appdata.rowMarkerAxes = axes('Parent', hFig,...
                        'Ylim', getLimit2(numel(obj.RowLabels)),...
                        'Position',axesPos(4,:),...
                        'HitTest', 'off',...
                        revAxesProps{:});
                    axis(appdata.rowMarkerAxes, 'off')
                end
                
                if colMarkerFlag
                    appdata.colMarkerAxes = axes('Parent', hFig,...
                                                 'Position', axesPos(5,:),...
                                                 'HitTest', 'off');
                    axis(appdata.colMarkerAxes, 'off')
                end
        end
        
        %== Plot heat map, calling HeatMap plot
        setappdata(hFig, 'DendrogramData', appdata)
        imAxes = plot@HeatMap(obj, hFig, varargin{:});
        %== Reposition the axes
        positionAxes(obj, imAxes)
        %== After plot heat amp
        appdata = getappdata(hFig, 'DendrogramData');
        
        %== Link Axes properties
        if ~isempty(appdata.colDendroAxes)
            set(hggetbehavior(appdata.colDendroAxes,'Zoom'),'Enable',false) ;
            linkaxes([imAxes appdata.colDendroAxes appdata.colMarkerAxes], 'x');
        end
        if ~isempty(appdata.rowDendroAxes)
            set(hggetbehavior(appdata.rowDendroAxes,'Zoom'), 'Enable',false);
            linkaxes([imAxes appdata.rowDendroAxes appdata.rowMarkerAxes], 'y');
        end
        
        if strcmpi(obj.Cluster, 'ALL')
            hXLimLink = linkprop([imAxes appdata.colDendroAxes appdata.colMarkerAxes], 'XLim');
            hYLimLink = linkprop([imAxes appdata.rowDendroAxes appdata.rowMarkerAxes], 'YLim');
            
            setappdata(imAxes,'ClustergramXLimLink',hXLimLink)
            setappdata(imAxes,'ClustergramYLimLink',hYLimLink)
        end
        
        %== Plot color markers
        if colMarkerFlag && ~strcmpi(obj.Cluster, 'COLUMN')
            [colMarkers, appdata.colMarkerAnnotation] = ...
                plotColorMarkers(obj, appdata.colMarkerAxes, hColLines, 2);
        end
        
        if rowMarkerFlag && ~strcmpi(obj.Cluster, 'ROW')
            rowTickText = getappdata(imAxes, 'YTickLabelTextHandles');
            if ~isempty(rowTickText)
                set(rowTickText, 'visible', 'off')
            end
            [rowMarkers, appdata.rowlMarkerAnnotation] = ...
                plotColorMarkers(obj, appdata.rowMarkerAxes, hRowLines, 1);
        end
        
        setappdata(hFig, 'DendrogramData', appdata)
        
        %== Save appdata
        if ishandle(obj.FigureHandle)
            appdata = getappdata(hFig, 'DendrogramData');
            appdata.colLineColor = [];
            appdata.rowGroupNames = [];
            appdata.colGroupNames =  [];
            appdata.rowLineColor = [];
            switch obj.Cluster
                case 'COLUMN'
                    appdata.rowLines = hRowLines;
                    appdata.rowGroupNames =  getGroupNames(obj.RowGroups);
                    appdata.rowLineColor = obj.DendroRowLineColor(obj.RowGroups,:);
                case 'ROW'
                    appdata.colLines = hColLines;
                    appdata.colLineColor = obj.DendroColLineColor(obj.ColGroups,:);
                    appdata.colGroupNames =  getGroupNames(obj.ColGroups);
                case 'ALL'
                    appdata.rowLines = hRowLines;
                    appdata.colLines = hColLines;
                    appdata.colLineColor = obj.DendroColLineColor(obj.ColGroups,:);
                    appdata.rowGroupNames =  getGroupNames(obj.RowGroups);
                    appdata.colGroupNames =  getGroupNames(obj.ColGroups);
                    appdata.rowLineColor = obj.DendroRowLineColor(obj.RowGroups,:);
            end
            set(appdata.rowDendroAxes, 'Tag', 'DendroRowAxes')
            set(appdata.colDendroAxes, 'Tag', 'DendroColAxes')
            set(appdata.colMarkerAxes, 'Tag', 'CMarkerColAxes')
            set(appdata.rowMarkerAxes, 'Tag', 'CMarkerRowAxes')
            
            appdata.colMarkers=colMarkers;
            appdata.rowMarkers = rowMarkers;
            
            setappdata(hFig, 'DendrogramData', appdata);

            obj.HMAxesHandle = imAxes;
        end
        
        %== Updtate colorbar position
        hmdata = getappdata(imAxes, 'HeatMapAxesData');
        if isempty(appdata.rowDendroAxes)
            left_pos = get(imAxes, 'Position');
        else
            left_pos = get(appdata.rowDendroAxes, 'Position');
        end
        hmdata.CBarLeftStart = left_pos(1);
        setappdata(imAxes, 'HeatMapAxesData', hmdata)
end
if bioinfoprivate.opttf(obj.Colorbar)
    plot@HeatMap(obj, imAxes, true); %plot colorbar
end

if nargout == 1
    varargout{1} = imAxes;
end
end %End of plot

%-------- Helper functions ---------
function hlines = plotDendrogram(xd, yd, cd, ax, dim, visFlag)
% Plot dendrogram to an axes and return the line handles

if dim == 1
    minv = min(min(xd));
    xd = xd -minv + 1;
elseif dim == 2
    minv = min(min(yd));
    yd = yd - minv + 1;
end

hlines = line(xd', yd', 'Parent', ax);
%== Set colors
for i = 1:size(cd, 1)
    set(hlines(i), 'Color', cd(i, :), 'Visible', visFlag)
    bioma.util.disableHGBehaviors(hlines(i))
end
end

function [hrecs, hanno]= plotColorMarkers(cg_obj, haxes, groupLineH, dim)
% Plot the color markers with annotation on top of the heatmap for column
% cluster (dimension 1)or to the right of the heatmap for row cluster
% (dimension 2). Dimension 1 color marker should have annotation on the
% marker is the space is fixed.  Dimension 2 color marker should have
% annotation to the right of the markers. Also Contextmenu for the color
% markers are: Change color, Remove.

hrecs = [];
hanno = [];
if dim == 2
    colorMarkers = cg_obj.ColumnGroupMarker;
    objGroups = cg_obj.ColGroups;
else
    colorMarkers = cg_obj.RowGroupMarker;
    objGroups = cg_obj.RowGroups;
end

cmgroups = zeros(length(colorMarkers), 1);
for i = 1:length(colorMarkers)
    cmgroups(i) = colorMarkers(i).GroupNumber;
end

markeridx = find(ismember(cmgroups, objGroups));

if isempty(markeridx)
    return;
end
nmarkers = numel(markeridx);

% Get the axes information
xl = get(haxes, 'Xlim');
xstart = xl(1);
width = diff(xl);
yl = get(haxes, 'Ylim');
height = diff(yl);
ystart = yl(1);
hrecs = zeros(nmarkers, 1);
hanno = zeros(nmarkers, 1);
for i = 1:nmarkers
    groupidx = find(objGroups == cmgroups(markeridx(i)));
    color = colorMarkers(markeridx(i)).Color;
    if isempty(color)
        color = 'b';
    end
    sel_groupidx = clusterPropagation(cg_obj, groupidx, dim);
    hcmgroups = groupLineH(sel_groupidx);
    
    if ~isempty(hcmgroups)
        if dim == 2 % col
            [startp, endp] = getGroupEnds(hcmgroups, 'xData');
            recPos = [xstart+(startp-1), ystart,...
                (endp-startp+1), height];
        else %row
            [startp, endp] = getGroupEnds(hcmgroups, 'yData');
            recPos = [xstart, ystart+(startp-1), width,...
                (endp-startp+1)];
            
            hanno(i) = text(xl(1)-width/5, ystart+(startp+endp-1)/2,...
                colorMarkers(markeridx(i)).Annotation,...
                'Parent', haxes,...
                'Clipping', 'off',...
                'VerticalAlignment', 'Middle',...
                'HorizontalAlignment', 'Left',...
                'FontSize', 10,...
                'FontWeight', 'bold',...
                'Visible', 'On');
            
        end
        hrecs(i) = rectangle('Parent', haxes, 'Position', recPos,...
            'FaceColor', color,...
            'LineStyle', 'none',...
            'DisplayName', colorMarkers(markeridx(i)).Annotation,...
            'Tag',num2str(markeridx(i)));
        set(get(get(hrecs(i),'Annotation'),'LegendInformation'),...
            'IconDisplayStyle','on');
    end
end
end

function [startp, endp] = getGroupEnds(hgroups, type)
% Helper function for plotColorMarker
xd = get(hgroups, type);
if size(xd, 1) > 1
    xd=cell2mat(xd);
    startp = fix(min(min(xd)));
    endp = fix(max(max(xd)));
else
    startp = fix(min(xd));
    endp = fix(max(xd));
end

end

function lim = getLimit(lineVals)
% Return the dendrogram axis limit determined by the dendrogram lines
lim = [min(min(lineVals)) max(max(lineVals))];
lim(2)=lim(2)+ceil(diff(lim))/75;
end

function lim = getLimit2(lineVals)
% Return the dendrogram axis limit determined by the dendrogram lines
ex=lineVals/75;
lim = [ex lineVals+ex];
end

function names = getGroupNames(gidx)
% Return cell array of group names
N = numel(gidx);
names = cell(N, 1);
for ind = 1:N
    names{ind} = ['Group ' num2str(gidx(ind))];
end
end
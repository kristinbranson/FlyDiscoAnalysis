classdef clustergram < HeatMap
%CLUSTERGRAM A heat map with dendrograms showing hierarchical clustering.
%   Clustergram performs a hierarchical clustering analysis of values in
%   the input data matrix and displays a heat map with row and column
%   dendrograms of the clustering. Typically the rows in the input matrix
%   are genes and the columns are samples.
%
%   Clustergram properties:
%       Cluster               - Dimension for clustering.
%       RowPDist              - Distance metric for compute pairwise distance between rows.
%       ColumnPDist           - Distance metric for compute pairwise distance between columns.
%       Linkage               - Linkage method to create the hierarchical cluster tree.
%       Dendrogram            - Color threshold property for dendrogram.
%       OptimalLeafOrder      - Logical flag to calculate optimal leaf ordering.
%       LogTrans              - Logical flag to log2 transform DATA from natural scale.
%       DisplayRatio          - Ratio between the spaces of dendrogram and heat map.
%       ColumnLabels          - Cell array of strings to label columns (x-axis).
%       RowLabels             - Cell array of strings to label rows (y-axis).
%       RowGroupMarker        - Structure array of information to annotate groups in row clusters.
%       ColumnGroupMarker     - Structure array of information to annotate groups in column clusters.
%       Standardize           - Direction in which the values are standardized.
%       ShowDendrogram        - Displaying dendrogram mode.
%       Colormap              - Colormap used to display the heat map.
%       DisplayRange          - The range of the values to be display.
%       Symmetric             - Logical flag to scale color symmetrically about zero.
%       ImputeFun             - Function name or handle for imputing missing values.
%       ColumnLabelsRotate    - Column labels orientation.
%       RowLabelsRotate       - Row labels orientation.
%       ColumnLabelsLocation  - Location of column labels.
%       RowLabelsLocation     - Location of row labels.
%       Annotate              - Logical flag to display value text in heat map.
%       AnnotPrecision        - Data value display precision.
%       AnnotColor            - Annotation text color
%       ColumnLabelsColor     - Structure array of information for coloring column labels.
%       RowLabelsColor        - Structure array of information for coloring row labels.
%       LabelsWithMarkers     - Logical flag to show colored markers for row/column labels.
%
%   Clustergram methods:
%       clustergram     - Create a clustergram object.
%       view            - Shows a clustergam object in a MATLAB figure.
%       plot            - Render a clustergran heat map and dendrograms.
%       clusterGroup    - Select a specified cluster group.
%       addXLabel       - Add clustergram x-axis (column) label.
%       addyLabel       - Add clustergram y-axis (row) label.
%       addTitle        - Add clustergram graph title.
%
%   Examples:
%       load filteredyeastdata;
%
%       % View the first 30 rows of genes.
%       cg = clustergram(yeastvalues(1:30,:));
%
%       % Add labels to the clustergram object cg.
%       set(cg, 'RowLabels', genes(1:30), 'ColumnLabels', times)
%
%       % Get properties of the clustergram object cg.
%       get(cg)
%
%       % Change clustering parameters
%       set(cg, 'Linkage', 'complete', 'Dendrogram', 3)
%
%       % Move the mouse over the dendrogram, left-click on a highlighted
%       % group, a datatip will show the group number and the names of the
%       % nodes under the group. Right-clicking on a highlighted group lets 
%       % you export the group in different ways.
%
%       % View all the data with a diverging red and blue colormap.
%       cg_all = clustergram(yeastvalues, 'Colormap', redbluecmap);
%
%       % Add annotation cluster group color markers to the clustergram
%       rm = struct('GroupNumber', {510, 593},...
%                   'Annotation', {'A', 'B'},...
%                   'Color', {'b', 'm'});
%
%       cm = struct('GroupNumber', {4,5},...
%                   'Annotation', {'Time1', 'Time2'},...
%                   'Color', {[1 1 0], [0.6 0.6 1]});
%
%       set(cg_all, 'RowGroupMarker', rm, 'ColumnGroupMarker', cm)
%
%   See also CLUSTERGRAMDEMO, DENDROGRAM, HEATMAP, LINKAGE, MAPCAPLOT,
%   PDIST, REDBLUECMAP, REDGREENCMAP, YEASTDEMO.

%   Copyright 2007-2009 The MathWorks, Inc.
%   $Revision: 1.1.6.10.4.19 $  $Date: 2010/12/09 23:16:36 $

properties(SetObservable=true, AbortSet=true)
    %CLUSTER Dimension for clustering.
    %    The Cluster property is a string (or a numeric value) specifying
    %    the dimension to cluster the data values. The dimension can be
    %    'COLUMN' (1) to cluster only the rows of the data along the
    %    columns, or 'ROW' (2) to cluster only the columns of the data
    %    along the rows, or the default 'ALL' (3) to cluster first the rows
    %    of data, then the columns.
    %
    %    See also CLUSTERGRAM.
    Cluster = 'ALL';
    
    %ROWPDIST Distance metric for compute pairwise distance between rows.
    %    The RowPDist property sets a distance metric for the function
    %    PDIST to use to compute the distance between rows. See the help
    %    for PDIST for more details of the available options. The default
    %    distance metric is 'Euclidean'. 
    %
    %    See also CLUSTERGRAM, PDIST.
    RowPDist = {'Euclidean'};
    
    %COLUMNPDIST Distance metric for compute pairwise distance between columns.
    %    The ColumnPDist property sets a distance metric for the function
    %    PDIST to use to compute the distance between columns. See the help
    %    for PDIST for more details of the available options. The default
    %    distance metric is 'Euclidean'. 
    %
    %    See also CLUSTERGRAM, PDIST.
    ColumnPDist = {'Euclidean'};
    
    %LINKAGE Linkage method to create the hierarchical cluster tree.
    %    The Linkage property sets the linkage method for the function
    %    LINKAGE to use to create the hierarchical cluster tree for both
    %    rows and columns.  See the help for LINKAGE for more details of
    %    the available options. It can be a string or a single-element cell
    %    array of strings used for both rows and columns, or a two-element
    %    cell array of strings, the first element is used for the linkage
    %    between rows, and the second element is used for the linkage
    %    between columns. To specify a linkage method for only one
    %    dimension, set the element for the other dimension to ''. The
    %    default method is 'average' for both dimensions. 
    %
    %    See also CLUSTERGRAM, LINKAGE.
    Linkage = {'Average'};
    
    %DENDROGRAM Color threshold property for dendrogram.
    %    The Dendrogram property sets the color threshold property for the
    %    function DENDROGRAM (the function used to create the dendrogram).
    %    See the help for DENDROGRAM for more details on the COLORTHRESHOLD
    %    option. The threshold can be a scalar used for both row and column
    %    dendrograms, or a two-element numeric vector or a cell array, the
    %    first element is used for the row dendrogram, and the second
    %    element is used for the column dendrogram. To specify a color
    %    threshold for only one dimension, set the element for the other
    %    dimension to '' or 0.
    %
    %    See also CLUSTERGRAM, DENDROGRAM.
    Dendrogram = {};
    
    %OPTIMALLEAFORDER Logical flag to calculate optimal leaf ordering.
    %   The OptimalLeafOrder property is a logical flag specifying the
    %   optimal leaf ordering calculation. When working with large data
    %   sets, calculating the optimal leaf ordering can be very time
    %   consuming and uses a large amount of memory. This option is
    %   disabled by default when the number of rows or columns is greater
    %   than 1000.
    %
    %    See also CLUSTERGRAM, OPTIMALLEAFORDER.
    OptimalLeafOrder = [];

    %LOGTRANS Logical flag to log2 transform DATA from natural scale.
    %   The LogTrans property is a logical flag specifying the log2
    %   transforming data from natural scale. By default, input data are
    %   assumed to be log2 based values.
    %
    %    See also CLUSTERGRAM.
    LogTrans = false;   
    
    %DISPLAYRATIO Ratio between the spaces of dendrogram and heat map.
    %   The DisplayRatio property sets the ratio of the space that the row
    %   and column dendrogram(s) occupy, relative to the width and height
    %   of the heat map. The ratio can be a single scalar value used as the
    %   ratio for both directions, or a two-element vector, the first
    %   element is used for the ratio of the row dendrogram width to the
    %   heat map width, and the second element is used for the ratio of the
    %   column dendrogram height to the heat map height. The second element
    %   is ignored for one-dimensional clustergrams. Default is 1/5.
    %
    %    See also CLUSTERGRAM.
    DisplayRatio = [1/5 1/5];

    %ROWGROUPMARKER Structure array of information to label cluster groups of rows.
    %   The RowGroupMarker property is a structure array of information to
    %   label the cluster groups of rows. The structure should contain
    %   these fields:
    %       GroupNumber
    %       Annotation
    %       Color
    %   GroupNumber is the row group number to label. Annotation is a
    %   string of text to annotate the row group and will be displayed next
    %   to the color marker. Color can be a string or three-element vector
    %   of RGB values specifying a color, which is used to label the row
    %   group. If this field is empty, default is 'blue'.
    %
    %    See also CLUSTERGRAM.
    RowGroupMarker = [];
    
    %COLUMNGROUPMARKER Structure array of information to label cluster groups of columns.
    %   The ColumnGroupMarker property is a structure array of information
    %   to label the cluster groups of columns. The structure should
    %   contain these fields:
    %       GroupNumber
    %       Annotation
    %       Color
    %   GroupNumber is the column group number to label. Annotation is a
    %   string of text to annotate the column group and will be displayed
    %   next to the color marker. Color can be a string or three-element
    %   vector of RGB values specifying a color, which is used to label the
    %   column group. If this field is empty, default is 'blue'.
    %
    %    See also CLUSTERGRAM.
    ColumnGroupMarker = [];   
    
    %SHOWDENDROGRAM Displaying dendrogram mode.
    %   The ShowDendrogram property determines whether to show dendrograms
    %   with heat map. The property value can be 'ON' (default) or 'OFF'.
    %
    %    See also CLUSTERGRAM.
    ShowDendrogram = 'on';
end

properties(GetAccess='private', SetAccess='private')
    FullWidth = 0.72;
    FullHeight = 0.75;
    
    RowCluster = [];
    ColCluster = [];
    DendroRowPerm = [];
    DendroColPerm = [];
    
    DendroRowLineX = [];
    DendroRowLineY = [];
    DendroRowLineColor = [];
    DendroColLineX = [];
    DendroColLineY = [];
    DendroColLineColor = [];

    ColGroups = [];
    RowGroups = [];
    ColNodes =[];
    RowNodes = [];
    NColGroups = 0;
    NRowGroups = 0;
end

properties(GetAccess='private', SetAccess='private', Dependent=true )
AxesPositions
end

methods
    function obj = clustergram(data, varargin)
        %CLUSTERGRAM Create a clustergram object.
        %
        %   C = CLUSTERGRAM(DATA) creates a clustergram C of hierarchical
        %   clustering of values in the matrix DATA and displays a heat map
        %   with row and column dendrograms of the clustering. The default
        %   clustering method is average linkage with Euclidean distance
        %   metric. DATA can be DataMatrix object or a MATLAB numeric
        %   matrix.
        %
        %   CLUSTERGRAM(...,'ROWLABELS',ROWLABELS) specifies the labels for
        %   the rows in the heat map. ROWLABELS can be a cell array of
        %   strings or a numeric array. 
        %
        %   CLUSTERGRAM(...,'COLUMNLABELS',COLUMNLABELS) specifies the
        %   labels for the columns in the heat map. COLUMNLABELS can be a
        %   cell array of strings or a numeric array. 
        %
        %   CLUSTERGRAM(...,'STANDARDIZE',DIM) specifies the dimension in
        %   which the values are standardized. The dimension can be
        %   'column' (1), 'row' (2), or the default 'none' (3).
        %
        %   CLUSTERGRAM(...,'CLUSTER',DIM) specifies the dimension for
        %   clustering the values in DATA. The dimension  can be 'COLUMN'
        %   (1), 'ROW' (2),  or ALL (3). The default dimension is 'ALL'.
        %
        %   CLUSTERGRAM(...,'ROWPDIST',DISTANCE) sets a distance metric for
        %   the function PDIST to use to compute the distance between rows.
        %   If the distance metric requires extra arguments, then these
        %   should be passed as a cell array. For example, to use the
        %   Minkowski distance with exponent P you would use {'minkowski',
        %   P}. See the help for PDIST for more details of the available
        %   options. The default distance metric is 'Euclidean'.
        %
        %   CLUSTERGRAM(...,'COLUMNPDIST',DISTANCE) sets a distance metric
        %   for the function PDIST to use to compute the distance between
        %   columns. It has the same options as ROWDIST. The default
        %   distance metric is 'Euclidean'.
        %
        %   CLUSTERGRAM(...,'LINKAGE', METHOD) sets the linkage method for
        %   the function LINKAGE to use to create the hierarchical cluster
        %   tree for both rows and columns. See the help for LINKAGE for
        %   more details of the available options.
        %
        %   CLUSTERGRAM(...,'DENDROGRAM',COLORTHR) sets the color threshold
        %   property for the function DENDROGRAM. See the help for
        %   DENDROGRAM for more details on the COLORTHRESHOLD option.
        %
        %   CLUSTERGRAM(...,'OPTIMALLEAFORDER',TF) disables the optimal
        %   leaf ordering calculation if TF is set to FALSE. This option is
        %   disabled by default when the number of rows or columns is
        %   greater than 1000. Set the value to TRUE to override this
        %   default.
        %
        %   CLUSTERGRAM(...,'COLORMAP',CMAP) allows you to specify the
        %   colormap used to display the clustergram heat map. It can be
        %   the name of a colormap, the function handle of a function that
        %   returns a colormap, or an M-by-3 array containing RGB values.
        %   The default is REDGREENCMAP.
        %
        %   CLUSTERGRAM(...,'DISPLAYRANGE', P) sets the display range of
        %   standardized values. P must be a positive scalar.  For example,
        %   if you specify REDGREENCMAP for the 'COLORMAP' property, pure
        %   red represents values equal to or greater than P, and pure
        %   green represents values equal to or less than -P. The default
        %   value is 3.
        %
        %   CLUSTERGRAM(...,'SYMMETRIC',TF) forces the color scale of the
        %   heat map to be symmetric about zero if TF is set to TRUE
        %   (default).
        %
        %   CLUSTERGRAM(...,'LOGTRANS', TR) will log2 transform DATA from
        %   natural scale. By default TR=FALSE and DATA are assumed to be
        %   log2 based values.
        %
        %   CLUSTERGRAM(...,'DISPLAYRATIO',R) sets the ratio of the space
        %   that the row and column dendrogram(s) occupy, relative to the
        %   width and height of the heat map. If R is a single scalar
        %   value, it is used as the ratio for both directions. If R is a
        %   two-element vector, the first element is used for the ratio of
        %   the row dendrogram width to the heat map width, and the second
        %   element is used for the ratio of the column dendrogram height
        %   to the heat map height. The second element is ignored for
        %   one-dimensional clustergrams. The default ratio is 1/5.
        %
        %   CLUSTERGRAM(...,'IMPUTEFUN',FUN) allows you to specify the name
        %   or function handle of a function that imputes missing data. FUN
        %   can also be a cell array with the first element being the
        %   function name or handle and other elements being the input
        %   property/value pairs for the function. The missing data points
        %   are colored gray in the heat map.
        %
        %   CLUSTERGRAM(...,'ROWGROUPMARKER', S) is an optional structure
        %   array for labeling the cluster groups in rows.  The structure
        %   should contain these fields:
        %       GroupNumber
        %       Annotation
        %       Color
        %   GroupNumber is the row group number to label. Annotation is a
        %   string of text to annotate the row group and will be displayed
        %   next to the color marker. Color can be a string or
        %   three-element vector of RGB values specifying a color, which is
        %   used to label the row group. If this field is empty, default is
        %   'blue'.
        %
        %   CLUSTERGRAM(...,'COLUMNGROUPMARKER', S) is an optional
        %   structure array for labeling the cluster groups of the columns.
        %   The structure should contain these fields:
        %       GroupNumber
        %       Annotation
        %       Color
        %   GroupNumber is the column group number to label. Annotation is
        %   a string of text to annotate the column group and will be
        %   displayed next to the color marker. Color can be a string or
        %   three-element vector of RGB values specifying a color, which is
        %   used to label the column groups. If this field is empty,
        %   default is 'blue'.
        %
        %   See also CLUSTERGRAM.
        
        if nargin < 1
            data = [];
        else
            % Check data to ensure 2 dimensional matrix
            if ndims(data)> 2 || ...
                    (~isnumeric(data) && ~isa(data, 'bioma.data.DataMatrix'))
                bioinfoprivate.bioclserror(mfilename, mfilename,...
                    'RequiresMatrixData', ...
                    'Input data must be a numeric matrix.')
            end
            inputval = varargin;
        end
        
        %== Construct object.
        obj = obj@HeatMap(data, false, 'RowLabelsLocation', 'right',...
                                       'DisplayRange', 3);
        if isempty(obj.Data)
            return;
        end
        
        %== Default optimalleafordering
        obj.OptimalLeafOrder = max(size(data)) <= 1000;
        
        % Check input options
        standardize = parseInputs(obj, inputval{:});
        
        % For copy object only. No need modify or view the data
        if obj.CopyOnly
            return;
        end
        
        %==Error on missing data and without Impute function
        if obj.MissingDataFlag && isempty(obj.ImputeFun)
            bioinfoprivate.bioclserror(mfilename, mfilename, 'MissingValue', ...
                ['Input data contains missing values.'...
                'Please provide an imputation function using the IMPUTEFUN option.'])
        end
        
        %== Validate parameters
        if isempty(obj.RowLabels)
            obj.RowLabels = cellstr(num2str((1:size(data,1))'));
        end
        
        if isempty(obj.ColumnLabels)
            obj.ColumnLabels = cellstr(num2str((1:size(data,2))'));
        end
        
        %== Set stanardize data option
        obj.Standardize = standardize;
        
        %== Compute cluster and dendrogram, update properties
        computeClusters(obj);
        
        %== Add listeners
        addPropertyListeners(obj);
        
        %== View the object
        obj.view;
    end % End of constructor
end % End of method block

methods
    function varargout = clusterGroup(obj, groupIdx, dim, varargin)
        %CLUSTERGROUP Select a specified cluster group.
        %
        %   CLUSTERGROUP(CG, GRPIDX, DIM) selects the specified cluster
        %   group GRPIDX in the clustergram window. It also highlights all
        %   the branches along the specified dimension DIM. GRPIDX must be
        %   a numeric group index. DIM can be 'COLUMN'(1) or 'ROW'(2).
        %
        %   GCG = CLUSTERGROUP(CG, GRPIDX, DIM) returns a new clustergram
        %   object, GCG, of the selected cluster group.
        %
        %   CLUSTERGROUP(..., 'COLOR', C) colors the dendrogram of
        %   specified cluster group with color C. C can be a three-element
        %   RGB vector or one of the predefined names. 
        %
        %   GINFO = CLUSTERGROUP(..., 'INFOONLY', TRUE) returns a structure
        %   containing information about the selected cluster group. 
        %   Default is FALSE.
        %
        %   See also CLUSTERGRAM.
        
        %== Check inputs
        bioinfochecknargin(nargin, 3, mfilename)
        
        if ~isnumeric(groupIdx) && isscalar(groupIdx)
            bioinfoprivate.bioclserror(mfilename, 'clusterGroup',...
                                      'InvalidDendroGroupIndex',...
                    'Invalid dendrogram group index. The group index must be a numeric scalar');
        end
        
        okdim = {'COLUMN','ROW'};
        if isnumeric(dim) && isscalar(dim)
            if dim < 1 || dim > 2
                bioinfoprivate.bioclserror(mfilename, 'clusterGroup',...
                                      'InvalidDimensionNumber',...
                    'Invalid dimension input. The dimension can be 1 (column) or 2 (row).');
            end

            dim = okdim{dim};
        end
        dim = bioma.util.validateStringArg(upper(dim), okdim,...
            'Dimension', 'clusterGroup', mfilename);
            
        switch dim
            case 'COLUMN'
                if ~ismember(groupIdx, obj.ColGroups)
                    bioinfoprivate.bioclserror(mfilename, 'clusterGroup',...
                                      'NotAGroupIndex',...
                    sprintf('Not such group index along column. The group range is %d and %d',...
                            min(obj.ColGroups), max(obj.ColGroups)));
                end
                fDim = 2;
            case 'ROW'
                 if ~ismember(groupIdx, obj.RowGroups)
                    bioinfoprivate.bioclserror(mfilename, 'clusterGroup',...
                                      'NotAGroupIndex',...
                    sprintf('Not such group index along row. The group range is %d and %d',...
                            min(obj.RowGroups), max(obj.RowGroups)));
                end
                fDim = 1;
        end
        
        inPV = parseDendroGroupInput(varargin{:});
        
        if inPV.InfoOnly
            if nargout > 0
                varargout{1} = getGroupInfo(obj, groupIdx, fDim);
            else
                return;
            end
        else
            if ishandle(obj.FigureHandle)
                view(obj,  groupIdx, fDim, inPV.Color)
            end
            
            if nargout > 0
                varargout{1} = getDendroGroupObject(obj, groupIdx, fDim);
            end
        end
    end
    
    function varargout = addTitle(obj, label, varargin) 
        %ADDTITLE  Add clustergram graph title.
        %
        %   ADDTITLE(CG, TITLESTR) adds the string, TITLESTR, at the top of
        %   the clustergram object CG graph display.
        %
        %   ADDTITLE(CG, TITLESTR, 'Property1', PropertyValue1,
        %   'Property2', PropertyValue2,...) sets the values of the
        %   specified properties of the title.
        %
        %   H = ADDTITLE(...) returns the handle to the text object used as
        %   the title.
        %
        %   See also CLUSTERGRAM, ADDXLABEL, ADDYLABEL.
        if nargin < 2
            ht = addTitle@HeatMap(obj);
        else
            ht = addTitle@HeatMap(obj, label, varargin{:});
        end
        
        if nargout > 0
            varargout{1} = ht;
        end
    end
    
    function varargout = addXLabel(obj, label, varargin)
        %ADDXLABEL  Add clustergram X-axis (column) label.
        %
        %   ADDXLABEL(CG, LABELSTR) adds the string, LABELSTR, beside the
        %   X-axis (column) of the clustergram object CG graph display.
        %
        %   ADDXLABEL(CG, LABELSTR, 'Property1', PropertyValue1,
        %   'Property2', PropertyValue2,...) sets the values of the
        %   specified properties of the X-axis label.
        %
        %   H = ADDXLABEL(...) returns the handle to the text object used
        %   as the X-axis label.
        %
        %   See also CLUSTERGRAM, ADDTITLE, ADDYLABEL.
        
        if nargin < 2
            ht = addXLabel@HeatMap(obj);
        else
            ht = addXLabel@HeatMap(obj, label, varargin{:});
        end 
        if nargout > 0
            varargout{1} = ht;
        end        
    end
    
    function varargout = addYLabel(obj, label, varargin)
        %ADDYLABEL  Add clustergram Y-axis (row) label.
        %
        %   ADDYLABEL(CG, LABELSTR) adds the string, LABELSTR, beside the
        %   Y-axis of the clustergram object CG graph display.
        %
        %   ADDYLABEL(CG, LABELSTR, 'Property1', PropertyValue1,
        %   'Property2', PropertyValue2,...) sets the values of the
        %   specified properties of the Y-axis label.
        %
        %   H = ADDYLABEL(...) returns the handle to the text object used
        %   as the y-axis label.
        %
        %   See also CLUSTERGRAM, ADDTITLE, ADDXLABEL.
        
        if nargin < 2
            ht = addYLabel@HeatMap(obj);
        else
            ht = addYLabel@HeatMap(obj, label, varargin{:});
        end 
        
        if nargout > 0
            varargout{1} = ht;
        end
    end
end

methods
    function positions = get.AxesPositions(obj)
        positions = getAxesPositions(obj);
    end
    
    function set.Cluster(obj, dim)
        okdim = {'COLUMN', 'ROW', 'ALL'};
        if isnumeric(dim) && isscalar(dim)
            if dim < 1 || dim > 3
                dim = 1;
            end
            
            dim = okdim{dim};
        end
        obj.Cluster = bioma.util.validateStringArg(upper(dim), okdim,...
            'Cluster', 'set', mfilename);
    end
    
    function set.RowPDist(obj, x)
        if iscell(x)
            obj.RowPDist = x;
        else
            obj.RowPDist = {x};
        end
    end
    
    function set.ColumnPDist(obj, x)
        if iscell(x)
            obj.ColumnPDist = x;
        else
            obj.ColumnPDist = {x};
        end
    end
    
    function set.Linkage(obj, x)
        if iscell(x)
            numx = numel(x);
            if numx ~= 1 && numx ~= 2
                bioinfoprivate.bioclserror(mfilename, 'set',...
                    'InvalidLinkageSize',...
                    'Linkage property must have 1 or 2 elements.');
            end
        end        
        obj.Linkage = x;
    end
    
    function set.Dendrogram(obj, x)
        if iscell(x) || isnumeric(x)
            numx = numel(x);
            if numx ~= 1 && numx ~= 2
                bioinfoprivate.bioclserror( mfilename, 'set',...
                    'InvalidDendrogramSize',...
                    'Dendrogram property must have 1 or 2 elements.');
            end
        end
        if ischar(x)
            obj.Dendrogram = {x};
        else
            obj.Dendrogram = x;
        end
    end
    
    function set.OptimalLeafOrder(obj, x)
        try
            obj.OptimalLeafOrder = bioinfoprivate.opttf(x, 'OptimalLeafOrder', 'set');
        catch ME
            bioinfoprivate.bioclsrethrow(mfilename, 'set', ME);
        end
    end
    
    function set.LogTrans(obj, x)
        try
            obj.LogTrans = bioinfoprivate.opttf(x,'LogTrans', 'set');
        catch ME
            bioinfoprivate.bioclsrethrow(mfilename, 'set', ME);
        end
    end
    
    function set.DisplayRatio(obj, x)
        numx = numel(x);
        if ~isa(x,'double')
            bioinfoprivate.bioclserror(mfilename, 'set', 'InvalidDisplayRatioInput',...
                'DISPLAYRATIO must be a double.');
        elseif numx ~= 1 && numx ~= 2
            bioinfoprivate.bioclserror(mfilename, 'set', 'InvalidDisplayRatioSize',...
                'DISPLAYRATIO must have 1 or 2 elements.');
        elseif ~all(x > 0 & x < 1);
            bioinfoprivate.bioclserror(mfilename, 'set', 'InvalidDisplayRatioValue',...
                'DISPLAYRATIO must be between 0 and 1.')
        end
        if numx == 1
            obj.DisplayRatio = [x x];
        else
            obj.DisplayRatio = x;
        end
    end
    
    function set.RowGroupMarker(obj, x)
        valid = checkColorMarkerStruct(x);
        if valid
            if numel(x) == 1 && numel(x.GroupNumber) > 1
               obj.RowGroupMarker = struct('GroupNumber', x.GroupNumber,...
                            'Annotation', x.Annotation,...
                            'Color', x.Color);
            else
                obj.RowGroupMarker = x;
            end
        else
            bioinfoprivate.bioclserror(mfilename, 'set', 'InvalidColorMarkerStructField',...
                'Invalid color marker structure field.')
        end
    end
    
    function set.ColumnGroupMarker(obj, x)
        valid = checkColorMarkerStruct(x);
        if valid
            if numel(x) == 1 && numel(x.GroupNumber) > 1
               obj.ColumnGroupMarker = struct('GroupNumber', x.GroupNumber,...
                            'Annotation', x.Annotation,...
                            'Color', x.Color);
                
            else
                obj.ColumnGroupMarker = x;
            end
        else
            bioinfoprivate.bioclserror(mfilename, 'set', 'InvalidColorMarkerStructField',...
                'Invalid color marker structure field.')
        end
    end
    
    function set.ShowDendrogram(obj, t)
        try
            tf = bioinfoprivate.opttf(t,'ShowDendrogram','set');
            if tf
                obj.ShowDendrogram = 'on';
            else
                obj.ShowDendrogram = 'off';
            end
        catch ME
            bioinfoprivate.bioclsrethrow(mfilename, 'set', ME);
        end
    end
end

methods(Hidden = true) 
    function data = getDisplayData(obj)
        if ~isempty(obj.DendroRowPerm) && ~isempty(obj.DendroColPerm)
            data = obj.PreOrderData(obj.DendroRowPerm, obj.DendroColPerm);
            
            switch obj.Cluster
                case 'COLUMN' % 1
                    data = data(obj.RowNodes, :);
                case 'ROW'
                    data = data(:, obj.ColNodes);
                case 'ALL'
                    data = data(obj.RowNodes, obj.ColNodes);
            end
        else
            data = obj.PreOrderData;
        end
    end
    
     function setDisplayData(obj, data) 
          setDisplayDataOnly(obj, data);
     end
    
    function data = getOriginalData(obj)
        % Overwrite by superclass method
        if ~isempty(obj.DendroRowPerm) && ~isempty(obj.DendroColPerm) && ~isempty(obj.PreOrderOriginalData)
            data = obj.PreOrderOriginalData(obj.DendroRowPerm, obj.DendroColPerm);
            
            switch obj.Cluster
                case 'COLUMN' % 1
                    data = data(obj.RowNodes, :);
                case 'ROW'
                    data = data(:, obj.ColNodes);
                case 'ALL'
                    data = data(obj.RowNodes, obj.ColNodes);
            end
        else
            data = obj.PreOrderOriginalData;
        end
    end
    
    function labels = getDimensionLabels(obj, dir)
        %Overwrite superclass HeatMap method
        switch dir
            case 1 % Row
                if ~isempty(obj.DendroRowPerm)
                    labels = obj.PreOrderRowLabels(obj.DendroRowPerm);
                    labels =  labels(obj.RowNodes);
                else
                    labels = obj.PreOrderRowLabels;
                end
            case 2 % Column
                if ~isempty(obj.DendroColPerm)
                    labels = obj.PreOrderColumnLabels(obj.DendroColPerm);
                    labels =  labels(obj.ColNodes);
                else
                    labels = obj.PreOrderColumnLabels;
                end
        end
    end
        
    function display(obj)
        n = numel(obj);
        if n > 1
            disp(obj)
        elseif n==0
            disp('    Empty array of clustergram objects')
        else
            [row, col] = size(obj.Data);
            switch obj.Cluster;
                case 'COLUMN'
                    if row < obj.NRowGroups+1 && row ~= 0
                        msg = sprintf('Clustergram object with %d(%d) rows of nodes.\n', row, obj.NRowGroups+1);
                    else
                        msg = sprintf('Clustergram object with %d rows of nodes.\n', row);
                    end
                case 'ROW'
                    if col < obj.NColGroups+1 && col ~= 0
                        msg = sprintf('Clustergram object with %d(%d) columns of nodes.\n', col, obj.NColGroups+1);
                    else
                        msg = sprintf('Clustergram object with %d columns of nodes.\n', col);
                    end
                case 'ALL'
                    if row == obj.NRowGroups+1 && col == obj.NColGroups+1 || row ==0 || col == 0
                       msg = sprintf('Clustergram object with %d rows of nodes and %d columns of nodes.\n',...
                                      row, col); 
                    else
                       msg = sprintf('Clustergram object with %d(%d) rows of nodes and %d(%d) columns of nodes.\n',...
                                    row, obj.NRowGroups+1, col, obj.NColGroups+1); 
                    end
            end
            disp(msg);
        end
    end
    
    function updateCluster(obj, src, evt) %#ok
        computeClusters(obj);
        view(obj);
    end
    
    function updateDisplay(obj, src, evt) %#ok
        switch src.Name
            case 'ShowDendrogram'
                view(obj, obj.ShowDendrogram)
            otherwise
                if ~obj.CopyOnly
                    view(obj);
                end
        end
    end
    
    function updateHMAxesProp(obj, src, evt)
        %== overwrite the HeatMap class method
        switch src.Name
            case 'ColumnLabelsLocation'
                obj.ColumnLabelsLocation = 'bottom';
                return;
            case 'RowLabelsLocation'
                obj.RowLabelsLocation = 'right';
                return;
            otherwise
                if ishandle(obj.FigureHandle)
                    updateHMAxesProp@HeatMap(obj, src, evt)
                end
        end
    end
    
    function computeClusters(obj)
        %COMPUTECLUSTERS Compute the clusters and updated clustergram properties.
        
        %obj.checkMissingData();
        if ~isempty(obj.ImputeFun)
          obj.imputeData(obj.PreOrderData);
        end
        
        %== Log2 transform data
        if obj.LogTrans
            obj.PreOrderData = log2(obj.PreOrderData);
        end
        
        %== Get clustergram properties by clustering direction
        [nRow, nCol] = size(obj.PreOrderData);
        data = obj.PreOrderData;
        pdistArgs = obj.RowPDist;
        dendroLocation = 'Left';
        dendroArgs = getDendroArgs(obj, 2);
        linkageArgs = getLinkageArgs(obj, 2);
        
        if strcmpi(obj.Cluster, 'ROW')
            data = obj.PreOrderData';
            pdistArgs = obj.ColumnPDist;
            dendroLocation = 'Top';
            dendroArgs = getDendroArgs(obj, 1);
            linkageArgs = getLinkageArgs(obj, 1);
        end
        %== Compute dendrogram along one direction
        [Z1, H1, ~, perm1] = computeDendrogram(data,...
            obj.OptimalLeafOrder,...
            dendroLocation,...
            pdistArgs,...
            linkageArgs,...
            dendroArgs);
        Z1 = Z1(:, 1:2);
        
        %== Compute dendrogram along other direction when needed
        if strcmpi(obj.Cluster, 'ALL')
            linkageArgs = getLinkageArgs(obj, 1);
            dendroArgs = getDendroArgs(obj, 1);
            [Z2, H2, ~, perm2] = computeDendrogram(data(perm1, :)',...
                obj.OptimalLeafOrder,...
                'Top',...
                obj.ColumnPDist,...
                linkageArgs,...
                dendroArgs);
            Z2 = Z2(:, 1:2);
        end
        %== Update clustergram properties
        switch obj.Cluster
            case 'COLUMN'
                obj.RowCluster = Z1;
                obj.DendroRowPerm = perm1;
                obj.DendroColPerm = 1:nCol;
                
                obj.NRowGroups = numel(H1);
                obj.RowNodes = 1:numel(H1)+1;
                obj.ColNodes = 1:nCol;
                obj.RowGroups = 1:numel(H1);
                [obj.DendroRowLineX,...
                    obj.DendroRowLineY,...
                    obj.DendroRowLineColor] = dendroLineInfo(H1);
            case 'ROW'
                obj.ColCluster = Z1;
                obj.DendroColPerm = perm1;
                obj.DendroRowPerm = 1:nRow;
                
                obj.NColGroups = numel(H1);
                obj.ColGroups = 1:numel(H1);
                obj.ColNodes = 1:numel(H1)+1;
                obj.RowNodes = 1:nRow;
                [obj.DendroColLineX,...
                    obj.DendroColLineY,...
                    obj.DendroColLineColor] = dendroLineInfo(H1);
            case 'ALL'
                obj.RowCluster = Z1;
                obj.DendroRowPerm = perm1;
                obj.NRowGroups = numel(H1);
                obj.RowGroups = 1:numel(H1);
                obj.RowNodes = 1:numel(H1)+1;
                [obj.DendroRowLineX,...
                    obj.DendroRowLineY,...
                    obj.DendroRowLineColor] = dendroLineInfo(H1);
                
                obj.ColCluster = Z2;
                obj.DendroColPerm = perm2;
                obj.NColGroups = numel(H2);
                obj.ColGroups = 1:numel(H2);
                obj.ColNodes =1:numel(H2)+1;
                [obj.DendroColLineX,...
                    obj.DendroColLineY,...
                    obj.DendroColLineColor] = dendroLineInfo(H2);
        end
        %== Close dendrogram figure opened by DENDROGRAM function
        close('ClustergramDendrogramFigure')
    end
    
    function [groupidx,nodeidx,groups,nodes] = clusterPropagation(obj, selidx, dim)
        % CLUSTERPROPAGATION Return the indices of sub groups and nodes.
        %
        % [GRPINDICES, NODEINDICES, GROUPS, NODES] =
        % CLUSTERPROPAGATION(OBJ, GRPIND, DIM) Propagates through the
        % cluster tree and find all the leaves below the selected group
        % GRPIND and returns found group and node indices in current group
        % and node properties of the clustergram object, and the group and
        % node indices in a vector of length ngroups + nnodes
        % [1:nnodes+ngroups]. Dim specifies the dimension of the group 1
        % for row and 2 for column.
        
        %== find nodes (children) below this group

        switch dim
            case 1
                perm = obj.DendroRowPerm;
                ngroups = obj.NRowGroups;
                cluster = obj.RowCluster;
                objGroups = obj.RowGroups;
                objNodes = obj.RowNodes;
            case 2
                perm = obj.DendroColPerm;
                ngroups = obj.NColGroups;
                cluster = obj.ColCluster;
                objGroups = obj.ColGroups;
                objNodes = obj.ColNodes;
        end
        
        groupidx = objGroups(selidx);
        
        nnodes = ngroups+1;
        children = false(nnodes + ngroups, 1);
        children(groupidx+nnodes) = true; % Selected
        for ind = ngroups:-1:1
            if children(ind+nnodes)
                children(cluster(ind,[1,2]))=true;
            end
        end
        
        groups = find(children(nnodes+1:nnodes+ngroups));
        nodes = children(1:nnodes);
        nodes = find(nodes(perm));
        
        groupidx = ismember(objGroups, groups);
        groupidx = find(groupidx);
        
        nodeidx = ismember(objNodes, nodes);
        nodeidx = find(nodeidx);
    end
    
    function newObj = createNewClustergramObj(obj)
        % Create a new clustergram object with the same properties.
        
        %== Create a new object of clustergram class
        newObj = clustergram(obj.PreOrderOriginalData,...
                             'CopyOnly', true,...
                             'RowLabels', obj.PreOrderRowLabels,...
                             'ColumnLabels', obj.PreOrderColumnLabels,...
                             'ShowDendrogram', obj.ShowDendrogram,...
                             'ImputeFun', obj.ImputeFun);
      
        %== Get meta-class of obj
        metaObj = metaclass(obj);
        %== Update new object's properties with obj's properties
        for iloop = 1 : length(metaObj.Properties)
            switch metaObj.Properties{iloop}.Name
                case 'RowLabels'
                case 'ColumnLabels'
                case 'CopyOnly'
                case 'AxesPositions'
                case 'ShowDendrogram'
                case 'ImputeFun'
                otherwise
                    newObj.(metaObj.Properties{iloop}.Name) = obj.(metaObj.Properties{iloop}.Name);
            end
        end
        
        newObj.FigureHandle = [];
        newObj.HMAxesHandle = [];
        addPropertyListeners(newObj)
    end
    %-------------------
    function positionAxes(obj, hHMAxes)
        %Reposition axes using HMAxes LooseInset so to show Xlabel, YLabel and
        %full tick labels
        appdata = getappdata(get(hHMAxes,'parent'), 'DendrogramData');
        if ~isfield(appdata, 'rowMarkerAxes')
            appdata.rowMarkerAxes = [];
            appdata.colDendroAxes = [];
            appdata.rowDendroAxes = [];
            appdata.colMarkerAxes = [];
        end
        lPos = get(hHMAxes, 'looseinset');
        
        hmPos = obj.AxesPositions(1,:);
        rowAPos = obj.AxesPositions(2,:);
        colAPos = obj.AxesPositions(3,:);
        rowMarkerAPos = obj.AxesPositions(4,:);
        colMarkerAPos = obj.AxesPositions(5,:);
        
        xDelta = hmPos(1) + hmPos(3)+lPos(3) - 1;
        if xDelta > 0 && isempty(appdata.rowMarkerAxes)
            %== Need to adjust the width and position
            ratio = obj.DisplayRatio(1);
            hmPos(1) = hmPos(1) - xDelta*ratio;
            hmPos(3) = hmPos(3) - xDelta*(1-ratio);
            rowAPos(3) = rowAPos(3) - xDelta*ratio;
            rowMarkerAPos(1) = rowMarkerAPos(1) - xDelta;
        end
        
        yDelta = lPos(2)-hmPos(2);
        if yDelta > 0
            ratio = obj.DisplayRatio(2);
            hmPos(2) = hmPos(2) + yDelta;
            hmPos(4) = hmPos(4) - yDelta*(1-ratio);
            
            colAPos(2) = colAPos(2) + yDelta*ratio;
            colAPos(4) = colAPos(4) - yDelta*ratio;
            colMarkerAPos(2) =  colMarkerAPos(2) + yDelta*ratio;
        end
        
        colAPos(1) = hmPos(1);
        colAPos(3) = hmPos(3);
        
        colMarkerAPos(1) = hmPos(1);
        colMarkerAPos(3) = hmPos(3);
        
        rowAPos(2) = hmPos(2);
        rowAPos(4) = hmPos(4);
        
        rowMarkerAPos(2) = hmPos(2);
        rowMarkerAPos(4) = hmPos(4);
        
        if ~isempty(hmPos) && all(hmPos(3:4) > 0)
            set(hHMAxes, 'Position', hmPos);
        end 
        
        if ~isempty(appdata.colDendroAxes) && ishandle(appdata.colDendroAxes) && all(colAPos(3:4)>0) 
            set(appdata.colDendroAxes, 'Position', colAPos);
        end
        
        if ~isempty(appdata.rowDendroAxes) && ishandle(appdata.rowDendroAxes) && all(rowAPos(3:4)>0) 
            set(appdata.rowDendroAxes, 'Position', rowAPos);
        end
        
        if ~isempty(appdata.colMarkerAxes) && ishandle(appdata.colMarkerAxes) && all(colMarkerAPos(3:4)>0)
            set(appdata.colMarkerAxes, 'Position', colMarkerAPos);
        end
        
        if ~isempty(appdata.rowMarkerAxes) && ishandle(appdata.rowMarkerAxes) && all(rowMarkerAPos(3:4)>0)
            rowTickText = getappdata(hHMAxes, 'YTickLabelTextHandles');
            if ~isempty(rowTickText)
                set(rowTickText, 'visible', 'off')
            end
            set(appdata.rowMarkerAxes, 'Position', rowMarkerAPos);
        end
        
        setappdata(get(hHMAxes,'parent'), 'DendrogramData', appdata);
    end
    %-------------
    function newcg_o = getDendroGroupObject(obj, groupIdx, dim)
        % Return a new clustergram object from the selected groups.
        % Propagate to get all the children under the group
        [sel_groupidx, ~, sel_groups, sel_nodes] =...
            clusterPropagation(obj, groupIdx, dim);
        
        %== create a new clustergram obj
        newcg_o = createNewClustergramObj(obj);
        
        if dim == 1 % Update Row data
            newcg_o.RowGroups = sel_groups;
            newcg_o.RowNodes = sel_nodes;
            
            if ~isempty(obj.RowGroupMarker)
                newcg_o.RowGroupMarker = getMarker(obj, obj.RowGroupMarker, sel_groups, dim);
            end
            
            if ishandle(obj.FigureHandle)
                appdata = getappdata(obj.FigureHandle, 'DendrogramData');
                newcg_o.DendroRowLineColor(sel_groups, :) = appdata.rowLineColor(sel_groupidx, :);
                if ~isempty(newcg_o.DendroColLineColor)
                    newcg_o.DendroColLineColor(obj.ColGroups, :) = appdata.colLineColor;
                end
            end
        elseif dim == 2
            newcg_o.ColGroups = sel_groups;
            newcg_o.ColNodes = sel_nodes;
            if ~isempty(obj.ColumnGroupMarker)
                newcg_o.ColumnGroupMarker = getMarker(obj, obj.ColumnGroupMarker, sel_groups, dim);
            end
            
            if ishandle(obj.FigureHandle)
                appdata = getappdata(obj.FigureHandle, 'DendrogramData');
                newcg_o.DendroColLineColor(sel_groups, :) = appdata.colLineColor(sel_groupidx, :);
                if ~isempty(newcg_o.DendroRowLineColor)
                    newcg_o.DendroRowLineColor(obj.RowGroups, :) = appdata.rowLineColor;
                end
            end
        end
    end
    %----------------
    
    function infoStruct = getGroupInfo(obj, groupIdx, dim)
        %Return a dendrogram group information structure with fields:
        % GroupNames
        % ColumnNodeNames
        % RowNodeNames
        % ExprValues
        appdata = getappdata(obj.FigureHandle, 'DendrogramData');
        [sel_groupidx, sel_nodeidx, ~, sel_nodes] =...
            clusterPropagation(obj, groupIdx, dim);
        
        rowOrder = obj.RowNodes;
        colOrder = obj.ColNodes;
        switch dim
            case 1
                infoStruct.GroupNames = appdata.rowGroupNames(sel_groupidx);
                infoStruct.RowNodeNames = obj.RowLabels(sel_nodeidx);
                infoStruct.ColumnNodeNames = obj.ColumnLabels;
                rowOrder =  sel_nodes;
            case 2
                infoStruct.GroupNames = appdata.colGroupNames(sel_groupidx);
                infoStruct.ColumnNodeNames = obj.ColumnLabels(sel_nodeidx);
                infoStruct.RowNodeNames = obj.RowLabels;
                colOrder = sel_nodes;
        end
        % Figure out if the query is from a subgroup clustergram, if it is
        % the originalData is adjusted.
        if all(size(obj.OriginalData) >= [max(rowOrder) max(colOrder)])
            infoStruct.ExprValues = obj.OriginalData(rowOrder, colOrder);
        else
            infoStruct.ExprValues = obj.OriginalData;
        end
    end
    
    function varargout = setdisp(obj, varargin)
            propertyNames = fieldnames(obj);
            
            propVals   = cell2struct(cell(size(propertyNames)), propertyNames, 1);
            propDescrs = cell2struct(cell(size(propertyNames)), propertyNames, 1);
            propVals.Cluster = '';
            propDescrs.Cluster = '[column | row | {all}]';
            propVals.RowPDist = [];
            propDescrs.RowPDist = '';
            propVals.ColumnPDist = [];
            propDescrs.ColumnPDist = '';
            propVals.Linkage = {};
            propDescrs.Linkage = '[single | complete | {average} | weighted | centroid | median | ward]';
            propVals.Dendrogram = [];
            propDescrs.Dendrogram = [];
            propVals.OptimalLeafOrder = [];
            propDescrs.OptimalLeafOrder = '[true | false]';
            propVals.LogTrans =  [];
            propDescrs.LogTrans = '[true | false]';
            propVals.DisplayRatio = [];
            propDescrs.DisplayRatio = [];
            propVals.RowGroupMarker = [];
            propDescrs.RowGroupMarker = 'A structure array.';
            propVals.ColumnGroupMarker = [];
            propDescrs.ColumnGroupMarker = 'A structure array.';
            propVals.ShowDendrogram = [];
            propDescrs.ShowDendrogram = '[true | false]';
            propVals.Standardize   = {};
            propDescrs.Standardize = '[column | {row} | none]';
            propVals.Symmetric   = {};
            propDescrs.Symmetric = '[true | false]';
            propVals.DisplayRange  = [];
            propDescrs.DisplayRange = 'Scalar';
            propVals.Colormap   = [];
            propDescrs.Colormap = [];
            propVals.ImputeFun   = [];
            propDescrs.ImputeFun = 'string -or- function handle -or- cell array';
            propVals.ColumnLabels   = {''};
            propDescrs.ColumnLabels = 'Cell array of strings, or an empty cell array';
            propVals.RowLabels   = {''};
            propDescrs.RowLabels = 'Cell array of strings, or an empty cell array';
            propVals.ColumnLabelsRotate   = [];
            propDescrs.ColumnLabelsRotate = [];
            propVals.RowLabelsRotate   = [];
            propDescrs.RowLabelsRotate = [];
            propVals.ColumnLabelsLocation   = [];
            propDescrs.ColumnLabelsLocation = '[ top | {bottom} ]';
            propVals.RowLabelsLocation   = [];
            propDescrs.RowLabelsLocation = '[ {left} | right ]';
            propVals.Annotate   = {};
            propDescrs.Annotate = '[on | {off}]';
            propVals.AnnotPrecision   = [];
            propDescrs.AnnotPrecision = [];
            propVals.AnnotColor   = {};
            propDescrs.AnnotColor = [];
            propVals.ColumnLabelsColor   = {};
            propDescrs.ColumnLabelsColor = 'A structure array.';
            propVals.RowLabelsColor   = {};
            propDescrs.RowLabelsColor = 'A structure array.';
            propVals.LabelsWithMarkers   = {};
            propDescrs.LabelsWithMarkers = '[true | false].';
            
            if nargin == 2
                try
                    name = bioinfoprivate.validateClassPropName(obj, varargin{1});
                    if nargout == 1
                        varargout{1} = propVals.(name);
                    else
                        disp(propDescrs.(name));
                    end
                catch ME
                    bioinfoprivate.bioclsrethrow(mfilename, 'set', ME);
                end
            else
                if nargout == 1
                    varargout{1} = propVals;
                else
                    disp(propDescrs);
                end
            end
        end
end % invisible methods block
end % End of classdef

%--------- Helper functions --------------------
function standardize = parseInputs(obj, varargin)
% Parse input PV pairs.
standardize = 'Row';
if nargin < 2
    return;
end
% Check for the right number of inputs
if rem(nargin-1, 2)== 1
    bioinfoprivate.bioclserror(mfilename, mfilename,...
                'IncorrectNumberOfArguments', ...
                'Incorrect number of arguments to clustergram.')
end
%== Backward compatibility checking
pvPairs = varargin;
obsoleteProps = {'PDist', 'Dimension', 'Standardize', 'CopyOnly',...
                 'SymmetricRange', 'Ratio', 'Impute', 'RowMarker', 'ColumnMarker'};
obsoleteIdx = true(length(pvPairs),1);
for iloop = 1:2:length(pvPairs)
    k = find(strncmpi(pvPairs{iloop}, obsoleteProps ,numel(pvPairs{iloop})));
    if ~isempty(k)
        p = pvPairs{iloop+1};
        switch k
            case 1 %PDist
                error(message('bioinfo:clustergram:PDistObsolete'))
            case 2 %Dimension
                error(message('bioinfo:clustergram:DimensionObsolete'))
            case 3 %Standardize
                standardize = p;
                pvPairs{iloop+1} = 'None';
            case 4 %CopyOnly
                obj.CopyOnly = p;
                obsoleteIdx(iloop) = false;
                obsoleteIdx(iloop + 1) = false;
            case 5 %SymmetricRange
                if numel(pvPairs{iloop}) == length('SymmetricRange')
                    error(message('bioinfo:clustergram:SymmetricRangeObsolete'))
                end
             case 6 %Ratio
                 if numel(pvPairs{iloop}) <= length('Ratio')
                     error(message('bioinfo:clustergram:RatioObsolete'))
                 end
            case 7 %Impute
                if numel(pvPairs{iloop}) <= length('Impute')
                    error(message('bioinfo:clustergram:ImputeObsolete'))
                end
             case 8 %RowMarker
                 error(message('bioinfo:clustergram:RowMarkerObsolete'))
             case 9 %ColumnMarker
                 error(message('bioinfo:clustergram:ColumnMarkerObsolete'))
        end
    end
end

pvPairs = bioinfoprivate.validateClassPVNames(obj, mfilename, pvPairs{obsoleteIdx});
if ~isempty(pvPairs)
    set(obj, pvPairs{:});
end
end

%--------------------
function [Z, lineH, T, Perm] = computeDendrogram(data, leafOrderFlag,...
                                                 dendroLoc, pdistArgs,...
                                                 linkageArgs, dendroArgs)
%Compute the hierarchical clustering and return dendrogram handles.

%== Create an invisible figure for dendrogram
hfig = figure('Visible','off', 'Name', 'ClustergramDendrogramFigure');

try
    if leafOrderFlag
        %= Calculate pairwise distances and linkage
        if numel(pdistArgs) == 1 && isnumeric(pdistArgs{1}),
          dist = pdistArgs{1};
        else
          dist = pdist(data, pdistArgs{:});
        end
        Z = linkage(dist, linkageArgs);
        order = optimalleaforder(Z, dist);
        clear('dist');
        
        if isempty(order) || (numel(unique(order)) ~= numel(order))
            warning(message('bioinfo:clustergram:optLeafOrderFailure', 'The optimal leaf ordering algorithm failed to converge.', 'Using default ordering.'))
            [lineH, T, Perm] = dendrogram(Z,0, dendroArgs{:},...
                                          'Orientation', dendroLoc); 
        else
            [lineH, T, Perm] = dendrogram(Z,0, dendroArgs{:},...
                                          'Orientation',dendroLoc, 'r', order);
        end
    else
        %= Calculate pairwise distances and linkage
        Z = linkage(data, linkageArgs, pdistArgs);
        [lineH, T, Perm] = dendrogram(Z,0, dendroArgs{:},...
                                      'Orientation',dendroLoc);
    end
catch ME
    delete(hfig)
    bioinfoprivate.bioclserror(mfilename, mfilename,'clusteringFailure',...
        'Encountered error while computing hierarchical clusters: \n%s', ME.message);
end
end

function [xd, yd, cd] = dendroLineInfo(hlines)
%Extract line xdata, ydata, and colors, and delete the lines.
n = numel(hlines);
xd = get(hlines, 'Xdata');
yd = get(hlines, 'Ydata');
cd = get(hlines, 'Color');

if n > 1
    xd = cell2mat(xd);
    yd = cell2mat(yd);
    cd = cell2mat(cd);
end
% Change single color to black
b = unique(cd, 'rows');
if size(b, 1) == 1
    cd = repmat([0 0 0], n, 1);
end

delete(hlines)
end

function dendroArgs = getDendroArgs(obj, dim)
% Return correct dendroArgs by dimension.
dendroArgs = {};
numdendroarg = numel(obj.Dendrogram);
if numdendroarg < 1
    return;
end

if dim == 1
    if numdendroarg > 1
        dendroArgs = obj.Dendrogram(2);
    else
        dendroArgs = obj.Dendrogram(1);
    end
elseif dim == 2
    dendroArgs = obj.Dendrogram(1);
end

% Dendrogram args only for ColorThreshold
if ~isempty(dendroArgs)
    if iscell(dendroArgs)
        dendroArgs = dendroArgs{:};
    end
    dendroArgs = {'ColorThreshold', dendroArgs};
end
end

function linkageArgs = getLinkageArgs(obj, dim)
% Return linkage args by dimension.
if iscell(obj.Linkage)
    numlinkagearg = numel(obj.Linkage);

    if dim == 1
        if numlinkagearg > 1
            linkageArgs = obj.Linkage{2};
        else
            linkageArgs = obj.Linkage{1};
        end
    elseif dim == 2
        linkageArgs = obj.Linkage{1};
    end
else
    linkageArgs = obj.Linkage;
end
end

function valid = checkColorMarkerStruct(cmStruct)
% Checks for color marker structure field names
if isempty(cmStruct)
    valid = true;
    return;
end

if ~isstruct(cmStruct)
    bioinfoprivate.bioclserror('clustergram', mfilename, 'InvalidColorMarkerInput',...
        'Color Marker must be a structure array or an empty value.');
end

valid =  isfield(cmStruct,'GroupNumber')&& ...
         isfield(cmStruct,'Annotation') && ...
         isfield(cmStruct,'Color');
end 
%--------------
function addPropertyListeners(obj)
addlistener(obj, 'Cluster', 'PostSet', @obj.updateCluster);
addlistener(obj, 'RowPDist', 'PostSet', @obj.updateCluster);
addlistener(obj, 'ColumnPDist', 'PostSet', @obj.updateCluster);
addlistener(obj, 'Linkage', 'PostSet', @obj.updateCluster);
addlistener(obj, 'Dendrogram', 'PostSet', @obj.updateCluster);
addlistener(obj, 'OptimalLeafOrder', 'PostSet', @obj.updateCluster);
addlistener(obj, 'LogTrans', 'PostSet', @obj.updateCluster);
addlistener(obj, 'DisplayRatio', 'PostSet', @obj.updateDisplay);
addlistener(obj, 'RowGroupMarker', 'PostSet', @obj.updateDisplay);
addlistener(obj, 'ColumnGroupMarker', 'PostSet', @obj.updateDisplay);
addlistener(obj, 'ShowDendrogram', 'PostSet', @obj.updateDisplay);
end
%--------------------
function axesPos = getAxesPositions(obj)
% Helper function to return positions for axes.
% fullW - fullWidth of the whole clustergram
% fullH - fullHeight of the whole clustergram
% dispRatio - dispRatio of the dendrograms [col(1), row(2)]
% clusterDim - Clustering dimension
% cmFlag - contains color marker if not return empty position

fullWidth = obj.FullWidth;
fullHeight = obj.FullHeight;
dispRatio = obj.DisplayRatio;
clusterDim = obj.Cluster;

cmFlag = [false false];

if ~isempty(obj.ColumnGroupMarker)
    cmFlag(1) = true;
end
if ~isempty(obj.RowGroupMarker)
    cmFlag(2) = true;
end

%== Initializing 
delta = 0.005; % gap between heatmap and dendrogram
vratio = 1/2; % 3/4for the starting point along Y
hratio = 1/3; % for the starting point along X

cmRatio = ([1/6 1/6.5].* dispRatio) .* cmFlag;
dispRatio = dispRatio - cmRatio;

%==
% Row 1 - HMAxes position
% Row 2 - Row Dendrogram
% Row 3 - Column Dendrogram
% Row 4 - Row Color marker axes
% Row 5 - Column Color marker Axes

axesPos = zeros(5,4);

switch clusterDim
    case 'COLUMN'
        dWidth = fullWidth * dispRatio(2);
        height = fullHeight*(1-dispRatio(1));
        dPos = [(1 - fullWidth)*hratio, (1 - fullHeight)*vratio, dWidth, height];
        
        imWidth = fullWidth-dWidth;
        imPos = [dPos(1) + dPos(3)+delta, dPos(2),imWidth,height];
        %For color marker
        cmWidth = fullWidth * cmRatio(2);
        cmPos = [imPos(1)+imWidth+delta imPos(2) cmWidth height];
        axesPos(1,:) = imPos;
        axesPos(2,:) = dPos;
        axesPos(4,:) = cmPos;
    case 'ROW'
        width = fullWidth;
        % For heatmap
        imHeight = fullHeight * (1 - dispRatio(1) - cmRatio(1));
        imPos = [(1 - fullWidth)/2.5, (1 - fullHeight)*vratio, width, imHeight];
        
        % For Color marker
        cmHeight = fullHeight * cmRatio(1);
        cmPos = [imPos(1) imPos(2)+imPos(4)+delta width cmHeight];
        % For dendrogram
        dHeight = fullHeight * dispRatio(1);
        dPos = [imPos(1) cmPos(2)+cmPos(4), width, dHeight];
        axesPos(1,:) = imPos;
        axesPos(3,:) = dPos;
        axesPos(5,:) = cmPos;
    case 'ALL'
        dWidth = fullWidth * dispRatio(2);
        height = fullHeight * (1-dispRatio(1)-cmRatio(1));
        
        dPos = [(1 - fullWidth)*hratio, (1 - fullHeight)*vratio, dWidth, height];
        
        imWidth = fullWidth-dWidth;
        imPos = [dPos(1) + dPos(3)+delta, dPos(2),imWidth,height];
        %For color marker
        cmWidth = fullWidth * cmRatio(2);
        cmPos = [imPos(1)+imWidth+delta imPos(2) cmWidth height];
        
        cmPos2 = [imPos(1) imPos(2)+imPos(4)+delta imPos(3) fullHeight*cmRatio(1)];
        dPos2 = [imPos(1), cmPos2(2)+cmPos2(4), imPos(3) fullHeight*dispRatio(1)];
        axesPos(1,:) = imPos;
        axesPos(2,:) = dPos;
        axesPos(3,:) = dPos2;
        axesPos(4,:) = cmPos;
        axesPos(5,:) = cmPos2;
end
end

%--------------------------------------------
function new_markers = getMarker(obj, old_markers, sel_groups, dim)
% For create new object and update its colormarkers
colorMarkers = old_markers;
cm_newidx = zeros(length(colorMarkers), 1);
for i = 1:length(colorMarkers)
    [~, ~, cm_groups] = clusterPropagation(obj, colorMarkers(i).GroupNumber, dim);
    if ~isempty(cm_groups)
        idx = find(ismember(cm_groups, sel_groups));
        if ~isempty(idx)
            cm_newidx(i) = i;
            colorMarkers(i).GroupNumber = cm_groups(idx(end));
        end
    end
end

new_markers = colorMarkers(cm_newidx ~= 0);
end

function inStruct = parseDendroGroupInput(varargin)
% Check for the right number of inputs
if rem(nargin,2)== 1
    bioinfoprivate.bioerror(mfilename,'IncorrectNumberOfArguments', ...
        'Incorrect number of arguments to %s.', 'clusterGroup')
end

% Allowed inputs
okargs = {'color', 'infoonly'};
inStruct.Color = [];           % Color dendrogram
inStruct.InfoOnly = false;      % get information only
for j=1:2:nargin
    [k, pval] = bioinfoprivate.pvpair(varargin{j}, varargin{j+1}, okargs, mfilename);
    switch(k)
        case 1 % Color
            if ischar(pval)
                try
                    inStruct.Color = bioinfoprivate.colorSpecLookUp(pval);
                catch ME
                    bioinfoprivate.bioclsrethrow(mfilename, 'clusterGroup', ME);
                end
            elseif isnumeric(pval) && isvector(pval) && max(size(pval))== 3 &&... 
                   max(pval) <= 1 && min(pval)>=0
                inStruct.Color = pval(:)';
            else
                bioinfoprivate.bioclserror(mfilename, 'clusterGroup', 'InvalidColor',...
                    'COLOR must be a three-element RGB vector or one of the predefined color names.');
            end 
        case 2 % InfoOnly
            inStruct.InfoOnly =  bioinfoprivate.opttf(pval,'InfoOnly','clusterGroup');
    end
end

end



function [handles,rawdata] = FlyBowlExamineVariables(vartype,varargin)

%% initialize outputs

% handles of everything
handles = struct('fig',1,'ax',[],...
  'data',[],'text',[],'manualpf',[],...
  'date',[],'selected',[],'selected1',[],....
  'xticklabels',[],'wait',[],...
  'menu',struct,'search',struct,'notes',struct);

% data pulled from sage
rawdata = [];

%% initialize state

% processed forms of rawdata for plotting
data = struct;
% experiment locations
data.expdirs = {};
data.expdir_bases = {};
data.nexpdirs = 0;
% line names
data.linenames = {};
data.lineidx = [];
data.nlines = 0;
% stats to plot
data.statnames = {};
data.nstats = 0;
% tables of stats per experiment
data.expstat = [];
data.stat = [];
data.normstat = [];
% registration data
data.registration_params = struct;
% date range
data.daterange = {};
% SAGE query info
data.queries = {};
% time data was pulled from SAGE
data.pull_data_datetime = [];
% data groups (e.g. line names)
data.groups = {};
% which group each experiment is in
data.groupidx = [];
% number of groups
data.ngroups = 0;
% normalization
data.mu = [];
data.sig = [];

% current GUI state
state = struct;
% whether we have made changes that need saving
state.needsave = false;
% list of annotations made for undoing
state.undolist = [];
% file to save annotations to
state.savefilename = '';
% name of mat file containing data
state.datafilename = '';
% selected group of data
state.datai_selected = [];
% selected statistic
state.stati_selected = [];
% data plots that are visible
state.idx_visible = [];
% number of days between today and last piece of data
state.mindaysprev = [];
% whether a subfunction is asking to quit
state.return = false;
% whether date range was input in function call
state.didinputdaterange = false;
% temporary files created
state.tempfilescreated = {};
% which groups are manual p, manual f, for quick ref
state.idx_manual_p = [];
state.idx_manual_f = [];
% which groups are automated f, for quick reg
state.idx_automated_f = [];

% parameters to gui
params = struct;
% type of variable we are plotting
params.vartype = vartype;
% possible groupings of data
params.groupnames = SetDefaultGroupNames();
% possible values for each group
params.groupvalues = cell(size(params.groupnames));
% analysis protocol, for reading settings files
params.analysis_protocol = 'current';
% settings root directory
params.settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
% name of dataloc parameters file
params.datalocparamsfilestr = 'dataloc_params.txt';
% root location of data
params.rootdatadir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data';
% date range stuff
params.date = struct;
% number of days to look at at a time
params.date.period = 7;
% max datenum to examine
params.date.maxdatenum = [];
% can be computed from maxdatenum, period
params.date.mindatenum = [];
params.date.maxdatestr = '';
params.date.mindatestr = '';
% current datenum when the gui was initialized
params.date.datenumnow = [];
% initial figure position
params.figpos = [];
% LDAP username
params.username = '';
% dataset to pull from SAGE flattened view
data.dataset = 'data';
% whether to load cached data
params.loadcacheddata = [];
% how to group the data for plotting
if strcmpi(params.vartype,'behavior'),
  params.plotgroup = 'line_name';
else
  params.plotgroup = 'none';
end
% data types to grab from SAGE (unset)
data.data_types = 0;
% extra stuff for restricting queries
params.query_leftovers = {};
% file locations
params.dataloc_params = struct;
% parameters in examine parameters file
params.examine_params = struct;
% flag values depend on var type

switch lower(params.vartype),
  case 'experiment',
    params.pvalue = 'P';
    params.fvalue = 'F';
    params.uvalue = 'U';
    params.manual_fn = 'manual_pf';
    params.manual_string = 'Manual P/F';
    params.manual_choices = {'Pass','Fail','Unknown'};
    params.docheckflags = false;
  case 'behavior',
    params.pvalue = 'N';
    params.fvalue = 'D';
    params.uvalue = 'U';
    params.manual_fn = 'manual_behavior';
    params.manual_string = 'Behavior';
    params.manual_choices = {'Normal','Different','Unknown'};
    params.docheckflags = true;
end

% whether to read missing data from file
params.readfromfile = true;

% annotations made
annot = struct;
% manual pf
annot.manual_pf = [];
% curation notes
annot.notes_curation = {};
% diagnostic info
annot.info = [];

% saved config
rc = struct('username','','savepath','','savedatapath','');

% data for plotting
plotdata = struct;
% colors per group
plotdata.colors = [];
% x location of each stat
plotdata.x = [];
% x-scale
plotdata.dx = [];
% y-scale
plotdata.dy = [];
% plot range
plotdata.xlim = [];
plotdata.ylim = [];

%% constants

constants = struct;

% max number of changes that can be undone
constants.maxnundo = 5;

% format for date options
constants.dateoptions_format = 'yyyy-mm-dd';

% format for Sage
constants.datetime_format = 'yyyymmddTHHMMSS';

% format for filenames
constants.filename_date_format = 'yyyymmdd';

% types of variables we can examine
constants.vartypes = {'experiment','behavior'};

% string metadata

constants.flag_metadata_fns = {'flag_redo','flag_review'};
constants.flag_transforms = [-.5,6];
constants.note_metadata_fns = {'notes_behavioral','notes_technical','notes_curation'};
constants.string_metadata_fns = {'line_name','manual_pf','automated_pf','experiment_name','experiment_protocol',...
  'experimenter','exp_datetime','bowl','camera','computer','harddrive','apparatus_id','cross_date','effector',...
  'environmental_chamber','gender','genotype','handler_sorting','handler_starvation','plate','rearing_protocol',...
  'rig','top_plate','file_system_path','manual_behavior','cross_handler'};
constants.manual_flag_fns = {'manual_pf','manual_behavior'};
constants.manual_flag_values = {{'U','P','F'},{'U','N','D'}};

%% default data_types to grab from Sage

if ~iscell(data.data_types),
  switch lower(params.vartype),
    case 'experiment',
      data.data_types = {'ufmf_diagnostics_summary_*',...
        'ctrax_diagnostics_*','registrationdata_*','sexclassifier_diagnostics_*',...
        'bias_diagnostics_*','temperature_diagnostics_*','bkgd_diagnostics_*'};
    case 'behavior',
      data.data_types = {'stats_perframe_*_meanmean_perexp','stats_perframe_*conditions','stats_perframe_*_nflies_analyzed'};
  end
end

%% parse inputs

if ~ismember(lower(params.vartype),constants.vartypes),
  error('Unknown variable type %s',params.vartype);
end

[params.analysis_protocol,params.settingsdir,params.datalocparamsfilestr,...
  handles.fig,params.date.period,params.date.maxdatenum,...
  params.figpos,params.date.datenumnow,params.username,params.rootdatadir,...
  data.dataset,params.loadcacheddata,state.datafilename,...
  params.groupnames,params.plotgroup,data.data_types,...
  params.docheckflags,params.flybowl_startdate,...
  params.SAGEdata,params.readfromfile,...
  params.query_leftovers] = ...
  myparse_nocheck(varargin,...
  'analysis_protocol',params.analysis_protocol,...
  'settingsdir',params.settingsdir,...
  'datalocparamsfilestr',params.datalocparamsfilestr,...
  'hfig',handles.fig,...
  'period',params.date.period,...
  'maxdatenum',params.date.maxdatenum,...
  'figpos',params.figpos,...
  'datenumnow',params.date.datenumnow,...
  'username',params.username,...
  'rootdatadir',params.rootdatadir,...
  'dataset',data.dataset,...
  'loadcacheddata',params.loadcacheddata,...
  'datafilename',state.datafilename,...
  'groupnames',params.groupnames,...
  'plotgroup',params.plotgroup,...
  'data_types',data.data_types,...
  'checkflags',params.docheckflags,...
  'flybowl_startdate','20110201',...
  'SAGEdata',[],...
  'readfromfile',params.readfromfile);

% number of date range options
params.date.maxperiodsprev = ceil((now-datenum(params.flybowl_startdate,'yyyymmdd'))/params.date.period);

%% read parameters

datalocparamsfile = fullfile(params.settingsdir,params.analysis_protocol,params.datalocparamsfilestr);
params.dataloc_params = ReadParams(datalocparamsfile);

paramsfileptr = sprintf('examine%svariablesparamsfilestr',lower(params.vartype));
examineparamsfile = fullfile(params.settingsdir,params.analysis_protocol,params.dataloc_params.(paramsfileptr));
params.examine_params = ReadParams(examineparamsfile);

examinefileptr = sprintf('examine%svariablesfilestr',lower(params.vartype));
examinefile = fullfile(params.settingsdir,params.analysis_protocol,params.dataloc_params.(examinefileptr));
data.examinestats = ReadExamineExperimentVariables(examinefile);
% abbr names of stats
SetStatNames();

registrationparamsfile = fullfile(params.settingsdir,params.analysis_protocol,params.dataloc_params.registrationparamsfilestr);
data.registration_params = ReadParams(registrationparamsfile);

%% load rc

if exist(params.examine_params.rcfile,'file'),
  rc = load(params.examine_params.rcfile);
  if ~isfield(rc,'savepath'),
    rc.savepath = '';
  end
  if ~isfield(rc,'savedatapath'),
    rc.savedatapath = '';
  end
end

%% select (default) date range

if isempty(params.SAGEdata),

SelectDateRange();

while true,

  %% ask questions:
  % user name, date range, whether to load data

  AskQuestions();
  if state.return,
    return;
  end

  %% load data from SAGE or cached mat file

  didgetdata = GetData();
  
  if state.return,
    return;
  end
  if didgetdata,
    break;
  end
  
end

else
  
  didgetdata = GetData(params.SAGEdata);
  if ~didgetdata,
    return;
  end

end

%% default file names

if strcmpi(vartype,'behavior'),
  prefix = 'ManualBehaviorAnnotation_FlyBowl';
  dataprefix = 'BehaviorData_FlyBowl';
else
  prefix = 'DataCuration_FlyBowl';
  dataprefix = 'ExperimentData_FlyBowl';
end
state.savefilename = fullfile(rc.savepath,sprintf('%s_%s_%sto%s.tsv',...
  prefix,rc.username,datestr(params.date.mindatenum,constants.filename_date_format),datestr(params.date.maxdatenum,constants.filename_date_format)));
if isempty(state.datafilename),
  state.datafilename = fullfile(rc.savedatapath,sprintf('%s_%sto%s.mat',...
    dataprefix,datestr(params.date.mindatenum,constants.filename_date_format),datestr(params.date.maxdatenum,constants.filename_date_format)));
end


%% add extra stuff

rawdata = AddExtraData(rawdata,data.statnames,data.registration_params);

%% modify groups of variables for searching to only be legal values

GetSearchGroups();

%% parse statistics

ParseStats();
annot.info = false(data.ngroups,data.nstats);

%% get mean and standard deviation for z-scoring

NormalizeStats();

%% set manual_pf = p, f marker, notes curation

InitializeManualPF();

%% create figure

CreateFigure();

%% plot data

PlotData();

%% create uicontrols

CreateUIControls();

%% rotate x-tick labels

RotateTickLabels();

%% set callbacks

SetCallbacks();

%% SUBFUNCTIONS %%

%% load data from SAGE or cached mat file

  function success = GetData(varargin)
    
    success = false;
    
    % try to load from cached data file
    if params.loadcacheddata,
      didloaddata = LoadData();
      if ~didloaddata,
        return;
      end
    elseif nargin >= 1,
      rawdata = varargin{1};
      save('TMP_FlyBowlExamineVariables_rawdata.mat','rawdata');      
    else
      
      % pull data from Sage
      data.queries = params.query_leftovers;
      data.queries(end+1:end+2) = {'daterange',data.daterange};
      data.queries(end+1:end+2) = {'data_type',data.data_types};
      data.queries(end+1:end+2) = {'flag_aborted',0};
      %data.queries(end+1:end+2) = {'automated_pf',{'P','U'}};
      data.queries(end+1:end+2) = {'experiment_name','FlyBowl_*'};
      data.pull_data_datetime = now;
      rawdata = SAGEGetBowlData(data.queries{:},'checkflags',params.docheckflags,'removemissingdata',true,'dataset',data.dataset,'analysis_protocol',params.analysis_protocol,...
        'settingsdir',params.settingsdir,'readfromfile',params.readfromfile);
      save('TMP_FlyBowlExamineVariables_rawdata.mat','rawdata');
      %tmp = load('TMP_FlyBowlExamineVariables_rawdata.mat','rawdata');
      %rawdata = tmp.rawdata;
      if isempty(rawdata),
        uiwait(warndlg(sprintf('No data for date range %s to %s',data.daterange{:}),'No data found'));
        return;
      end
    end

    % process data
    
    % sort by date
    date = {rawdata.exp_datetime};
    [~,order] = sort(date);
    rawdata = rawdata(order);
    
    % file system paths
    data.nexpdirs = numel(rawdata);
    data.expdir_bases = cell(1,data.nexpdirs);
    data.expdirs = cell(1,data.nexpdirs);
    for tmpi = 1:data.nexpdirs,
      data.expdir_bases{tmpi} = regexprep(rawdata(tmpi).experiment_name,'^FlyBowl_','');
      data.expdirs{tmpi} = fullfile(params.rootdatadir,data.expdir_bases{tmpi});
      rawdata(tmpi).file_system_path = data.expdirs{tmpi};
    end
    
    % make sure the following fields exist
    
    % notes_curation
    if ~isfield(rawdata,'notes_curation'),
      for tmpi = 1:data.nexpdirs,
        rawdata(tmpi).notes_curation = '';
      end
    else
      % this seems to be NULL
      for tmpi = 1:data.nexpdirs,
        if strcmpi(rawdata(tmpi).notes_curation,'NULL'),
          rawdata(tmpi).notes_curation = '';
        end
      end
      
    end
    
    % manual_behavior
    if ~isfield(rawdata,'manual_behavior'),
      for tmpi = 1:data.nexpdirs,
        rawdata(tmpi).manual_behavior = 'U';
      end
    end

    % manual_curator
    if ~isfield(rawdata,'manual_curator'),
      for tmpi = 1:data.nexpdirs,
        rawdata(tmpi).manual_curator = '';
      end
    else
      % this seems to be NULL
      for tmpi = 1:data.nexpdirs,
        if strcmpi(rawdata(tmpi).manual_curator,'NULL'),
          rawdata(tmpi).manual_curator = '';
        end
      end
      
    end
    
    % manual_curation_date
    if ~isfield(rawdata,'manual_curation_date'),
      for tmpi = 1:data.nexpdirs,
        rawdata(tmpi).manual_curation_date = '';
      end
    else
      % this seems to be NULL
      for tmpi = 1:data.nexpdirs,
        if strcmpi(rawdata(tmpi).manual_curation_date,'NULL'),
          rawdata(tmpi).manual_curation_date = '';
        end
      end
    end
    
    % get line name info
    [data.linenames,~,data.lineidx] = unique({rawdata.line_name});
    data.nlines = numel(data.linenames);
    
    % sort into plot groups
    if strcmpi(params.plotgroup,'none'),
      data.groups = data.expdir_bases;
      data.ngroups = data.nexpdirs;
      data.groupidx = 1:data.nexpdirs;
    else
      if ~isfield(rawdata,params.plotgroup),
        error('Plot group %s is not a field of data',params.plotgroup);
      end
      [data.groups,~,data.groupidx] = unique({rawdata.(params.plotgroup)});
      data.ngroups = numel(data.groups);
    end
    
    success = true;
    
  end

%% load data from file

  function didloaddata = LoadData()
    
    % fields to find in cached data:
    % rawdata
    % queries
    % pull_data_datetime
    % daterange
    
    didloaddata = false;
    
    try
      
      % load everything from file
      tmp = load(state.datafilename);

      % get rawdata
      if ~isfield(tmp,'rawdata'),
        error('rawdata variable not in %s',state.datafilename);
      end
      rawdata = tmp.rawdata;
      didloaddata = true;
      
      % try to get queries
      if isfield(tmp,'queries'),
        data.queries = tmp.queries;
      else
        % set to datafilename if no queries available
        data.queries = state.datafilename;
      end
      
      % try to get pull_data_datatime
      if isfield(tmp,'pull_data_datetime'),
        data.pull_data_datetime = tmp.pull_data_datetime;
      else
        % set to file timestamp if not available
        tmp = dir(state.datafilename);
        data.pull_data_datetime = tmp.datenum;
      end
      
      % try to get date range
      if isfield(tmp,'daterange'),
        data.daterange = tmp.daterange;
      else
        % use dates in experiments if not available
        tmpdates = sort({rawdata.exp_datetime});
        data.daterange = [tmpdates(1),tmpdates(end)];
      end
      
      % set date parameters based on date range
      params.date.maxdatenum = datenum(data.daterange{2},constants.datetime_format);
      params.date.mindatenum = datenum(data.daterange{1},constants.datetime_format);
      params.date.period = params.date.maxdatenum - params.date.mindatenum;
      state.mindaysprev = params.date.datenumnow-params.date.maxdatenum;

    catch ME
      uiwait(warndlg(sprintf('Could not load rawdata from %s:\n%s',state.datafilename,getReport(ME)),'Error loading data'));
    end
    
  end


%% set callbacks

  function SetCallbacks()

    set(handles.fig,'CloseRequestFcn',@close_fig_Callback);
    set(handles.ax,'ButtonDownFcn',@ButtonDownFcn);
    set(handles.fig,'WindowButtonMotionFcn',@MotionFcn,'BusyAction','queue');

  end

%% rotate x-tick labels

  function RotateTickLabels()

    handles.xticklabels = rotateticklabel(handles.ax,90);
    % make sure the ticks don't overlap the x-axis
    ex = get(handles.xticklabels(1),'Extent');
    y1 = ex(2)+ex(4);
    offy = y1-plotdata.ylim(1);
    for i = 1:numel(handles.xticklabels),
      pos = get(handles.xticklabels(i),'Position');
      pos(2) = pos(2) - offy;
      set(handles.xticklabels(i),'Position',pos);
    end
    set(handles.xticklabels,'Interpreter','none');
    
  end

%% create uicontrols

  function CreateUIControls()
    
    % text box
    set(handles.ax,'Units','normalized');
    axpos = get(handles.ax,'Position');
    textpos = [axpos(1),axpos(2)+axpos(4),axpos(3)/2,1-axpos(2)-axpos(4)-.01];
    handles.text = annotation('textbox',textpos,'BackgroundColor','k','Color','g',...
      'String','Experiment info','Interpreter','none');
    
    % manual pf
    
    set(handles.text,'Units','Pixels');
    textpos_px = get(handles.text,'Position');
    set(handles.text,'Units','normalized');
    margin = 5;
    w1 = 80;
    w2 = 100;
    w3 = 100;
    h1 = 20;
    h2 = 30;
    c1 = ((params.figpos(4)-margin) + textpos_px(2))/2;
    manualpf_textpos = [textpos_px(1)+textpos_px(3)+margin,c1-h1/2,w1,h1];
    handles.manualpf = struct;
    handles.manualpf.text = uicontrol(handles.fig,'Style','text','Units','Pixels',...
      'Position',manualpf_textpos,'String',[params.manual_string,':'],...
      'BackgroundColor',get(handles.fig,'Color'),'Visible','off');
    manualpf_popuppos = [manualpf_textpos(1)+manualpf_textpos(3)+margin,...
      c1-h2/2,w2,h2];
    handles.manualpf.popup = uicontrol(handles.fig,'Style','popupmenu','Units','Pixels',...
      'Position',manualpf_popuppos,'String',params.manual_choices,...
      'Value',3,'Visible','off','Callback',@manualpf_popup_Callback);
    manualpf_pushbuttonpos = [manualpf_popuppos(1)+manualpf_popuppos(3)+margin,...
      c1-h2/2,w3,h2];
    handles.manualpf.pushbutton = uicontrol(handles.fig,'Style','pushbutton','Units','Pixels',...
      'Position',manualpf_pushbuttonpos,'String','Add Note...',...
      'Callback',@manualpf_pushbutton_Callback,...
      'Visible','off');
    manualpf_pushbutton_info_pos = [manualpf_pushbuttonpos(1)+manualpf_pushbuttonpos(3)+margin,...
      c1-h2/2,w3,h2];
    handles.manualpf.pushbutton_info = uicontrol(handles.fig,'Style','pushbutton','Units','Pixels',...
      'Position',manualpf_pushbutton_info_pos,'String','Add Info...',...
      'Callback',@manualpf_pushbutton_info_Callback,...
      'Visible','off');
    set([handles.manualpf.text,handles.manualpf.popup,handles.manualpf.pushbutton,handles.manualpf.pushbutton_info],'Units','normalized');

    % date
    
    set(handles.ax,'Units','Pixels');
    axpos_px = get(handles.ax,'Position');
    set(handles.ax,'Units','normalized');
    w1 = 20;
    h1 = 20;
    w2 = 150;
    margin = 5;
    nextpos = [axpos_px(1)+axpos_px(3)-w1,axpos_px(2)+axpos_px(4)+margin,w1,h1];
    currpos = [nextpos(1)-margin-w2,nextpos(2),w2,h1];
    prevpos = [currpos(1)-margin-w1,nextpos(2),w1,h1];
    handles.date = struct;
    if state.mindaysprev <= 0,
      s = 'this week';
    elseif state.mindaysprev < 7,
      s = 'last week';
    else
      s = sprintf('%d weeks ago',ceil(state.mindaysprev/7));
    end
    handles.date.curr = uicontrol(handles.fig,'Style','text','Units','Pixels',...
      'Position',currpos,'String',s);
    handles.date.next = uicontrol(handles.fig,'Style','pushbutton','Units','Pixels',...
      'Position',nextpos,'String','>','Callback',@nextdate_Callback);
    handles.date.prev = uicontrol(handles.fig,'Style','pushbutton','Units','Pixels',...
      'Position',prevpos,'String','<','Callback',@prevdate_Callback);
    if state.mindaysprev < .0001,
      set(handles.date.next,'Enable','off');
    else
      set(handles.date.prev,'Enable','on');
    end
    set([handles.date.curr,handles.date.prev,handles.date.next],'Units','normalized');
    
    % create menus
    handles.menu = struct;
    handles.menu.file = uimenu('Label','File','Parent',handles.fig);
    handles.menu.options = uimenu('Label','Options','Parent',handles.fig);
    handles.menu.set = uimenu('Label','Set','Parent',handles.fig);
    handles.menu.info = uimenu('Label','Info','Parent',handles.fig);
    handles.menu.plot_manual_p = uimenu(handles.menu.options,'Label','Plot manual_pf = p',...
      'Checked','on','Callback',@plot_manual_p_Callback);
    handles.menu.plot_manual_f = uimenu(handles.menu.options,'Label','Plot manual_pf = f',...
      'Checked','on','Callback',@plot_manual_f_Callback);
    handles.menu.set_rest_manual_p = uimenu(handles.menu.set,'Label','Set manual_pf == u -> p',...
      'Callback',@set_rest_manual_p_Callback);
    handles.menu.set_rest_manual_f = uimenu(handles.menu.set,'Label','Set manual_pf == u -> f',...
      'Callback',@set_rest_manual_f_Callback);
    handles.menu.undo = uimenu(handles.menu.set,'Label','Undo',...
      'Callback',@undo_Callback,'Enable','off','Accelerator','z');
    handles.menu.save = uimenu(handles.menu.file,'Label','Save Spreadsheet...',...
      'Callback',@save_Callback,'Accelerator','s');
    handles.menu.load = uimenu(handles.menu.file,'Label','Load Spreadsheet...',...
      'Callback',@load_Callback,'Accelerator','l');
    handles.menu.save_data = uimenu(handles.menu.file,'Label','Save Data...',...
      'Callback',@savedata_Callback);
    handles.menu.load_data = uimenu(handles.menu.file,'Label','Load Data...',...
      'Callback',@loaddata_Callback);
    handles.menu.open = uimenu(handles.menu.file,'Label','View Experiment',...
      'Callback',@open_Callback,'Enable','off','Accelerator','o');
    handles.menu.search = uimenu(handles.menu.info,'Label','Search...',...
      'Callback',@search_Callback,'Accelerator','f');
    
    daterangeprint = {datestr(params.date.mindatenum,constants.dateoptions_format),datestr(params.date.maxdatenum,constants.dateoptions_format)};
    handles.menu.daterange = uimenu(handles.menu.info,'Label',...
      sprintf('Date range: %s - %s ...',daterangeprint{:}));
    handles.menu.username = uimenu(handles.menu.info,'Label',...
      sprintf('User: %s',rc.username));
    
  end


%% set manual_pf = p, f marker

  function InitializeManualPF()
    
    badidx = cellfun(@isempty,{rawdata.(params.manual_fn)});
    for tmpj = find(badidx),
      warning('No data for %s, experiment %s',params.manual_fn,strtrim(rawdata(tmpj).experiment_name));
    end
    annot.manual_pf = repmat(params.uvalue,[1,data.ngroups]);
    state.idx_automated_pf = false([1,data.ngroups]);
    if ~strcmpi(params.plotgroup,'none'),    
      for tmpi = 1:data.ngroups,
        s1 = lower([rawdata(data.groupidx==tmpi).(params.manual_fn)]);
        if ismember(lower(params.pvalue),s1) && ismember(lower(params.fvalue),s1),
          warning('Experiments for %s = %s are marked as both %s (%d/%d) and %s (%d/%d). Using unknown (%s).',...
            params.plotgroup,data.groups{tmpi},params.pvalue,nnz(s1==lower(params.pvalue)),numel(s1),params.fvalue,nnz(s1==lower(params.fvalue)),numel(s1),params.uvalue);
          annot.manual_pf(tmpi) = params.uvalue;
        elseif ismember(lower(params.pvalue),s1),
          annot.manual_pf(tmpi) = params.pvalue;
        elseif ismember(lower(params.fvalue),s1'),
          annot.manual_pf(tmpi) = params.fvalue;
        else
          annot.manual_pf(tmpi) = params.uvalue;
        end
        
        % same for automated_pf
        if isfield(rawdata,'automated_pf'),
          s1 = lower([rawdata(data.groupidx==tmpi).automated_pf]);
          state.idx_automated_f(tmpi) = all(s1 == 'f');
        end

        
      end
    else
      annot.manual_pf(~badidx) = [rawdata.(params.manual_fn)];
      state.idx_automated_f = lower([rawdata.automated_pf]) == 'f';
    end
    state.idx_manual_p = lower(annot.manual_pf) == lower(params.pvalue);
    state.idx_manual_f = lower(annot.manual_pf) == lower(params.fvalue);
    
    if ~strcmp(params.plotgroup,'none'),
      annot.notes_curation = cell(1,data.ngroups);
      for tmpi = 1:data.ngroups
        annot.notes_curation{tmpi} = unique({rawdata(data.groupidx==tmpi).notes_curation});
      end
    else
      annot.notes_curation = {rawdata.notes_curation};
    end

    
  end


%% plot data
  function PlotData()
    
    if data.ngroups == 1,
      off = 0;
    else
      off = (2*((1:data.ngroups)-(data.ngroups+1)/2)/(data.ngroups-1))*params.examine_params.offx;
    end
    
    handles.data = nan(1,data.ngroups);
    plotdata.x = nan(data.ngroups,data.nstats);
    for datai = 1:data.ngroups,
      plotdata.x(datai,:) = (1:data.nstats)+off(datai);
      handles.data(datai) = plot(handles.ax,plotdata.x(datai,:),data.normstat(datai,:),'o',...
        'color',plotdata.colors(datai,:),'markerfacecolor',plotdata.colors(datai,:),...
        'markersize',6,'HitTest','off');
    end
    state.idx_visible = true(1,data.ngroups);

    % set markers
    set(handles.data(state.idx_manual_p),'Marker','+');
    set(handles.data(state.idx_manual_f|(state.idx_automated_f&~state.idx_manual_p)),'Marker','x');

    
    plotdata.xlim = [0,data.nstats+1];
    miny = min(data.normstat(:));
    maxy = max(data.normstat(:));
    dy = maxy - miny;
    if dy == 0,
      maxy = miny + .001;
    end
    plotdata.ylim = [miny-.01*dy,maxy+.01*dy];
    plotdata.dx = diff(plotdata.xlim)/params.figpos(3);
    plotdata.dy = diff(plotdata.ylim)/params.figpos(4);
    set(handles.ax,'XLim',plotdata.xlim,'YLim',plotdata.ylim,'XTick',1:data.nstats,'XTickLabel',data.statnames,'XGrid','on');
    
    ylabel(handles.ax,'Stds from mean');
    
    % selected experiment
    handles.selected = plot(0,0,'o','color','k','Visible','off','HitTest','off','MarkerSize',10,'MarkerFaceColor','k');
    handles.selected1 = plot(0,0,'o','color','r','Visible','off','HitTest','off','MarkerSize',10,'MarkerFaceColor','r');
    
  end
    
%% create figure
  function CreateFigure()
    
    if ishandle(handles.fig),
      close(handles.fig);
    end
    figure(handles.fig);
    clf(handles.fig,'reset');
    if isempty(params.figpos),
      params.figpos = params.examine_params.figpos;
    end
    set(handles.fig,'Units','Pixels','Position',params.figpos,'MenuBar','none','ToolBar','figure');
    handles.ax = axes('Parent',handles.fig,'Units','Normalized','Position',params.examine_params.axespos);
    % plot 0
    plot(handles.ax,[0,data.nstats+1],[0,0],'k-','HitTest','off');
    hold(handles.ax,'on');
    plotdata.colors = jet(data.ngroups)*.7;
    drawnow;
    
  end


%% set stat names

  function SetStatNames()

    nstats = numel(data.examinestats);
    data.statnames = cell(1,nstats);
    for tmpi = 1:nstats,
      
      data.statnames{tmpi} = sprintf('%s_',data.examinestats{tmpi}{:});
      data.statnames{tmpi} = data.statnames{tmpi}(1:end-1);
      data.statnames{tmpi} = strrep(data.statnames{tmpi},'stats_perframe_','');
      data.statnames{tmpi} = strrep(data.statnames{tmpi},'meanmean_perexp_','');
      data.statnames{tmpi} = strrep(data.statnames{tmpi},'flyany_frame','');
      data.statnames{tmpi} = strrep(data.statnames{tmpi},'ufmf_diagnostics_summary','ufmf');
      data.statnames{tmpi} = strrep(data.statnames{tmpi},'ufmf_diagnostics_stream','ufmf');
      data.statnames{tmpi} = strrep(data.statnames{tmpi},'temperature_diagnostics','temp');
      data.statnames{tmpi} = strrep(data.statnames{tmpi},'bias_diagnostics','bias');
      data.statnames{tmpi} = strrep(data.statnames{tmpi},'bkgd_diagnostics','bkgd');
      data.statnames{tmpi} = strrep(data.statnames{tmpi},'ctrax_diagnostics','ctrax');
      data.statnames{tmpi} = strrep(data.statnames{tmpi},'registrationdata','reg');
      data.statnames{tmpi} = strrep(data.statnames{tmpi},'sexclassifier_diagnostics','sex');
      data.statnames{tmpi} = strrep(data.statnames{tmpi},'stats_perframe_','');
      data.statnames{tmpi} = strrep(data.statnames{tmpi},'flyany_frame','');
      data.statnames{tmpi} = strrep(data.statnames{tmpi},'_perexp','');
      
    end
    
  end

%% get mean and standard deviation for z-scoring
  function NormalizeStats()
    statsfileptr = sprintf('examine%svariablesstatsfilestr',params.vartype);
    if isempty(params.examine_params.(statsfileptr)),
      data.mu = nanmean(data.stat,1);
      data.sig = nanstd(data.stat,1,1);
    else
      examinevariablesstatsfile = fullfile(params.settingsdir,params.analysis_protocol,...
        params.examine_params.(statsfileptr));
      normstats = load(examinevariablesstatsfile);
      data.mu = nan(1,data.nstats);
      data.sig = nan(1,data.nstats);
      
      % data format used with behavior
      if isfield(normstats,'meanstatsperfly'),
        normstats = normstats.meanstatsperfly;
        for tmpi = 1:data.nstats,
          canload = true;
          if isempty(regexp(data.examinestats{tmpi}{1},'^stats_perframe','once')),
            canload = false;
          end
          if canload,
            if numel(data.examinestats{tmpi}) ~= 3,
              canload = false;
            end
          end
          if canload,
            fn = sprintf('%s_%s',data.examinestats{tmpi}{1},data.examinestats{tmpi}{3});
            fn = regexprep(fn,'^stats_perframe_','');
            canload = isfield(normstats,fn);
          end
          if canload,
            data.mu(tmpi) = normstats.(fn).meanmean;
            data.sig(tmpi) = normstats.(fn).stdmean;
          else
            data.mu(tmpi) = 0;
            data.sig(tmpi) = 1;
            warning('Not normalizing%s',sprintf(' %s',data.examinestats{tmpi}{:}));
          end
        end
        
        % data format used with experiment
      elseif isfield(normstats,'data.mu') && isfield(normstats,'data.sig'),
        
        for tmpi = 1:data.nstats,
          
          pathcurr = data.examinestats{tmpi};
          statcurr = normstats.data.mu;
          for tmpj = 1:numel(pathcurr),
            statcurr = statcurr.(pathcurr{tmpj});
          end
          data.mu(tmpi) = statcurr;
          statcurr = normstats.data.sig;
          for tmpj = 1:numel(pathcurr),
            statcurr = statcurr.(pathcurr{tmpj});
          end
          data.sig(tmpi) = statcurr;
          
        end
        
      end
    end
    
    data.z = data.sig;
    data.z(data.z == 0) = 1;
    data.normstat = bsxfun(@rdivide,bsxfun(@minus,data.stat,data.mu),data.z);

    
  end

%% parse statistics

  function ParseStats()
    
    data.nstats = numel(data.examinestats);
    data.stat = nan(data.ngroups,data.nstats);
    data.expstat = nan(data.nexpdirs,data.nstats);

    for tmpi = 1:data.nstats,
  
      if ~isfield(rawdata,data.examinestats{tmpi}{1}),
        warning('Examine stat %sis not in rawdata',sprintf('%s ',data.examinestats{tmpi}{:}));
        continue;
      end

      % special case: flags
      if ismember(data.examinestats{tmpi}{1},constants.flag_metadata_fns),
        % flags are now binary
        v = (double([rawdata.(data.examinestats{tmpi}{1})])-constants.flag_transforms(1))*constants.flag_transforms(2);
        data.expstat(:,tmpi) = v;
        % inefficient if group = experiment, but only done once
        for tmpj = 1:data.ngroups,
          data.stat(tmpj,tmpi) = nanmean(v(data.groupidx==tmpj));
        end

      % special case: notes
      elseif ismember(data.examinestats{tmpi}{1},constants.note_metadata_fns),
        
        v = cellfun(@(s) ~isempty(s) && ~strcmpi(s,'None') && ~strcmpi(s,'NULL'),{rawdata.(data.examinestats{tmpi}{1})});
        v = double(v)*2*3-3; % make -3, 3
        data.expstat(:,tmpi) = v;
        % inefficient if group = experiment, but only done once
        for tmpj = 1:data.ngroups,
          data.stat(tmpj,tmpi) = nanmean(v(data.groupidx==tmpj));
        end

      % special case: strings
      elseif ismember(data.examinestats{tmpi}{1},constants.string_metadata_fns),
        [uniquevals,~,v] = unique({rawdata.(data.examinestats{tmpi}{1})});
        if numel(uniquevals) == 1,
          v(:) = 0;
        else
          v = v - 1;
          v = (v/max(v)*2-1)*3;
        end
        data.expstat(:,tmpi) = v;
        % inefficient if group = experiment, but only done once
        for tmpj = 1:data.ngroups,
          data.stat(tmpj,tmpi) = nanmean(v(data.groupidx==tmpj));
        end
    
        % numbers
      else
    
        datacurr = {rawdata.(data.examinestats{tmpi}{1})};
        for k = 2:numel(data.examinestats{tmpi}),
          fn = data.examinestats{tmpi}{k};
          for tmpj = 1:numel(datacurr),
            if isstruct(datacurr{tmpj}) && isfield(datacurr{tmpj},fn),
              datacurr{tmpj} = datacurr{tmpj}.(fn);
            else
              datacurr{tmpj} = [];
            end
          end
        end
        badidx = cellfun(@isempty,datacurr);
        if any(badidx),
          for tmpj = find(badidx),
            datacurr{tmpj} = nan;
            warning('No data for stat %s experiment %s',sprintf('%s,',data.examinestats{tmpi}{:}),strtrim(rawdata(tmpj).experiment_name));
          end
        end
        badidx = ~cellfun(@isnumeric,datacurr);
        if any(badidx),
          for tmpj = find(badidx),
            datacurr{tmpj} = nan;
            warning('Non-numeric data for stat %s experiment %s',sprintf('%s,',data.examinestats{tmpi}{:}),strtrim(rawdata(tmpj).experiment_name));
          end
        end
        v = cell2mat(datacurr);
        data.expstat(:,tmpi) = v;
        % inefficient if group = experiment, but only done once
        for tmpj = 1:data.ngroups,
          data.stat(tmpj,tmpi) = nanmean(v(data.groupidx==tmpj));
        end
      end
    end
    
  end

%% modify groups of variables for seraching to only be legal values

  function GetSearchGroups()

    params.groupnames = intersect(params.groupnames,fieldnames(rawdata));
    params.groupvalues = cell(size(params.groupnames));
    for tmpi = 1:numel(params.groupnames),
      fn = params.groupnames{tmpi};
      if ismember(fn,constants.manual_flag_fns),
        tmpj = find(strcmp(fn,constants.manual_flag_fns),1);
        values = constants.manual_flag_values{tmpj};
      else
        values = {rawdata.(fn)};
        if all(cellfun(@ischar,values)),
          values = unique(values);
        else
          values = unique([values{:}]);
        end
      end
      params.groupvalues{tmpi} = values;
    end
  end

%% ask questions:
% user name, date range, whether to load rawdata

  function AskQuestions()
    
    havequestions = isempty(params.username) || isempty(params.loadcacheddata) || ...
      (params.loadcacheddata && (isempty(state.datafilename) || ~exist(state.datafilename,'file'))) || ...
      (~params.loadcacheddata && ~state.didinputdaterange);
    if havequestions,
      
      % username input
      if ~isempty(params.username),
        rc.username = params.username;
      end
      
      % date ranges allowed
      if state.didinputdaterange,
        daterange_strings = {sprintf('%s to %s',params.date.mindatestr,params.date.maxdatestr)};
        mindatenum_choices = datenum(params.date.mindatestr,constants.datetime_format);
        maxdatenum_choices = datenum(params.date.maxdatestr,constants.datetime_format);
      else
        % first day of week choices
        mindatenum_choices = fliplr([params.date.mindatenum-params.date.period*params.date.maxperiodsprev:params.date.period:params.date.mindatenum-params.date.period,params.date.mindatenum:params.date.period:params.date.datenumnow]);
        maxdatenum_choices = mindatenum_choices + params.date.period;
        mindatestr_choices = datestr(mindatenum_choices,constants.dateoptions_format);
        maxdatestr_choices = datestr(maxdatenum_choices,constants.dateoptions_format);
        daterange_strings = cellstr(cat(2,mindatestr_choices,...
          repmat(' to ',[numel(maxdatenum_choices),1]),...
          maxdatestr_choices));
      end
      
      % default date range
      weekidx = find(maxdatenum_choices == params.date.maxdatenum,1);
      [success,newweekidx,newusername,newdatafilename] = ExamineInitializeParams(daterange_strings,weekidx,rc.username,params.loadcacheddata,state.datafilename);
      if ~success,
        state.return = true;
        return;
      end
      if ~isempty(newusername),
        rc.username = newusername;
        params.username = newusername;
      end
      params.loadcacheddata = ~isempty(newdatafilename);
      if params.loadcacheddata,
        state.datafilename = newdatafilename;
        [rc.savedatapath] = fileparts(state.datafilename);
      end
      params.date.maxdatenum = maxdatenum_choices(newweekidx);
      params.date.mindatenum = mindatenum_choices(newweekidx);
      params.date.mindatestr = datestr(params.date.mindatenum,constants.datetime_format);
      params.date.maxdatestr = datestr(params.date.maxdatenum,constants.datetime_format);
      data.daterange = {params.date.mindatestr,params.date.maxdatestr};
      %maxdaysprev = params.date.datenumnow-params.date.mindatenum;
      state.mindaysprev = params.date.datenumnow-params.date.maxdatenum;
      
    end
    
  end

%% select date range
  function SelectDateRange()
    
    % if current time not set, use now
    if isempty(params.date.datenumnow),
      params.date.datenumnow = now;
    end
    
    % choose params.date.maxdatenum, by default, to be last Saturday
    state.didinputdaterange = ~isempty(params.date.maxdatenum);
    if ~state.didinputdaterange,
      params.date.maxdatenum = params.date.datenumnow;
      daycurr = weekday(params.date.maxdatenum);
      daywant = 7;
      params.date.maxdatenum = floor(params.date.maxdatenum)+daywant-daycurr-7;
    end
    
    % date range: params.date.maxdatenum-params.date.period to params.date.maxdatenum
    params.date.mindatenum = params.date.maxdatenum - params.date.period;
    params.date.mindatestr = datestr(params.date.mindatenum,constants.datetime_format);
    params.date.maxdatestr = datestr(params.date.maxdatenum,constants.datetime_format);
    data.daterange = {params.date.mindatestr,params.date.maxdatestr};
    
    % how many days from now is the params.date.period examined?
    state.mindaysprev = params.date.datenumnow-params.date.maxdatenum;
    
  end

%% CALLBACKS %%

%% update stat selected-based stuff

  function UpdateStatSelected()
    
    %set(hfig,'Interruptible','off');
    s = printfun(state.datai_selected,state.stati_selected);
    set(handles.text,'String',s); 
    set(handles.selected1,'XData',plotdata.x(state.datai_selected,state.stati_selected),...
      'YData',data.normstat(state.datai_selected,state.stati_selected),'visible','on');
    set(handles.menu.open,'Enable','on');
    set(handles.xticklabels,'Color','k','FontWeight','normal');
    set(handles.xticklabels(state.stati_selected),'Color','r','FontWeight','bold');

  end

%% update all info dependent on stat, experiment selected

  function UpdateSelected()
      
    UpdateStatSelected();
    
    set(handles.selected,'XData',plotdata.x(state.datai_selected,:),...
      'YData',data.normstat(state.datai_selected,:),'Visible','on');
    set(handles.menu.open,'Enable','on');
    
    manual_pf_curr = annot.manual_pf(state.datai_selected);
    s = get(handles.manualpf.popup,'String');
    vcurr = find(strncmpi(manual_pf_curr,s,1),1);
    if isempty(vcurr),
      error('Unknown manual_pf %s',manual_pf_curr);
    end
    set(handles.manualpf.popup,'Value',vcurr);
    
  end

%% print function: text to put in textbox about selected data group

  function s = printfun(datai,stati)

    expdiris = find(data.groupidx == datai);
    nexpdirscurr = numel(expdiris);
    
    if nexpdirscurr == 0,
      s = '';
      return;
    end
    
    % first line: group/experiment name
    s = data.groups(datai);

    if isempty(stati),
      if strcmpi(params.plotgroup,'none'),
        % if no stat seleced and no grouping, nothing else to say
      else
        % if no stat selected and there is grouping, then output the experiment
        % names
        s = [s,data.expdir_bases(expdiris)];
      end
      return;
    end
    
    % special case: notes, strings
    if ismember(data.examinestats{stati}{1},constants.note_metadata_fns) || ...
        ismember(data.examinestats{stati}{1},constants.string_metadata_fns),
      % convert cellstrs to single lines
      if iscell(rawdata(expdiris(1)).(data.examinestats{stati}{1})),
        s1 = cellfun(@(ss) sprintf('%s ',ss{:}),...
          rawdata(expdiris).(data.examinestats{stati}{1}),...
          'UniformOutput',false);
      else
        s1 = {rawdata(expdiris).(data.examinestats{stati}{1})};
      end
      % only one output value, then just show it
      if numel(unique(s1)) == 1,
        s{end+1} = sprintf('%s = %s',data.statnames{stati},s1{1});
      else
        s{end+1} = sprintf('%s =',data.statnames{stati});
        for expdirii = 1:numel(expdiris),
          expdiri = expdiris(expdirii);
          s{end+1} = sprintf('%s: %s',data.expdir_bases{expdiri},s1{expdirii}); %#ok<AGROW>
        end
      end
    else
      % numeric outputs
      v = data.stat(datai,stati);
      vz = data.normstat(datai,stati);
      v_exp = data.expstat(expdiris,stati);
      % fix flag normalization
      if ismember(data.examinestats{stati}{1},constants.flag_metadata_fns),
        v = v/constants.flag_transforms(2)+constants.flag_transforms(1);
      end
      % only one output value, then just show it
      if numel(unique(v_exp)) == 1,
        s{end+1} = sprintf('%s = %s = %s std',data.statnames{stati},num2str(v),num2str(vz));
      else
        s{end+1} = sprintf('mean(%s) = %s = %s std',data.statnames{stati},num2str(v),num2str(vz));
        for expdirii = 1:numel(expdiris),
          expdiri = expdiris(expdirii);
          s{end+1} = [data.expdir_bases{expdiri},': ',num2str(v_exp(expdirii))]; %#ok<AGROW>
        end
      end
    end
  end

%% button down function

  function ButtonDownFcn(varargin)
        
    try
      
      if ~isempty(state.stati_selected),
        set(handles.xticklabels(state.stati_selected),'Color','k','FontWeight','normal');
      end

      
      % get current point
      tmp = get(handles.ax,'CurrentPoint');
      xclicked = tmp(1,1);
      yclicked = tmp(1,2);
      tmpidx = find(state.idx_visible);
      tmpn = numel(tmpidx);
      [d,closest] = min( ((reshape(plotdata.x(state.idx_visible,:),[1,tmpn*data.nstats])-xclicked)/plotdata.dx).^2+...
        ((reshape(data.normstat(state.idx_visible,:),[1,tmpn*data.nstats])-yclicked)/plotdata.dy).^2 );
      d = sqrt(d);
      
      if d > params.examine_params.maxdistclick,
        set(handles.selected,'Visible','off');
        set(handles.selected1,'Visible','off');
        set(handles.menu.open,'Enable','off');
        state.stati_selected = [];
        state.datai_selected = [];
        set(handles.text,'String','');
        set([handles.manualpf.text,handles.manualpf.popup,...
          handles.manualpf.pushbutton,handles.manualpf.pushbutton_info],'Visible','off');
        return;
      end
      
      SelectionType = get(handles.fig,'SelectionType');
      set([handles.manualpf.text,handles.manualpf.popup,handles.manualpf.pushbutton,handles.manualpf.pushbutton_info],'Visible','on');

      [state.datai_selected,state.stati_selected] = ind2sub([tmpn,data.nstats],closest);
      state.datai_selected = tmpidx(state.datai_selected);

      UpdateSelected();
      
      if strcmp(SelectionType,'open'),

        open_Callback();

      end
      
    catch ME,

      fprintf('Error evaluating buttondownfcn, disabling:\n');
      set(handles.ax,'ButtonDownFcn','');
      rethrow(ME);
    end
    
  end

%% mouse moved callback

  function MotionFcn(hObject,event) %#ok<INUSD>
    
    if isempty(state.datai_selected),
      return;
    end
    
    try
      tmp = get(handles.ax,'CurrentPoint');
      xhover = tmp(1,1);
      if xhover < 0 || xhover > data.nstats+1,
        return;
      end
      state.stati_selected = min(data.nstats,max(1,round(xhover)));
      UpdateStatSelected();
    catch ME,
      fprintf('Error evaluating motionfcn, disabling:\n');
      set(handles.fig,'WindowButtonMotionFcn','');
      rethrow(ME);
    end
    
  end

%% save rc file

  function SaveConfigFile()

    try
      save(params.examine_params.rcfile,'-struct','rc');
    catch ME,
      warning('Error saving rc file:\n %s',getReport(ME));
    end

  end
%% clean up

  function close_fig_Callback(hObject,event)
    
    if state.needsave,
      answer = questdlg(sprintf('Save %s and notes?',params.manual_fn));
      switch lower(answer),
        case 'cancel',
          return;
        case 'yes',
          save_Callback(hObject,event);
      end
    end


    SaveConfigFile();
    
    % delete temporary files
    try
    for tmpi = 1:numel(state.tempfilescreated);
      if exist(state.tempfilescreated{tmpi},'file'),
        delete(state.tempfilescreated{tmpi});
      end
    end
    catch ME
      warning('Error deleting temporary files:\n %s',getReport(ME));
    end
    
    % delete search box
    if isfield(handles.search,'dialog') && ishandle(handles.search.dialog),
      delete(handles.search.dialog);
    end
    
    % delete notes dialog
    if isfield(handles.notes,'dialog') && ishandle(handles.notes.dialog),
      close(handles.notes.dialog);
    end

    if ishandle(hObject),
      delete(hObject);
    end
    
  end


%% open experiment/group view

  function open_Callback(hObject,event) %#ok<INUSD>
    % open experiment
    if isempty(state.datai_selected),
      return;
    end

    expdiris = find(data.groupidx == state.datai_selected);
    if isempty(expdiris),
      return;
    end
    
    % only one experiment, then just open the folder
    if numel(expdiris) == 1,
      if ispc,
        winopen(data.expdirs{expdiris});
      else
        web(data.expdirs{expdiris},'-browser');
      end
      return;
    end
    
    % create an html page with links to all experiments
    filenamecurr = fullfile(tempdir,[data.groups{state.datai_selected},'.html']);
    fid = fopen(filenamecurr,'w');
    if fid >= 0,
      fprintf(fid,'<html>\n<title>%s %s</title>\n<body>\n',params.plotgroup,data.groups{state.datai_selected});
      fprintf(fid,'<h1>%s</h1>\n',data.groups{state.datai_selected});
      fprintf(fid,'<ul>\n');
      for expdiri = expdiris(:)',
        fprintf(fid,'  <li><a href="file://%s">%s</a></li>\n',data.expdirs{expdiri},data.expdir_bases{expdiri});
      end
      fprintf(fid,'</ul>\n');
      fprintf(fid,'</body>\n</html>\n');
      fclose(fid);
    end
    if ~exist(filenamecurr,'file'),
      warning('Could not open temporary file %s',filenamecurr);
      return;
    end
    state.tempfilescreated{end+1} = filenamecurr;
    % open this page
    web(filenamecurr,'-browser');
  end

%% manual_pf = p callback

  function plot_manual_p_Callback(hObject,event) %#ok<INUSD>
    
    v = get(hObject,'Checked');
    set(handles.data(state.idx_manual_p),'Visible',v);
    if strcmpi(v,'on'),
      set(hObject,'Checked','off');
      state.idx_visible(state.idx_manual_p) = false;
      set(handles.data(~state.idx_visible),'Visible','off');
    else
      set(hObject,'Checked','on');
      state.idx_visible(state.idx_manual_p) = true;
      set(handles.data(state.idx_visible),'Visible','on');
    end
    
  end

%% manual_pf = f callback

  function plot_manual_f_Callback(hObject,event) %#ok<INUSD>
    
    v = get(hObject,'Checked');
    set(handles.data(state.idx_manual_f|(state.idx_automated_f&~state.idx_manual_p)),'Visible',v);
    if strcmpi(v,'on'),
      set(hObject,'Checked','off');
      state.idx_visible(state.idx_manual_f|(state.idx_automated_f&~state.idx_manual_p)) = false;
      set(handles.data(~state.idx_visible),'Visible','off');
    else
      set(hObject,'Checked','on');
      state.idx_visible(state.idx_manual_f|(state.idx_automated_f&~state.idx_manual_p)) = true;
      set(handles.data(state.idx_visible),'Visible','on');
    end
    
  end
  
  %% go to next date range

  function nextdate_Callback(hObject,event)
    
    if state.needsave,
      answer = questdlg('Save manual_pf and notes?');
      switch lower(answer),
        case 'cancel',
          return;
        case 'yes',
          save_Callback(hObject,event);
      end
    end
    
    if isfield(handles.search,'dialog') && ...
        ishandle(handles.search.dialog),...
        delete(handles.search.dialog);
    end
    
    params.loadcacheddata = false;
    params.date.maxdatenum = params.date.maxdatenum + params.date.period;
    
    RecursiveCall();

  end

%% go to previous date range

  function prevdate_Callback(hObject,event)
    
    if state.needsave,
      answer = questdlg('Save manual_pf and notes?');
      switch lower(answer),
        case 'cancel',
          return;
        case 'yes',
          save_Callback(hObject,event);
      end
    end
    
    if isfield(handles.search,'dialog') && ...
        ishandle(handles.search.dialog),...
        delete(handles.search.dialog);
    end
    
    params.loadcacheddata = false;
    params.date.maxdatenum = params.date.maxdatenum - params.date.period;
    
    RecursiveCall();

  end

%% call this function again, maybe with different state

  function RecursiveCall()

    SaveConfigFile();
    [handles,rawdata] = FlyBowlExamineVariables(params.vartype,...
      'analysis_protocol',params.analysis_protocol,...
      'settingsdir',params.settingsdir,...
      'datalocparamsfilestr',params.datalocparamsfilestr,...
      'hfig',handles.fig,...
      'period',params.date.period,...
      'maxdatenum',params.date.maxdatenum,...
      'figpos',get(handles.fig,'Position'),...
      'datenumnow',params.date.datenumnow,...
      'username',params.username,...
      'rootdatadir',params.rootdatadir,...
      'dataset',data.dataset,...
      'loadcacheddata',params.loadcacheddata,...
      'datafilename',state.datafilename,...
      'groupnames',params.groupnames,...
      'plotgroup',params.plotgroup,...
      'data_types',data.data_types);

  end

%% set manual pf callback
  
  function manualpf_popup_Callback(hObject,event) %#ok<INUSD>
    
    if isempty(state.datai_selected),
      warning('No experiment selected');
      return;
    end

    addToUndoList();
    
    s = get(hObject,'String');
    vcurr = get(hObject,'Value');
    manual_pf_full = s{vcurr};
    % first letter
    manual_pf_new = manual_pf_full(1);
    expdiris = find(data.groupidx == state.datai_selected);
    for expdiri = expdiris,
      rawdata(expdiri).(params.manual_fn) = manual_pf_new;
      % also store curator, curation time if experiment variables
      if strcmpi(params.vartype,'experiment'),
        rawdata(expdiri).manual_curator = params.username;
        rawdata(expdiri).manual_curation_date = datestr(now,constants.datetime_format);
      end
    end
    annot.manual_pf(state.datai_selected) = manual_pf_new;

    switch lower(manual_pf_new),
      
      case lower(params.pvalue),
        
        % update indices
        state.idx_manual_p(state.datai_selected) = true;
        state.idx_manual_f(state.datai_selected) = false;

        % update marker
        set(handles.data(state.datai_selected),'Marker','+');

        % update visible
        show_p = get(handles.menu.plot_manual_p,'Checked');
        state.idx_visible(state.datai_selected) = strcmpi(show_p,'on');
        set(handles.data(state.datai_selected),'Visible',show_p);

      case lower(params.fvalue),

        % update indices
        state.idx_manual_p(state.datai_selected) = false;
        state.idx_manual_f(state.datai_selected) = true;

        % update marker
        set(handles.data(state.datai_selected),'Marker','x');
        
        % update visible
        show_f = get(handles.menu.plot_manual_f,'Checked');
        state.idx_visible(state.datai_selected) = strcmpi(show_f,'on');
        set(handles.data(state.datai_selected),'Visible',show_f);
        
      case lower(params.uvalue),
        
        % update indices
        state.idx_manual_p(state.datai_selected) = false;
        state.idx_manual_f(state.datai_selected) = false;
        
        % update marker
        set(handles.data(state.datai_selected),'Marker','o');

        % update visible
        state.idx_visible(state.datai_selected) = true;
        set(handles.data(state.datai_selected),'Visible','on');
        
      otherwise
        
        error('Unknown manual_pf value %s',manual_pf_new);
        
    end
    state.needsave = true;

  end

%% notes callbacks
 function manualpf_pushbutton_Callback(hObject,event) %#ok<INUSD>
    
   handles.notes.dialog = dialog('Name','Add notes','WindowStyle','Normal','Resize','on');
   done_pos = [.29,.02,.2,.1];
   cancel_pos = [.51,.02,.2,.1];
   notes_pos = [.02,.14,.96,.84];

   if iscell(annot.notes_curation{state.datai_selected}),
     notes_curation_curr = annot.notes_curation{state.datai_selected};
   else
     notes_curation_curr = regexp(annot.notes_curation{state.datai_selected},'\\n','split');
   end

   handles.notes.pushbutton_done = uicontrol(handles.notes.dialog,'Style','pushbutton',...
     'Units','normalized','Position',done_pos,...
     'String','Done','Callback',@notes_done_Callback);
   handles.notes.pushbutton_cancel = uicontrol(handles.notes.dialog,'Style','pushbutton',...
     'Units','normalized','Position',cancel_pos,...
     'String','Cancel','Callback',@notes_cancel_Callback);
   handles.notes.edit_notes = uicontrol(handles.notes.dialog,'Style','edit',...
     'Units','normalized','Position',notes_pos,...
     'Min',0,'Max',25,...
     'String',notes_curation_curr,...
     'HorizontalAlignment','left',...
     'BackgroundColor','w');
   uiwait(handles.notes.dialog);
   
 end

  function notes_done_Callback(hObject,event) %#ok<INUSD>
    
    addToUndoList();

    annot.notes_curation{state.datai_selected} = get(handles.notes.edit_notes,'String');
    expdiris = find(data.groupidx == state.datai_selected);
    for expdiri = expdiris,
      rawdata(expdiri).notes_curation = annot.notes_curation{state.datai_selected};
    end
    
    % if notes_curation is one of the fields plotted, we will need to
    % update the plot data
    tmpi = find(strcmpi(data.statnames,'notes_curation'),1);
    if ~isempty(tmpi),
      s = annot.notes_curation{state.datai_selected};
      v = ~isempty(s) && ~strcmpi(s,'None') && ~strcmpi(s,'NULL');
      v = double(v)*2*3-3; % make -3, 3
      data.stat(state.datai_selected,tmpi) = v;
      data.normstat(state.datai_selected,tmpi) = ...
        (data.stat(state.datai_selected,tmpi) - data.mu(tmpi))/data.z(tmpi);
      set(handles.data(state.datai_selected),'YData',data.normstat(state.datai_selected,:));
      set(handles.selected,'YData',data.normstat(state.datai_selected,:));
      set(handles.selected1,'YData',data.normstat(data.datai_selected,state.stati_selected));
    end
    
    state.needsave = true;
    
    close(handles.notes.dialog);
    
  end

  function notes_cancel_Callback(hObject,event) %#ok<INUSD>
    
    close(handles.notes.dialog);
    
  end

%% info callback

 function manualpf_pushbutton_info_Callback(hObject,event) %#ok<INUSD>
   
   ncstats = 4;
   nrstats = ceil(data.nstats/ncstats);
   checkbox_h = 20;
   checkbox_w = 200;
   border = 10;
   donebutton_h = 20;
   donebutton_w = 60;
   dialog_h = nrstats*checkbox_h+border*3+donebutton_h;
   dialog_w = ncstats*checkbox_w+border*2;
   
   handles.info.dialog = dialog('Name','Add information','WindowStyle','Normal',...
     'Resize','on','Units','pixels','CloseRequestFcn',@info_done_Callback);
   SetFigureSize(handles.info.dialog,dialog_w,dialog_h);
   handles.info.checkboxes = nan(1,data.nstats);
   
   for tmpi = 1:data.nstats,
     [tmpr,tmpc] = ind2sub([nrstats,ncstats],tmpi);
     checkbox_pos = [border+(tmpc-1)*checkbox_w,dialog_h-(border+tmpr*checkbox_h),...
       checkbox_w,checkbox_h];
     handles.info.checkboxes(tmpi) = uicontrol(handles.info.dialog,'Style','checkbox',...
       'String',data.statnames{tmpi},'Position',checkbox_pos,...
       'Value',annot.info(state.datai_selected,tmpi));
   end
   handles.info.donebutton = uicontrol(handles.info.dialog,'Style','pushbutton',...
     'String','Done','Position',...
     [dialog_w/2-donebutton_w/2,border,donebutton_w,donebutton_h],...
     'Callback',@info_done_Callback);
   uiwait(handles.info.dialog);
   
 end

  function info_done_Callback(hObject,event) %#ok<INUSD>
    
    addToUndoList();
    
    for tmpi = 1:data.nstats,
      annot.info(state.datai_selected,tmpi) = get(handles.info.checkboxes(tmpi),'Value') == 1;
    end
    state.needsave = true;
    delete(handles.info.dialog);
    
  end

%% set rest callbacks

  function set_rest_manual_p_Callback(hObject,event) %#ok<INUSD>
    set_rest_manual_Callback(params.pvalue);
  end

  function set_rest_manual_f_Callback(hObject,event) %#ok<INUSD>
    set_rest_manual_Callback(params.fvalue);
  end

  function set_rest_manual_Callback(s)

    addToUndoList();
    
    b = questdlg(sprintf('Set all experiments plotted with %s = %s to %s?',params.manual_fn,params.uvalue,s),...
      sprintf('Set %s for rest?',params.manual_fn));
    if ~strcmpi(b,'Yes'),
      return;
    end
    
    idx = find(lower(annot.manual_pf) == lower(params.uvalue));
    
    if isempty(idx),
      return;
    end
    
    state.needsave = true;
    
    for expdiri = find(ismember(data.groupidx,idx)),
      rawdata(expdiri).(params.manual_fn) = s;
      if strcmpi(params.vartype,'experiment'),
        rawdata(expdiri).manual_curator = params.username;
        rawdata(expdiri).manual_curation_date = datestr(now,constants.datetime_format);
      end
    end
    annot.manual_pf(idx) = s;
    
    if lower(s) == lower(params.pvalue),

      % update indices
      state.idx_manual_p(idx) = true;
      state.idx_manual_f(idx) = false;
      
      % update marker
      set(handles.data(idx),'Marker','+');
      
      % update visible
      plot_p = get(handles.menu.plot_manual_p,'Checked');
      state.idx_visible(idx) = strcmpi(plot_p,'on');
      set(handles.data(idx),'Visible',plot_p);
      
    else

      % update indices
      state.idx_manual_p(idx) = false;
      state.idx_manual_f(idx) = true;
      
      % update marker
      set(handles.data(idx),'Marker','x');
      
      % update visible
      plot_f = get(handles.menu.plot_manual_f,'Checked');
      state.idx_visible(idx) = strcmpi(plot_f,'on');
      set(handles.data(idx),'Visible',plot_f);
      
    end
    
  end

%% save spreadsheet

  function save_Callback(hObject,event) %#ok<INUSD>
    
    while true
      [savefilename1,savepath1] = uiputfile(state.savefilename,sprintf('Save %s Curation tsv file',params.manual_fn));
      if ~ischar(savefilename1),
        return;
      end
      state.savefilename = fullfile(savepath1,savefilename1);
      rc.savepath = savepath1;
      %[savepath2,savefilename2] = fileparts(savefilename);
      %savefilename2 = fullfile(savepath2,[savefilename2,'_diagnosticinfo.mat']);
      
      fid = fopen(state.savefilename,'w');
      if fid < 0,
        warndlg(sprintf('Could not open file %s for writing. Make sure it is not open in another program.',state.savefilename),'Could not save');
        continue;
      end
      break;
    end
    if strcmpi(params.plotgroup,'none'),
      fprintf(fid,'#line\t');
    else
      fprintf(fid,'#%s\t',params.plotgroup);
    end
    if strcmpi(params.vartype,'experiment'),
      fprintf(fid,'experiment\t%s\tmanual_curator\tmanual_curation_date\tnotes_curation\tdiagnostic_fields\n',params.manual_fn);
    else
      fprintf(fid,'experiment\t%s\tnotes_curation\tdiagnostic_fields\n',params.manual_fn);
    end
    for tmpi = 1:data.ngroups,
      for tmpj = find(data.groupidx == tmpi),
        if iscell(rawdata(tmpj).notes_curation),
          notes_curation_curr = sprintf('%s\\n',rawdata(tmpj).notes_curation{:});
          % remove extra \n
          notes_curation_curr = notes_curation_curr(1:end-2);
        else
          notes_curation_curr = rawdata(tmpj).notes_curation;
        end
        % only \n?
        if numel(notes_curation_curr) >= 2 && strcmp(notes_curation_curr(end-1:end),'\n'),
          notes_curation_curr = notes_curation_curr(1:end-2);
        end
        % \n"?
        if numel(notes_curation_curr) >= 3 && strcmp(notes_curation_curr(end-2:end),'\n"'),
          notes_curation_curr = [notes_curation_curr(1:end-3),'"'];
        end
        if strcmpi(params.plotgroup,'none'),
          fprintf(fid,'%s\t',rawdata(tmpj).line_name);
        else
          fprintf(fid,'%s\t',data.groups{tmpi});
        end          
        fprintf(fid,'%s\t%s\t',...
          rawdata(tmpj).experiment_name,...
          rawdata(tmpj).(params.manual_fn));
        % print curation time, curator if experiment variables
        if strcmpi(params.vartype,'experiment'),
          fprintf(fid,'%s\t%s\t',...
            rawdata(tmpj).manual_curator,...
            rawdata(tmpj).manual_curation_date);
        end
        % also print info
        idx = annot.info(tmpi,:);
        if ~any(idx),
          tmps = '';
        else
          tmps = sprintf('%s,',data.statnames{idx});
          tmps = tmps(1:end-1);
        end
        fprintf(fid,'%s\t%s\n',notes_curation_curr,tmps);
      end
    end
    fclose(fid);
    state.needsave = false;
    
  end

%% load curation spreadsheet callback

  function load_Callback(hObject,event) %#ok<INUSD>
        
    if state.needsave,
      res = questdlg('Load state from file? All changes will be lost.');
      if ~strcmpi(res,'yes'),
        return;
      end
    end
    addToUndoList();
    
    [loadfilename,loadpath] = uigetfile(state.savefilename,'Load manual_behavior tsv');
    if ~ischar(loadfilename),
      return;
    end
    loadfilename = fullfile(loadpath,loadfilename);
    
    annotcurr = state.undolist(end);
    
    fid = fopen(loadfilename,'r');
    groupisset = false(1,data.ngroups);
    experiment_names = {rawdata.experiment_name};
    while true,
      
      ss = fgetl(fid);
      % end of file
      if ~ischar(ss), break; end
      
      % comments
      if isempty(ss) || ~isempty(regexp(ss,'^\s*$','once')) || ...
          ~isempty(regexp(ss,'^\s*#','once')),
        continue;
      end
      
      % split at tabs
      m = regexp(ss,'\t','split');
%       % no group name
%       if strcmpi(params.plotgroup,'none'),
%         m = [{''},m];
%       end
      if numel(m) < 5 || numel(m) > 7,
        warning('Skipping line %s: wrong number of fields',ss);
        continue;
      end
      if numel(m) < 7,
        m = [m,repmat({''},[1,7-numel(m)])]; %#ok<AGROW>
      end

      % remove ""s
      for tmpi = 1:numel(m),
        if regexp(m{tmpi},'^".*"$','once'),
          m{tmpi} = m{tmpi}(2:end-1);
        end
      end
      
      %groupname = m{1};
      experiment_name = m{2};
      manual_pf_curr = m{3};
      manual_curator = m{4};
      manual_curation_date = m{5};
      notes_curation_curr = m{6};
      notes_curation_curr = regexp(notes_curation_curr,'\\n','split');
      info_s = regexp(m{7},',','split');
      info_s = setdiff(info_s,{''});
      info_curr = ismember(data.statnames,info_s);
      unknown_stats = ~ismember(info_s,data.statnames);
      if any(unknown_stats),
        warning(['Unknown info stats: ',sprintf('%s ',info_s{unknown_stats})]);
      end
      
      tmpi = find(strcmp(experiment_names,experiment_name),1);
      expi = tmpi;
      if isempty(tmpi),
        fprintf('experiment %s not currently examined, skipping\n',experiment_name);
        continue;
      end
      
      if ~strcmpi(params.plotgroup,'none'),
        tmpi = data.groupidx(tmpi);
      end
      groupname = data.groups{tmpi};
      
      if groupisset(tmpi) && lower(annotcurr.manual_pf(tmpi)) ~= lower(params.uvalue),
        if lower(manual_pf_curr(1)) == lower(params.uvalue),
          continue;
        elseif manual_pf_curr(1) ~= annotcurr.manual_pf(tmpi),
          warning('Conflicting labels for %s %s, using loaded value %s by %s at %s',params.plotgroup,groupname,annotcurr.manual_pf(tmpi),manual_curator,manual_curation_date);
        end
      end
      groupisset(tmpi) = true;
      
      annotcurr.manual_pf(tmpi) = manual_pf_curr(1);
      annotcurr.notes_curation{tmpi} = notes_curation_curr;
      annotcurr.info(tmpi,:) = info_curr;
      rawdata(expi).manual_curator = manual_curator;
      rawdata(expi).manual_curation_date = manual_curation_date;
      rawdata(expi).(params.manual_fn) = manual_pf_curr(1);
      rawdata(expi).notes_curation = notes_curation_curr;
    end
    fclose(fid);

    state.needsave = true;
    setManualPFState(annotcurr);
    
  end

%% search callback

  function search_Callback(hObject,event) %#ok<INUSD>
    
    if ~isfield(handles.search,'dialog') || ~ishandle(handles.search.dialog),
      line_height = 20;
      prompt_width = 50;
      answer_width = 200;
      button_width = 75;
      button_height = 30;
      gap_width = 5;
      gap_height = 10;
      nl = 3;
      dlg_height = (nl-1)*line_height + button_height + (nl+1)*gap_height;
      dlg_width = prompt_width + answer_width + 3*gap_width;
      handles.search.dialog = dialog('Name','Find experiments','WindowStyle','normal',...
        'Units','pixels');
      SetFigureSize(handles.search.dialog,dlg_width,dlg_height);
      handles.search.field_text = uicontrol(handles.search.dialog,'Style','text',...
        'Units','pixels',...
        'Position',[gap_width,dlg_height-(line_height+gap_height),prompt_width,line_height],...
        'String','Field: ',...
        'HorizontalAlignment','right');
      groupi = find(strcmp(params.groupnames,'line_name'),1);
      if isempty(groupi),
        groupi = 1;
      end
      handles.search.field_popupmenu = uicontrol(handles.search.dialog,'Style','popupmenu',...
        'Units','pixels',...
        'Position',[prompt_width+2*gap_width,dlg_height-(line_height+gap_height),answer_width,line_height],...
        'String',params.groupnames,...
        'Value',groupi,...
        'HorizontalAlignment','right',...
        'Callback',@search_field_Callback);
      handles.search.value_text = uicontrol(handles.search.dialog,'Style','text',...
        'Units','pixels',...
        'Position',[gap_width,dlg_height-2*(line_height+gap_height),prompt_width,line_height],...
        'String','Value: ',...
        'HorizontalAlignment','right');
      if iscell(params.groupvalues{groupi}),
        group_s = params.groupvalues{groupi};
      else
        group_s = num2str(params.groupvalues{groupi}(:));
      end
      handles.search.value_popupmenu = uicontrol(handles.search.dialog,'Style','popupmenu',...
        'Units','pixels',...
        'Position',[prompt_width+2*gap_width,dlg_height-2*(line_height+gap_height),answer_width,line_height],...
        'String',group_s,...
        'Value',1,...
        'HorizontalAlignment','right',...
        'Callback',@search_value_Callback);
      handles.search.findnext_pushbutton = uicontrol(handles.search.dialog,'Style','pushbutton',...
        'Units','pixels',...
        'Position',[dlg_width/2-button_width/2,dlg_height-(2*line_height+button_height+3*gap_height),button_width,button_height],...
        'String','Find next',...
        'HorizontalAlignment','center',...
        'Callback',@search_findnext_Callback);
      
    else
      figure(handles.search.dialog);
    end
    
    handles.search.resulti = 0;
    
  end

  function search_field_Callback(varargin)
    
    groupi = get(handles.search.field_popupmenu,'Value');
    if iscell(params.groupvalues{groupi}),
      group_s = params.groupvalues{groupi};
    else
      group_s = num2str(params.groupvalues{groupi}(:));
    end
    set(handles.search.value_popupmenu,'String',group_s,'Value',1);
    handles.search.resulti = 0;
    
  end

  function search_value_Callback(varargin)
    
    handles.search.resulti = 0;
    
  end

  function search_findnext_Callback(varargin)
    
    groupi = get(handles.search.field_popupmenu,'Value');
    groupname = params.groupnames{groupi};
    groupvaluei = get(handles.search.value_popupmenu,'Value');
    % need to do the search
    if handles.search.resulti == 0,
      if iscell(params.groupvalues{groupi}),
        results_exp = strcmp({rawdata.(groupname)},params.groupvalues{groupi}{groupvaluei});
      else
        results_exp = find([rawdata.(groupname)] == params.groupvalues{groupi}(groupvaluei));
      end
      handles.search.results = unique(data.groupidx(results_exp));
      if isempty(handles.search.results),
        uiwait(warndlg(sprintf('No experiments found with %s = %s',groupname,params.groupvalues{groupi}{groupvaluei})));
        return;
      end
    end
    
    handles.search.resulti = handles.search.resulti + 1;
    if handles.search.resulti > numel(handles.search.results),
      handles.search.resulti = 1;
    end
    
    state.datai_selected = handles.search.results(handles.search.resulti);
    fprintf('Selecting group %s\n',data.groups{state.datai_selected});
    
    UpdateSelected();
    
  end

%% reset state

  function setManualPFState(s)
  
 
    annot.manual_pf = s.manual_pf;
    annot.notes_curation = s.notes_curation;
    for expdiri = 1:data.nexpdirs,
      rawdata(expdiri).(params.manual_fn) = annot.manual_pf(data.groupidx(expdiri));
      rawdata(expdiri).notes_curation = annot.notes_curation{data.groupidx(expdiri)};
    end
    state.idx_manual_p = lower(annot.manual_pf) == lower(params.pvalue);
    state.idx_manual_f = lower(annot.manual_pf) == lower(params.fvalue);
    state.idx_visible = true(1,data.ngroups);
    if strcmpi(get(handles.menu.plot_manual_p,'Checked'),'off'),
      state.idx_visible(state.idx_manual_p) = false;
    end
    if strcmpi(get(handles.menu.plot_manual_f,'Checked'),'off'),
      state.idx_visible(state.idx_manual_f|(state.idx_automated_f&~state.idx_manual_p)) = false;
    end
    set(handles.data,'Marker','o');
    set(handles.data(state.idx_manual_p),'Marker','+');
    set(handles.data(state.idx_manual_f|(state.idx_automated_f&~state.idx_manual_p)),'Marker','x');
    set(handles.data(state.idx_visible),'Visible','on');
    set(handles.data(~state.idx_visible),'Visible','off');
    
    if ~isempty(state.datai_selected),
      manual_pf_curr = annot.manual_pf(state.datai_selected);
      ss = get(handles.manualpf.popup,'String');
      vcurr = find(strncmpi(manual_pf_curr,ss,1),1);
      if isempty(vcurr),
        error('Unknown %s %s',params.manual_fn,manual_pf_curr);
      end
      set(handles.manualpf.popup,'Value',vcurr);
    end
    
    annot.info = s.info;
    
    
    % if notes_curation is one of the fields plotted, we will need to
    % update the plot data
    % TODO    
    
  end

%% undo callback

  function undo_Callback(hObject,event) %#ok<INUSD>
    
    if isempty(state.undolist), return; end
    setManualPFState(state.undolist(end));
    state.undolist = state.undolist(1:end-1);
    if isempty(state.undolist),
      set(handles.menu.undo,'Enable','off');
    end
    
  end

%% add state to undo list

  function addToUndoList()

    if numel(state.undolist) >= constants.maxnundo,
      state.undolist = state.undolist(end-constants.maxnundo+1:end);
    end
    state.undolist = structappend(state.undolist,annot);
    set(handles.menu.undo,'Enable','on');
    
  end

%% save data callback

  function savedata_Callback(hObject,event) %#ok<INUSD>

    while true
      [savefilename1,savepath1] = uiputfile(state.datafilename,'Save cached data');
      if ~ischar(savefilename1),
        return;
      end
      state.datafilename = fullfile(savepath1,savefilename1);
      rc.savedatapath = savepath1;
      try
        save(state.datafilename,'-struct','data');
        save('-append',state.datafilename,'rawdata');
        save('-append',state.datafilename,'params');
      catch ME
        warndlg(sprintf('Could not save to file %s: %s',state.datafilename,getReport(ME)),'Could not save');
        continue;
      end
      break;
    end
 
  end

  function loaddata_Callback(varargin)
    
    if state.needsave,
      res = questdlg('Load new data? All changes will be lost.');
      if ~strcmpi(res,'yes'),
        return;
      end
    end
    
    [loadfilename1,loadpath1] = uigetfile(state.datafilename,'Load cached data');
    if ~ischar(loadfilename1),
      return;
    end
    state.datafilename = fullfile(loadpath1,loadfilename1);
    rc.savedatapath = loadpath1;
    params.loadcacheddata = true;
    params.datafilename = state.datafilename;
    RecursiveCall();
    
  end

end

function [handles,data] = FlyBowlExamineBehaviorVariables(varargin)

handles = struct;
data = [];
groupnames = SetDefaultGroupNames();

[analysis_protocol,settingsdir,datalocparamsfilestr,hfig,period,maxdatenum,...
  figpos,datenumnow,sage_params_path,sage_db,username,rootdatadir,loadcacheddata,leftovers] = ...
  myparse_nocheck(varargin,...
  'analysis_protocol','current',...
  'settingsdir','/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings',...
  'datalocparamsfilestr','dataloc_params.txt',...
  'hfig',1,...
  'period',7,...
  'maxdatenum',[],...
  'figpos',[],...
  'datenumnow',[],...
  'sage_params_path','',...
  'sage_db',[],...
  'username','',...
  'rootdatadir','/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data',...
  'loadcacheddata','',...
  'groupnames',groupnames);

maxnundo = 5;
maxperiodsprev = 10;
tempfilescreated = {};

didinputweek = ~isempty(maxdatenum) && ~isempty(datenumnow);
if ~didinputweek,
  datenumnow = now;
end

if isempty(maxdatenum),
  maxdatenum = datenumnow;
  daycurr = weekday(maxdatenum);
  daywant = 7;
  maxdatenum = floor(maxdatenum)+daywant-daycurr-7;

end
%format = 'yyyy-mm-ddTHH:MM:SS';
format = 'yyyymmddTHHMMSS';
mindatenum = maxdatenum - period;
mindatestr = datestr(mindatenum,format);
maxdatestr = datestr(maxdatenum,format);
daterange = {mindatestr,maxdatestr};
maxdaysprev = datenumnow-mindatenum;
mindaysprev = datenumnow-maxdatenum;

% first day of week choices
mindatenum_choices = fliplr([mindatenum-period*maxperiodsprev:period:mindatenum-period,mindatenum:period:datenumnow]);
maxdatenum_choices = mindatenum_choices + period;
mindatestr_choices = datestr(mindatenum_choices,'yyyy-mm-dd');
maxdatestr_choices = datestr(maxdatenum_choices,'yyyy-mm-dd');

%% read parameters

datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
dataloc_params = ReadParams(datalocparamsfile);
examineparamsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.examinebehaviorvariablesparamsfilestr);
examine_params = ReadParams(examineparamsfile);
examinefile = fullfile(settingsdir,analysis_protocol,dataloc_params.examinebehaviorvariablesfilestr);
examinestats = ReadExamineBehaviorVariables(examinefile);

if exist(examine_params.rcfile,'file'),
  rc = load(examine_params.rcfile);
else
  rc = struct('username','','savepath','');
end

%% get user name

havequestions = isempty(username) || ~didinputweek;

if havequestions,
  
  if ~isempty(username),
    rc.username = username;
  end

  if didinputweek,
    daterange_strings = {sprintf('%s to %s',mindatestr,maxdatestr)};
  else
    daterange_strings = cellstr(cat(2,mindatestr_choices,...
      repmat(' to ',[numel(maxdatenum_choices),1]),...
      maxdatestr_choices));
  end
  weekidx = find(maxdatenum_choices == maxdatenum,1);
  [success,newweekidx,newusername] = ExamineInitializeParams(daterange_strings,weekidx,rc.username);
  if ~success,
    return;
  end
  if ~isempty(newusername),
    rc.username = newusername;
  end
  if ~didinputweek,
    maxdatenum = maxdatenum_choices(newweekidx);
    mindatenum = mindatenum_choices(newweekidx);
    mindatestr = datestr(mindatenum,format);
    maxdatestr = datestr(maxdatenum,format);
    daterange = {mindatestr,maxdatestr};
    maxdaysprev = datenumnow-mindatenum;
    mindaysprev = datenumnow-maxdatenum;
  end
  
end

if ~isfield(rc,'savepath'),
  rc.savepath = '';
end
savefilename = fullfile(rc.savepath,sprintf('ManualBehaviorAnnotation_FlyBowl_%s_%sto%s.tsv',...
  rc.username,datestr(mindatenum,'yyyymmdd'),datestr(maxdatenum,'yyyymmdd')));
needsave = false;

%% string metadata

flag_metadata_fns = {'flag_redo','flag_review'};
note_metadata_fns = {'notes_behavioral','notes_technical','notes_curation'};
string_metadata_fns = {'line_name','manual_pf','automated_pf','experiment_name','experiment_protocol',...
  'experimenter','exp_datetime','bowl','camera','computer','harddrive','apparatus_id','cross_date','effector',...
  'environmental_chamber','gender','genotype','handler_sorting','handler_starvation','plate','rearing_protocol',...
  'rig','top_plate','file_system_path','manual_behavior'};
flag_options = {{'None',''},{'Rearing problem','Flies look sick',...
  'See behavioral notes','See technical notes'}};

%% connect to SAGE for storing manualnd
% 
% if isempty(sage_db) && ~isempty(sage_params_path),
%   try
%     sage_db = ConnectToSage(sage_params_path);
%   catch ME,
%     getReport(ME);
%     warning('Could not connect to Sage');
%     sage_db = [];
%   end
% else
%   sage_db = [];
% end

%% get data
%data_types = {'stats_perframe_*'};
%data_types = examinestats;

didloaddata = false;
if ~isempty(loadcacheddata),
  try
    load(loadcacheddata,'data');
    didloaddata = true;
  catch ME
    warning('Could not load data from %s:\n%s',loadcacheddata,getReport(ME));
  end
end

if ~didloaddata,
  
  queries = leftovers;
  queries(end+1:end+2) = {'daterange',daterange};
  queries(end+1:end+2) = {'data_type','stats_perframe_*'};
  queries(end+1:end+2) = {'flag_aborted',0};
  queries(end+1:end+2) = {'automated_pf','P'};
  queries(end+1:end+2) = {'experiment_name','FlyBowl_*'};
  data = SAGEGetBowlData(queries{:},'removemissingdata',true,'rootdir',rootdatadir);
  if isempty(data),
    uiwait(warndlg(sprintf('No data for date range %s to %s',daterange{:}),'No data found'));
    if didinputweek,
      fprintf('Trying previous week...\n');
      maxdatenum = maxdatenum - period;
      handles = FlyBowlExamineBehaviorVariables(...
        'analysis_protocol',analysis_protocol,...
        'settingsdir',settingsdir,...
        'datalocparamsfilestr',datalocparamsfilestr,...
        'hfig',hfig,...
        'figpos',figpos,...
        'period',period,...
        'maxdatenum',maxdatenum,...
        'datenumnow',datenumnow,...
        'sage_params_path',sage_params_path,...
        'sage_db',sage_db,...
        'username',rc.username,...
        'rootdatadir',rootdatadir,...
        'loadcacheddata',loadcacheddata,...
        'groupnames',groupnames);      
      return;
    else
      fprintf('Choose a different week.\n');
      handles = FlyBowlExamineBehaviorVariables(...
        'analysis_protocol',analysis_protocol,...
        'settingsdir',settingsdir,...
        'datalocparamsfilestr',datalocparamsfilestr,...
        'hfig',hfig,...
        'figpos',figpos,...
        'period',period,...
        'maxdatenum',[],...
        'datenumnow',datenumnow,...
        'sage_params_path',sage_params_path,...
        'sage_db',sage_db,...
        'username',rc.username,...
        'rootdatadir',rootdatadir,...
        'loadcacheddata',loadcacheddata,...
        'groupnames',groupnames);
      return;
    end
  end
  
end
% sort by date
date = {data.exp_datetime};
[~,order] = sort(date);
data = data(order);
nexpdirs = numel(data);
expdir_bases = {data.experiment_name};
expdir_bases = cellfun(@(s) regexprep(s,'^FlyBowl_',''),expdir_bases,'UniformOutput',false);

%% for now, add on extra diagnostics here. we will store these later
need_file_system_path = ~isfield(data,'file_system_path');
need_manual_behavior = ~isfield(data,'manual_behavior');
need_notes_curation = ~isfield(data,'notes_curation');
for i = 1:nexpdirs,
  if need_file_system_path,
    data(i).file_system_path = fullfile(rootdatadir,expdir_bases{i});
  end
  if need_manual_behavior,
    data(i).manual_behavior = 'U';
  end
  if need_notes_curation,
    data(i).notes_curation = '';
  end
end
expdirs = {data.file_system_path};

%% groups of variables for searching

groupnames = intersect(groupnames,fieldnames(data));
groupvalues = cell(size(groupnames));
for i = 1:numel(groupnames),
  fn = groupnames{i};
  values = {data.(fn)};
  if all(cellfun(@ischar,values)),
    values = unique(values);
  else
    values = unique([values{:}]);
  end
  groupvalues{i} = values;
end

%% get behavior variables

nstats = numel(examinestats);
% average over lines

[linenames,~,lineidx] = unique({data.line_name});
nlines = numel(linenames);
stat = nan(nlines,nstats);
expstat = nan(nexpdirs,nstats);

for i = 1:nstats,

  if ~isfield(data,examinestats{i}{1}),
    continue;
  end
  
  % special case: flags
  if ismember(examinestats{i}{1},flag_metadata_fns),
    % flags are now binary
    % make -3, 3
    v = double([data.(examinestats{i}{1})])*2*3-3;
    for j = 1:nlines,
      stat(j,i) = nanmean(v(lineidx==j));
    end
    expstat(:,i) = v;
    
%   if ismember(examinestats{i}{1},flag_metadata_fns),
%     
%     % which set of flags
%     v = zeros(1,nexpdirs);
%     for j = 1:numel(flag_options),
%       v(ismember({data.(examinestats{i}{1})},flag_options{j})) = j;
%     end
%     v = (v/numel(flag_options)*2-1)*3;
%     stat(:,i) = v;

    
  % special case: notes
  elseif ismember(examinestats{i}{1},note_metadata_fns),
    
    v = cellfun(@(s) ~isempty(s) && ~strcmpi(s,'None'),{data.(examinestats{i}{1})});
    v = double(v)*2*3-3; % make -3, 3
    % take average over lines
    for j = 1:nlines,
      stat(j,i) = nanmean(v(lineidx==j));
    end
    expstat(:,i) = v;

  % special case: strings
  elseif ismember(examinestats{i}{1},string_metadata_fns),
    [uniquevals,~,v] = unique({data.(examinestats{i}{1})});
    if numel(uniquevals) == 1,
      v(:) = 0;
    else
      v = v - 1;
      v = (v/max(v)*2-1)*3;
    end
    for j = 1:nlines,
      stat(j,i) = nanmean(v(lineidx==j));
    end
    expstat(:,i) = v;

  % numbers
  else
    
    datacurr = {data.(examinestats{i}{1})};
    for k = 2:numel(examinestats{i}),
      datacurr = cellfun(@(s) s.(examinestats{i}{k}),datacurr,'UniformOutput',false);
    end
    badidx = cellfun(@isempty,datacurr);
    for j = find(badidx),
      datacurr{j} = nan;
      warning('No data for stat %s experiment %s',sprintf('%s,',examinestats{i}{:}),strtrim(data(j).experiment_name));
    end
    datacurr = cell2mat(datacurr);
    for j = 1:nlines,
      stat(j,i) = nanmean(datacurr(lineidx==j));
    end
    expstat(:,i) = datacurr;
    
  end
end

%% get mean and standard deviation for z-scoring

if isempty(examine_params.examinebehaviorvariablesstatsfilestr),
  mu = nanmean(stat,1);
  sig = nanstd(stat,1,1);  
else
  examinebehaviorvariablesstatsfile = fullfile(settingsdir,analysis_protocol,...
    examine_params.examinebehaviorvariablesstatsfilestr);
  normstats = load(examinebehaviorvariablesstatsfile,'meanstatsperfly');
  normstats = normstats.meanstatsperfly;
  mu = nan(1,nstats);
  sig = nan(1,nstats);
  for i = 1:nstats,
    canload = true;
    if isempty(regexp(examinestats{i}{1},'^stats_perframe','once')),
      canload = false;
    end
    if canload,
      if numel(examinestats{i}) ~= 3,
        canload = false;
      end
    end
    if canload,
      fn = sprintf('%s_%s',examinestats{i}{1},examinestats{i}{3});
      fn = regexprep(fn,'^stats_perframe_','');
      canload = isfield(normstats,fn);
    end
    if canload,
      mu(i) = normstats.(fn).meanmean;
      sig(i) = normstats.(fn).stdmean;
    else
      mu(i) = 0;
      sig(i) = 1;
%       fprintf('Not normalizing');
%       fprintf(' %s',examinestats{i}{:});
%       fprintf('\n');
    end
  end
end

%% abbr names of stats

statnames = SetStatNames(examinestats);

%% z-score stats

z = sig;
z(z == 0) = 1;
normstat = bsxfun(@rdivide,bsxfun(@minus,stat,mu),z);

%% create figure

if ishandle(hfig),
  close(hfig);
end
figure(hfig);
clf(hfig,'reset');
if isempty(figpos),
  figpos = examine_params.figpos;
end
set(hfig,'Units','Pixels','Position',figpos,'MenuBar','none','ToolBar','figure');
hax = axes('Parent',hfig,'Units','Normalized','Position',examine_params.axespos);
% plot 0
plot(hax,[0,nstats+1],[0,0],'k-','HitTest','off');
hold(hax,'on');
colors = jet(nlines)*.7;
drawnow;

%% plot diagnostics

if nlines == 1,
  off = 0;
else
  off = (2*((1:nlines)-(nlines+1)/2)/(nlines-1))*examine_params.offx;
end

h = nan(1,nlines);
x = nan(nlines,nstats);
for linei = 1:nlines,
  x(linei,:) = (1:nstats)+off(linei);
  h(linei) = plot(hax,x(linei,:),normstat(linei,:),'o',...
    'color',colors(linei,:),'markerfacecolor',colors(linei,:),...
    'markersize',6,'HitTest','off');
end

xlim = [0,nstats+1];
miny = min(normstat(:));
maxy = max(normstat(:));
dy = maxy - miny;
if dy == 0,
  maxy = miny + .001;
end
ylim = [miny-.01*dy,maxy+.01*dy];
dx = diff(xlim)/figpos(3);
dy = diff(ylim)/figpos(4);
set(hax,'XLim',xlim,'YLim',ylim,'XTick',1:nstats,'XTickLabel',statnames,'XGrid','on');

ylabel(hax,'Stds from mean');

%% set manual_behavior = n, d marker

%gray_p_color = [.7,.7,.7];
%gray_f_color = [.5,.5,.5];
badidx = cellfun(@isempty,{data.manual_behavior});
for j = find(badidx),
  warning('No data for manual_behavior, experiment %s',strtrim(data(j).experiment_name));
end
manual_behavior = repmat('U',[1,nlines]);
for i = 1:nlines,
  s1 = lower([data(lineidx==i).manual_behavior]);
  if ismember('n',s1) && ismember('d',s1),
    warning('Experiments for line %s are marked as both normal (%d/%d) and different (%d/%d). Using unknown.',...
      linenames(i),nnz(s1=='n'),numel(s1),nnz(s1=='d'),numel(s1));
    manual_behavior(i) = 'U';
  elseif ismember('n',s1),
    manual_behavior(i) = 'N';
  elseif ismember('d',s1'),
    manual_behavior(i) = 'D';
  else
    manual_behavior(i) = 'U';
  end
end
idx_manual_n = lower(manual_behavior) == 'n';
idx_manual_d = lower(manual_behavior) == 'd';
idx_visible = true(1,nlines);
set(h(idx_manual_d),'Marker','+');
set(h(idx_manual_n),'Marker','x');

notes_curation = cell(1,nlines);
for i = 1:nlines,
  notes_curation{i} = unique({data(lineidx==i).notes_curation});
end

%% selected experiment

hselected = plot(0,0,'o','color','k','Visible','off','HitTest','off','MarkerSize',10,'MarkerFaceColor','k');
hselected1 = plot(0,0,'o','color','r','Visible','off','HitTest','off','MarkerSize',10,'MarkerFaceColor','r');
linei_selected = [];
stati_selected = [];


%% Examine menu

hmenu = struct;
hmenu.file = uimenu('Label','File','Parent',hfig);
hmenu.options = uimenu('Label','Options','Parent',hfig);
hmenu.set = uimenu('Label','Set','Parent',hfig);
hmenu.info = uimenu('Label','Info','Parent',hfig);
hmenu.plot_manual_n = uimenu(hmenu.options,'Label','Plot manual_behavior = normal',...
  'Checked','on','Callback',@plot_manual_n_Callback);
hmenu.plot_manual_d = uimenu(hmenu.options,'Label','Plot manual_behavior = different',...
  'Checked','on','Callback',@plot_manual_d_Callback);
hmenu.set_rest_manual_n = uimenu(hmenu.set,'Label','Set manual_behavior == unknown -> normal',...
  'Callback',@set_rest_manual_n_Callback);
hmenu.set_rest_manual_d = uimenu(hmenu.set,'Label','Set manual_behavior == unknown -> different',...
  'Callback',@set_rest_manual_d_Callback);
hmenu.undo = uimenu(hmenu.set,'Label','Undo',...
  'Callback',@undo_Callback,'Enable','off','Accelerator','z');
hmenu.save = uimenu(hmenu.file,'Label','Save...',...
  'Callback',@save_Callback,'Accelerator','s');
hmenu.load = uimenu(hmenu.file,'Label','Load Spreadsheet...',...
  'Callback',@load_Callback,'Accelerator','l');
hmenu.open = uimenu(hmenu.file,'Label','View Experiment',...
  'Callback',@open_Callback,'Enable','off','Accelerator','o');
hmenu.search = uimenu(hmenu.info,'Label','Search...',...
  'Callback',@search_Callback,'Accelerator','f');
hsearch = struct;

daterangeprint = {datestr(mindatenum,'yyyy-mm-dd'),datestr(maxdatenum,'yyyy-mm-dd')};
hmenu.daterange = uimenu(hmenu.info,'Label',...
  sprintf('Date range: %s - %s ...',daterangeprint{:}));
hmenu.username = uimenu(hmenu.info,'Label',...
  sprintf('User: %s',rc.username));
set(hfig,'CloseRequestFcn',@close_fig_Callback);

%% text box

set(hax,'Units','normalized');
axpos = get(hax,'Position');
textpos = [axpos(1),axpos(2)+axpos(4),axpos(3)/2,1-axpos(2)-axpos(4)-.01];
htext = annotation('textbox',textpos,'BackgroundColor','k','Color','g',...
  'String','Experiment info','Interpreter','none');

%% manual nd

set(htext,'Units','Pixels');
textpos_px = get(htext,'Position');
set(htext,'Units','normalized');
margin = 5;
w1 = 80;
w2 = 100;
w3 = 100;
h1 = 20;
h2 = 30;
c1 = ((figpos(4)-margin) + textpos_px(2))/2;
manual_behavior_textpos = [textpos_px(1)+textpos_px(3)+margin,c1-h1/2,w1,h1];
hmanual_behavior = struct;
hmanual_behavior.text = uicontrol(hfig,'Style','text','Units','Pixels',...
  'Position',manual_behavior_textpos,'String','Behavior:',...
  'BackgroundColor',get(hfig,'Color'),'Visible','off');
manual_behavior_popuppos = [manual_behavior_textpos(1)+manual_behavior_textpos(3)+margin,...
  c1-h2/2,w2,h2];
hmanual_behavior.popup = uicontrol(hfig,'Style','popupmenu','Units','Pixels',...
  'Position',manual_behavior_popuppos,'String',{'Normal','Different','Unknown'},...
  'Value',3,'Visible','off','Callback',@manual_behavior_popup_Callback);
manual_behavior_pushbuttonpos = [manual_behavior_popuppos(1)+manual_behavior_popuppos(3)+margin,...
  c1-h2/2,w3,h2];
hmanual_behavior.pushbutton = uicontrol(hfig,'Style','pushbutton','Units','Pixels',...
  'Position',manual_behavior_pushbuttonpos,'String','Add Note...',...
  'Callback',@manual_behavior_pushbutton_Callback,...
  'Visible','off');
manual_behavior_pushbutton_info_pos = [manual_behavior_pushbuttonpos(1)+manual_behavior_pushbuttonpos(3)+margin,...
  c1-h2/2,w3,h2];
hmanual_behavior.pushbutton_info = uicontrol(hfig,'Style','pushbutton','Units','Pixels',...
  'Position',manual_behavior_pushbutton_info_pos,'String','Add Info...',...
  'Callback',@manual_behavior_pushbutton_info_Callback,...
  'Visible','off');
set([hmanual_behavior.text,hmanual_behavior.popup,hmanual_behavior.pushbutton,hmanual_behavior.pushbutton_info],'Units','normalized');

hnotes = struct;
hinfo = struct;
info = false(nlines,nstats);


%% date

set(hax,'Units','Pixels');
axpos_px = get(hax,'Position');
set(hax,'Units','normalized');
w1 = 20;
h1 = 20;
w2 = 150;
margin = 5;
nextpos = [axpos_px(1)+axpos_px(3)-w1,axpos_px(2)+axpos_px(4)+margin,w1,h1];
currpos = [nextpos(1)-margin-w2,nextpos(2),w2,h1];
prevpos = [currpos(1)-margin-w1,nextpos(2),w1,h1];
hdate = struct;
if mindaysprev <= 0,
  s = 'this week';
elseif mindaysprev < 7,
  s = 'last week';
else
  s = sprintf('%d weeks ago',ceil(mindaysprev/7));
end
hdate.curr = uicontrol(hfig,'Style','text','Units','Pixels',...
  'Position',currpos,'String',s);
hdate.next = uicontrol(hfig,'Style','pushbutton','Units','Pixels',...
  'Position',nextpos,'String','>','Callback',@nextdate_Callback);
hdate.prev = uicontrol(hfig,'Style','pushbutton','Units','Pixels',...
  'Position',prevpos,'String','<','Callback',@prevdate_Callback);
if mindaysprev < .0001,
  set(hdate.next,'Enable','off');
else
  set(hdate.prev,'Enable','on');
end
set([hdate.curr,hdate.prev,hdate.next],'Units','normalized');

%% rotate x-tick labels

hx = rotateticklabel(hax,90);
% make sure the ticks don't overlap the x-axis
ex = get(hx(1),'Extent');
y1 = ex(2)+ex(4);
offy = y1-ylim(1);
for i = 1:numel(hx),
  pos = get(hx(i),'Position');
  pos(2) = pos(2) - offy;
  set(hx(i),'Position',pos);
end
set(hx,'Interpreter','none');

%% set buttondownfcn

set(hax,'ButtonDownFcn',@ButtonDownFcn);

% hchil = findobj(hfig,'Units','Pixels');
% set(hchil,'Units','normalized');

%% set motion function

set(hfig,'WindowButtonMotionFcn',@MotionFcn,'BusyAction','queue');

%% undo list

undolist = [];

%% return values

handles = struct;
handles.hx = hx;
handles.hax = hax;
handles.hfig = hfig;
handles.h = h;
handles.hselected = hselected;
handles.hselected1 = hselected1;
handles.htext = htext;
handles.hdate = hdate;
% handles.hcmenu_manual_behavior = hcmenu_manual_behavior;
% handles.hcmenu = hcmenu;

%% set stat names

  function statnames = SetStatNames(examinestats)

    statnames = cell(1,nstats);
    
    for tmpi = 1:nstats,
      
      statnames{tmpi} = sprintf('%s_',examinestats{tmpi}{:});
      statnames{tmpi} = statnames{tmpi}(1:end-1);
      statnames{tmpi} = strrep(statnames{tmpi},'stats_perframe_','');
      statnames{tmpi} = strrep(statnames{tmpi},'meanmean_perexp_','');
      statnames{tmpi} = strrep(statnames{tmpi},'flyany_frame','');
      statnames{tmpi} = strrep(statnames{tmpi},'ufmf_diagnostics_summary','ufmf');
      statnames{tmpi} = strrep(statnames{tmpi},'ufmf_diagnostics_stream','ufmf');
      statnames{tmpi} = strrep(statnames{tmpi},'temperature_diagnostics','temp');
      statnames{tmpi} = strrep(statnames{tmpi},'bias_diagnostics','bias');
      statnames{tmpi} = strrep(statnames{tmpi},'bkgd_diagnostics','bkgd');
      statnames{tmpi} = strrep(statnames{tmpi},'ctrax_diagnostics','ctrax');
      statnames{tmpi} = strrep(statnames{tmpi},'registrationdata','reg');
      statnames{tmpi} = strrep(statnames{tmpi},'sexclassifier_diagnostics','sex');
      statnames{tmpi} = strrep(statnames{tmpi},'stats_perframe_','');
      statnames{tmpi} = strrep(statnames{tmpi},'flyany_frame','');
      statnames{tmpi} = strrep(statnames{tmpi},'_perexp','');
      
    end
    
  end

%% print function

  function s = printfun(linei,stati)
    
    expdiris = find(lineidx == linei);
    s = cell(1,2+numel(expdiris));
    s{1} = linenames{linei};
    % special case: notes, strings
    if isempty(stati),
      s = s(1);
    else
      if ismember(examinestats{stati}{1},note_metadata_fns) || ...
          ismember(examinestats{stati}{1},string_metadata_fns),
        s{2} = sprintf('%s =',statnames{stati});
        for expdirii = 1:numel(expdiris),
          expdiri = expdiris(expdirii);
          if iscell(data(expdiri).(examinestats{stati}{1})),
            s{2+expdirii} = [expdir_bases{expdiri},': ',sprintf('%s ',data(expdiri).(examinestats{stati}{1}){:})];
          else
            s{2+expdirii} = [expdir_bases{expdiri},': ',data(expdiri).(examinestats{stati}{1})];
          end
        end
      else
        s{2} = sprintf('mean(%s) = %s = %s std',statnames{stati},...
          num2str(stat(linei,stati)),num2str(normstat(linei,stati)));
        for expdirii = 1:numel(expdiris),
          expdiri = expdiris(expdirii);
          s{2+expdirii} = [expdir_bases{expdiri},': ',num2str(expstat(expdiri,stati))];
        end
        
      end
    end
    
  end

%% button down function

  function ButtonDownFcn(varargin)
        
    try
      
      if ~isempty(stati_selected),
        set(hx(stati_selected),'Color','k','FontWeight','normal');
      end

      
      % get current point
      tmp = get(hax,'CurrentPoint');
      xclicked = tmp(1,1);
      yclicked = tmp(1,2);
      tmpidx = find(idx_visible);
      tmpn = numel(tmpidx);
      [d,closest] = min( ((reshape(x(idx_visible,:),[1,tmpn*nstats])-xclicked)/dx).^2+...
        ((reshape(normstat(idx_visible,:),[1,tmpn*nstats])-yclicked)/dy).^2 );
      d = sqrt(d);
      
      if d > examine_params.maxdistclick,
        set(hselected,'Visible','off');
        set(hselected1,'Visible','off');
        set(hmenu.open,'Enable','off');
        stati_selected = [];
        linei_selected = [];
        set(htext,'String','');
        set([hmanual_behavior.text,hmanual_behavior.popup,hmanual_behavior.pushbutton,hmanual_behavior.pushbutton_info],'Visible','off');
        return;
      end
      
      SelectionType = get(hfig,'SelectionType');
      set([hmanual_behavior.text,hmanual_behavior.popup,hmanual_behavior.pushbutton,hmanual_behavior.pushbutton_info],'Visible','on');

      [linei_selected,stati_selected] = ind2sub([tmpn,nstats],closest);
      linei_selected = tmpidx(linei_selected);

      UpdateSelected();
      
      if strcmp(SelectionType,'open'),

        open_Callback();

      end
      
    catch ME,

      fprintf('Error evaluating buttondownfcn, disabling:\n');
      set(hax,'ButtonDownFcn','');
      rethrow(ME);
    end
    
  end

%% update selected line stuff

  function UpdateSelected()
    
    UpdateStatSelected();
    
    set(hselected,'XData',x(linei_selected,:),'YData',normstat(linei_selected,:),'Visible','on');
    set(hmenu.open,'Enable','on');
    
    % TODO!
    manual_behavior_curr = manual_behavior(linei_selected);
    s = get(hmanual_behavior.popup,'String');
    vcurr = find(strncmpi(manual_behavior_curr,s,1),1);
    if isempty(vcurr),
      error('Unknown manual_behavior %s',manual_behavior_curr);
    end
    set(hmanual_behavior.popup,'Value',vcurr);
    
  end
      
%% update stat selected-based stuff

  function UpdateStatSelected()
    
    %set(hfig,'Interruptible','off');
    s = printfun(linei_selected,stati_selected);
    set(htext,'String',s);
    set(hselected1,'XData',x(linei_selected,stati_selected),...
      'YData',normstat(linei_selected,stati_selected),'visible','on');
    set(hmenu.open,'Enable','on');
    set(hx,'Color','k','FontWeight','normal');
    set(hx(stati_selected),'Color','r','FontWeight','bold');
    %set(hfig,'Interruptible','on');

  end

%% manual_behavior = n callback

  function plot_manual_n_Callback(hObject,event) %#ok<INUSD>
    
    v = get(hObject,'Checked');
    set(h(idx_manual_n),'Visible',v);
    if strcmpi(v,'on'),
      set(hObject,'Checked','off');
      idx_visible(idx_manual_n) = false;
      set(h(~idx_visible),'Visible','off');
    else
      set(hObject,'Checked','on');
      idx_visible(idx_manual_n) = true;
      set(h(idx_visible),'Visible','on');
    end
    
  end

%% manual_behavior = d callback

  function plot_manual_d_Callback(hObject,event) %#ok<INUSD>
    
    v = get(hObject,'Checked');
    set(h(idx_manual_d),'Visible',v);
    if strcmpi(v,'on'),
      set(hObject,'Checked','off');
      idx_visible(idx_manual_d) = false;
      set(h(~idx_visible),'Visible','off');
    else
      set(hObject,'Checked','on');
      idx_visible(idx_manual_d) = true;
      set(h(idx_visible),'Visible','on');
    end
    
  end

%% go to next date range

  function nextdate_Callback(hObject,event)
    
    if needsave,
      answer = questdlg('Save manual_behavior and notes?');
      switch lower(answer),
        case 'cancel',
          return;
        case 'yes',
          save_Callback(hObject,event);
      end
    end
    
    maxdatenum = maxdatenum + period;
    handles = FlyBowlExamineBehaviorVariables(...
      'analysis_protocol',analysis_protocol,...
      'settingsdir',settingsdir,...
      'datalocparamsfilestr',datalocparamsfilestr,...
      'hfig',hfig,...
      'period',period,...
      'maxdatenum',maxdatenum,...
      'figpos',get(hfig,'Position'),...
      'datenumnow',datenumnow,...
      'sage_params_path',sage_params_path,...
      'sage_db',sage_db,...
      'username',rc.username,...
      'rootdatadir',rootdatadir,...
      'loadcacheddata',loadcacheddata,...
      'groupnames',groupnames);

  end

%% go to previous date range

  function prevdate_Callback(hObject,event)
    
    if needsave,
      answer = questdlg('Save manual_behavior and notes?');
      switch lower(answer),
        case 'cancel',
          return;
        case 'yes',
          save_Callback(hObject,event);
      end
    end
   
    maxdatenum = maxdatenum - period;
    handles = FlyBowlExamineBehaviorVariables(...
      'analysis_protocol',analysis_protocol,...
      'settingsdir',settingsdir,...
      'datalocparamsfilestr',datalocparamsfilestr,...
      'hfig',hfig,...
      'period',period,...
      'maxdatenum',maxdatenum,...
      'figpos',get(hfig,'Position'),...
      'datenumnow',datenumnow,...
      'sage_params_path',sage_params_path,...
      'sage_db',sage_db,...
      'username',rc.username,...
      'rootdatadir',rootdatadir,...
      'loadcacheddata',loadcacheddata,...
      'groupnames',groupnames);

  end

%% set manual pf callback
  
  function manual_behavior_popup_Callback(hObject,event) %#ok<INUSD>
    
    if isempty(linei_selected),
      warning('No line selected');
      return;
    end

    addToUndoList();
    
    s = get(hObject,'String');
    vcurr = get(hObject,'Value');
    manual_behavior_full = s{vcurr};
    manual_behavior_new = manual_behavior_full(1);
    for expdiri = find(lineidx==linei_selected),
      data(expdiri).manual_behavior = manual_behavior_new;
    end
    %SetManualPF(data(expdiri_selected),sage_db);
    manual_behavior(linei_selected) = manual_behavior_new;

    switch lower(manual_behavior_new),
      
      case 'n',
        
        % update indices
        idx_manual_n(linei_selected) = true;
        idx_manual_d(linei_selected) = false;

        % update marker
        set(h(linei_selected),'Marker','+');

        % update visible
        set(h(linei_selected),'Visible',get(hmenu.plot_manual_n,'Checked'));

      case 'd',

        % update indices
        idx_manual_n(linei_selected) = false;
        idx_manual_d(linei_selected) = true;

        % update marker
        set(h(linei_selected),'Marker','x');
        
        % update visible
        set(h(linei_selected),'Visible',get(hmenu.plot_manual_d,'Checked'));
        
      case 'u',
        
        % update indices
        idx_manual_n(linei_selected) = false;
        idx_manual_d(linei_selected) = false;
        
        % update marker
        set(h(linei_selected),'Marker','o');

        % update visible
        idx_visible(linei_selected) = true;
        set(h(linei_selected),'Visible','on');
        
      otherwise
        
        error('Unknown manual_behavior value %s',manual_behavior_new);
        
    end
    needsave = true;

    
  end

%% add note callback
 function manual_behavior_pushbutton_Callback(hObject,event) %#ok<INUSD>
    
   hnotes.dialog = dialog('Name','Add notes','WindowStyle','Normal','Resize','on');
   done_pos = [.29,.02,.2,.1];
   cancel_pos = [.51,.02,.2,.1];
   notes_pos = [.02,.14,.96,.84];
   
   % TODO
   if iscell(notes_curation{linei_selected}),
     notes_curation_curr = notes_curation{linei_selected};
   else
     notes_curation_curr = regexp(notes_curation{linei_selected},'\\n','split');
   end

   hnotes.pushbutton_done = uicontrol(hnotes.dialog,'Style','pushbutton',...
     'Units','normalized','Position',done_pos,...
     'String','Done','Callback',@notes_done_Callback);
   hnotes.pushbutton_cancel = uicontrol(hnotes.dialog,'Style','pushbutton',...
     'Units','normalized','Position',cancel_pos,...
     'String','Cancel','Callback',@notes_cancel_Callback);
   hnotes.edit_notes = uicontrol(hnotes.dialog,'Style','edit',...
     'Units','normalized','Position',notes_pos,...
     'Min',0,'Max',25,...
     'String',notes_curation_curr,...
     'HorizontalAlignment','left',...
     'BackgroundColor','w');
   uiwait(hnotes.dialog);
   
 end

  function notes_done_Callback(hObject,event) %#ok<INUSD>
    
    addToUndoList();

    notes_curation{linei_selected} = get(hnotes.edit_notes,'String');
    expdiris = find(lineidx == linei_selected);
    for expdiri = expdiris,
      data(expdiri).notes_curation = notes_curation{linei_selected};
    end
    
    tmpi = find(strcmpi(examinestats,'notes_curation'),1);
    if ~isempty(tmpi),
      s = notes_curation{linei_selected};
      v = ~isempty(s) && ~strcmpi(s,'None');
      v = double(v)*2*3-3; % make -3, 3
      stat(linei_selected,tmpi) = v;
      normstat(linei_selected,tmpi) = ...
        (stat(linei_selected,tmpi) - mu(tmpi))/z(tmpi);
      set(h(linei_selected),'YData',normstat(linei_selected,:));
      set(hselected,'YData',normstat(linei_selected,:));
      set(hselected1,'YData',normstat(linei_selected,stati_selected));
    end
    
    needsave = true;
    
    close(hnotes.dialog);
    
  end

  function notes_cancel_Callback(hObject,event) %#ok<INUSD>
    
    close(hnotes.dialog);
    
  end

  function set_rest_manual_n_Callback(hObject,event) %#ok<INUSD>
    set_rest_manual_Callback('N');
  end

  function set_rest_manual_d_Callback(hObject,event) %#ok<INUSD>
    set_rest_manual_Callback('D');
  end

  function set_rest_manual_Callback(s)

    addToUndoList();
    
    b = questdlg(sprintf('Set all experiments plotted with manual_behavior = U to %s?',s),...
      'Set manual_behavior for rest?');
    if ~strcmpi(b,'Yes'),
      return;
    end
    
    idx = find(lower(manual_behavior) == 'u');
    
    if isempty(idx),
      return;
    end
    
    needsave = true;
    
    for tmpi = idx,
      for tmpj = find(lineidx == tmpi),
        data(tmpj).manual_behavior = s;
        %SetManualPF(data(tmpi),sage_db);
      end
    end
    manual_behavior(idx) = s;
    
    if lower(s) == 'n',

      % update indices
      idx_manual_n(idx) = true;
      idx_manual_d(idx) = false;
      
      % update marker
      set(h(idx),'Marker','+');
      
      % update visible
      set(h(idx),'Visible',get(hmenu.plot_manual_n,'Checked'));
      
    else

      % update indices
      idx_manual_n(idx) = false;
      idx_manual_d(idx) = true;
      
      % update marker
      set(h(idx),'Marker','x');
      
      % update visible
      set(h(idx),'Visible',get(hmenu.plot_manual_d,'Checked'));
      
    end
    
  end

  function save_Callback(hObject,event) %#ok<INUSD>
    
    while true
      [savefilename1,savepath1] = uiputfile(savefilename,'Save manual_behavior tsv');
      if ~ischar(savefilename1),
        return;
      end
      savefilename = fullfile(savepath1,savefilename1);
      rc.savepath = savepath1;
      %[savepath2,savefilename2] = fileparts(savefilename);
      %savefilename2 = fullfile(savepath2,[savefilename2,'_diagnosticinfo.mat']);
      
      fid = fopen(savefilename,'w');
      if fid < 0,
        warndlg(sprintf('Could not open file %s for writing. Make sure it is not open in another program.',savefilename),'Could not save');
        continue;
      end
      break;
    end
    fprintf(fid,'#line_name\texperiment_name\tmanual_behavior\tnotes_curation\tdiagnostic_fields\n');
    % TODO
    for tmpi = 1:nlines,
      for tmpj = find(lineidx == tmpi),
        if iscell(data(tmpj).notes_curation),
          notes_curation_curr = sprintf('%s\\n',data(tmpj).notes_curation{:});
        else
          notes_curation_curr = data(tmpj).notes_curation;
        end
        fprintf(fid,'%s\t%s\t%s\t%s\t',...
          data(tmpj).line_name,...
          data(tmpj).experiment_name,...
          data(tmpj).manual_behavior,...
          notes_curation_curr);
        % also print info
        idx = info(tmpi,:);
        if ~any(idx),
          tmps = '';
        else
          tmps = sprintf('%s,',statnames{idx});
          tmps = tmps(1:end-1);
        end
        fprintf(fid,'%s\n',tmps);
      end
    end
    fclose(fid);
    %save(savefilename2,'info','data','examinestats','linenames');
    needsave = false;
    
  end

  function close_fig_Callback(hObject,event)
    
    if needsave,
      answer = questdlg('Save manual_behavior and notes?');
      switch lower(answer),
        case 'cancel',
          return;
        case 'yes',
          save_Callback(hObject,event);
      end
    end

    try
      save(examine_params.rcfile,'-struct','rc');
    catch ME,
      warning('Error saving rc file:\n %s',getReport(ME));
    end
    
    % delete temporary files
    try
    for tmpi = 1:numel(tempfilescreated);
      if exist(tempfilescreated{tmpi},'file'),
        delete(tempfilescreated{tmpi});
      end
    end
    catch ME
      warning('Error deleting temporary files:\n %s',getReport(ME));
    end
    
    if ishandle(hObject),
      delete(hObject);
    end
    
  end

%% add state to undo list

  function addToUndoList()

    if numel(undolist) >= maxnundo,
      undolist = undolist(end-maxnundo+1:end);
    end
    undolist = structappend(undolist,...
      struct('manual_behavior',{[data.manual_behavior]},...
      'notes_curation',{{data.notes_curation}},...
      'line_manual_behavior',{manual_behavior},...
      'line_notes_curation',{notes_curation},...
      'info',{info}));
    set(hmenu.undo,'Enable','on');
  end


%% reset state

  function setManualPFState(s)
  
    manual_behavior = s.line_manual_behavior;
    notes_curation = s.line_notes_curation;
    for tmpi = 1:nexpdirs,
      data(tmpi).manual_behavior = s.manual_behavior(tmpi);
      data(tmpi).notes_curation = s.notes_curation{tmpi};
    end
    idx_manual_n = lower(manual_behavior) == 'n';
    idx_manual_d = lower(manual_behavior) == 'd';
    idx_visible = true(1,nlines);
    if strcmpi(get(hmenu.plot_manual_n,'Checked'),'off'),
      idx_visible(idx_manual_n) = false;
    end
    if strcmpi(get(hmenu.plot_manual_d,'Checked'),'off'),
      idx_visible(idx_manual_d) = false;
    end
    set(h,'Marker','o');
    set(h(idx_manual_n),'Marker','+');
    set(h(idx_manual_d),'Marker','x');
    set(h(idx_visible),'Visible','on');
    set(h(~idx_visible),'Visible','off');
    
    if ~isempty(linei_selected),
      manual_behavior_curr = manual_behavior(linei_selected);
      ss = get(hmanual_behavior.popup,'String');
      vcurr = find(strncmpi(manual_behavior_curr,ss,1),1);
      if isempty(vcurr),
        error('Unknown manual_behavior %s',manual_behavior_curr);
      end
      set(hmanual_behavior.popup,'Value',vcurr);
    end
    
    info = s.info;
    
  end

  function undo_Callback(hObject,event) %#ok<INUSD>
    
    if isempty(undolist), return; end
    setManualPFState(undolist(end));
    undolist = undolist(1:end-1);
    if isempty(undolist),
      set(hmenu.undo,'Enable','off');
    end
    
  end

%% add info callback
 function manual_behavior_pushbutton_info_Callback(hObject,event) %#ok<INUSD>

   ncstats = 4;
   nrstats = ceil(nstats/ncstats);
   checkbox_h = 20;
   checkbox_w = 200;
   border = 10;
   donebutton_h = 20;
   donebutton_w = 60;
   dialog_h = nrstats*checkbox_h+border*3+donebutton_h;
   dialog_w = ncstats*checkbox_w+border*2;
   
   hinfo.dialog = dialog('Name','Add information','WindowStyle','Normal',...
     'Resize','on','Units','pixels','CloseRequestFcn',@info_done_Callback);
   SetFigureSize(hinfo.dialog,dialog_w,dialog_h);
   hinfo.checkboxes = nan(1,nstats);
   
   for tmpi = 1:nstats,
     [tmpr,tmpc] = ind2sub([nrstats,ncstats],tmpi);
     checkbox_pos = [border+(tmpc-1)*checkbox_w,dialog_h-(border+tmpr*checkbox_h),...
       checkbox_w,checkbox_h];
     hinfo.checkboxes(tmpi) = uicontrol(hinfo.dialog,'Style','checkbox',...
       'String',statnames{tmpi},'Position',checkbox_pos,...
       'Value',info(linei_selected,tmpi));
   end
   hinfo.donebutton = uicontrol(hinfo.dialog,'Style','pushbutton',...
     'String','Done','Position',...
     [dialog_w/2-donebutton_w/2,border,donebutton_w,donebutton_h],...
     'Callback',@info_done_Callback);
   uiwait(hinfo.dialog);
   
 end

  function info_done_Callback(hObject,event) %#ok<INUSD>
    
    for tmpi = 1:nstats,
      info(linei_selected,tmpi) = get(hinfo.checkboxes(tmpi),'Value') == 1;
    end
    needsave = true;
    delete(hinfo.dialog);
    
  end

  function MotionFcn(hObject,event) %#ok<INUSD>
    
    if isempty(linei_selected),
      return;
    end
    
    try
      tmp = get(hax,'CurrentPoint');
      xhover = tmp(1,1);
      if xhover < 0 || xhover > nstats+1,
        return;
      end
      stati_selected = min(nstats,max(1,round(xhover)));
      UpdateStatSelected();
    catch ME,
      fprintf('Error evaluating motionfcn, disabling:\n');
      set(hfig,'WindowButtonMotionFcn','');
      rethrow(ME);
    end
    
  end

%% open experiment callback

  function open_Callback(hObject,event) %#ok<INUSD>
    % open experiment
    if isempty(linei_selected),
      return;
    end
    
    % create an html page with links to all experiments
    filenamecurr = fullfile(tempdir,[linenames{linei_selected},'.html']);
    fid = fopen(filenamecurr,'w');
    if fid >= 0,
      fprintf(fid,'<html>\n<title>%s</title>\n<body>\n',linenames{linei_selected});
      fprintf(fid,'<h1>Line %s</h1>\n',linenames{linei_selected});
      fprintf(fid,'<ul>\n');
      expdiris = find(lineidx == linei_selected);
      for expdiri = expdiris(:)',
        fprintf(fid,'  <li><a href="file://%s">%s</a></li>\n',expdirs{expdiri},expdir_bases{expdiri});
      end
      fprintf(fid,'</ul>\n');
      fprintf(fid,'</body>\n</html>\n');
      fclose(fid);
    end
    if ~exist(filenamecurr,'file'),
      warning('Could not open temporary file %s',filenamecurr);
      return;
    end
    tempfilescreated{end+1} = filenamecurr;
    % open this page
    web(filenamecurr,'-browser');
  end

%% load curation spreadsheet callback

  function load_Callback(hObject,event) %#ok<INUSD>
    
    numel_info = 5;
    numel_noinfo = numel_info-1;
    
    if needsave,
      res = questdlg('Load state from file? All changes will be lost.');
      if ~strcmpi(res,'yes'),
        return;
      end
    end
    addToUndoList();
    
    [loadfilename,loadpath] = uigetfile(savefilename,'Load manual_behavior tsv');
    if ~ischar(loadfilename),
      return;
    end
    loadfilename = fullfile(loadpath,loadfilename);
    
    state = undolist(end);
    experiment_names = {data.experiment_name};

    [loadpath2,loadfilename2] = fileparts(loadfilename);
    loadfilename2 = fullfile(loadpath2,[loadfilename2,'_diagnosticinfo.mat']);
    loadedinfo = false;
    if exist(loadfilename2,'file'),
      try
        tmpinfo = load(loadfilename2,'info','data','examinestats','linenames');
        tmpinfo.experiment_names = {tmpinfo.data.experiment_name};
        newstatnames = SetStatNames(tmpinfo.examinestats);
        if ~isempty(setdiff(statnames,newstatnames)) || ...
            ~isempty(setdiff(newstatnames,statnames)),
          error('examined stats different, not loading diagnostic info');
        end
        [isintersect,idxcurr] = ismember(tmpinfo.linenames,linenames);
        state.info(idxcurr(isintersect),:) = tmpinfo.info(isintersect,:);
      catch ME,
        warning('Could not load info from %s:\n%s',loadfilename2,getReport(ME));
      end
    end
    
    fid = fopen(loadfilename,'r');
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
      if numel(m) < numel_noinfo || numel(m) > numel_info,
        warning('Skipping line %s: wrong number of fields',s);
      end
      
      isinfo = numel(m) == numel_info;
      linename = m{1};
      experiment_name = m{2};
      manual_behavior_curr = m{3};
      notes_curation_curr = m{4};
      notes_curation_curr = regexp(notes_curation_curr,'\\n','split');
      if isinfo,
        info_s = regexp(m{5},',','split');
        info_s = setdiff(info_s,{''});
        info_curr = ismember(statnames,info_s);
        unknown_stats = ~ismember(info_s,statnames);
        if any(unknown_stats),
          warning(['Unknown info stats: ',sprintf('%s ',info_s{unknown_stats})]);
        end
      end
      
      tmpi = find(strcmp(linenames,linename),1);
      if isempty(tmpi),
        fprintf('line %s not currently examined, skipping\n',linename);
        continue;
      end
      
      state.manual_behavior(tmpi) = manual_behavior_curr(1);
      state.notes_curation{tmpi} = notes_curation_curr;
      tmpi = find(strcmp(linename,linenames),1);
      if isempty(tmpi),
        % this should never happen
        fprintf('line %s not currently examined, skipping\n',linename);
      end
      state.line_manual_behavior(tmpi) = manual_behavior_curr(1);
      state.line_notes_curation{tmpi} = notes_curation_curr;
      if isinfo,
        state.info(tmpi,:) = info_curr;
      end
    end
    fclose(fid);

    needsave = true;
    setManualPFState(state);
    
  end

%% search callback

  function search_Callback(hObject,event) %#ok<INUSD>
    
    if ~isfield(hsearch,'dialog') || ~ishandle(hsearch.dialog),
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
      hsearch.dialog = dialog('Name','Find experiments','WindowStyle','normal',...
        'Units','pixels');
      SetFigureSize(hsearch.dialog,dlg_width,dlg_height);
      hsearch.field_text = uicontrol(hsearch.dialog,'Style','text',...
        'Units','pixels',...
        'Position',[gap_width,dlg_height-(line_height+gap_height),prompt_width,line_height],...
        'String','Field: ',...
        'HorizontalAlignment','right');
      groupi = find(strcmp(groupnames,'line_name'),1);
      if isempty(groupi),
        groupi = 1;
      end
      hsearch.field_popupmenu = uicontrol(hsearch.dialog,'Style','popupmenu',...
        'Units','pixels',...
        'Position',[prompt_width+2*gap_width,dlg_height-(line_height+gap_height),answer_width,line_height],...
        'String',groupnames,...
        'Value',groupi,...
        'HorizontalAlignment','right',...
        'Callback',@search_field_Callback);
      hsearch.value_text = uicontrol(hsearch.dialog,'Style','text',...
        'Units','pixels',...
        'Position',[gap_width,dlg_height-2*(line_height+gap_height),prompt_width,line_height],...
        'String','Value: ',...
        'HorizontalAlignment','right');
      if iscell(groupvalues{groupi}),
        group_s = groupvalues{groupi};
      else
        group_s = num2str(groupvalues{groupi}(:));
      end
      hsearch.value_popupmenu = uicontrol(hsearch.dialog,'Style','popupmenu',...
        'Units','pixels',...
        'Position',[prompt_width+2*gap_width,dlg_height-2*(line_height+gap_height),answer_width,line_height],...
        'String',group_s,...
        'Value',1,...
        'HorizontalAlignment','right',...
        'Callback',@search_value_Callback);
      hsearch.findnext_pushbutton = uicontrol(hsearch.dialog,'Style','pushbutton',...
        'Units','pixels',...
        'Position',[dlg_width/2-button_width/2,dlg_height-(2*line_height+button_height+3*gap_height),button_width,button_height],...
        'String','Find next',...
        'HorizontalAlignment','center',...
        'Callback',@search_findnext_Callback);
      
    else
      figure(hsearch.dialog);
    end
    
    hsearch.resulti = 0;
    
  end

  function search_field_Callback(varargin)
    
    groupi = get(hsearch.field_popupmenu,'Value');
    if iscell(groupvalues{groupi}),
      group_s = groupvalues{groupi};
    else
      group_s = num2str(groupvalues{groupi}(:));
    end
    set(hsearch.value_popupmenu,'String',group_s,'Value',1);
    hsearch.resulti = 0;
    
  end

  function search_value_Callback(varargin)
    
    hsearch.resulti = 0;
    
  end

  function search_findnext_Callback(varargin)
    
    groupi = get(hsearch.field_popupmenu,'Value');
    groupname = groupnames{groupi};
    groupvaluei = get(hsearch.value_popupmenu,'Value');
    % need to do the search
    if hsearch.resulti == 0,
      if iscell(groupvalues{groupi}),
        results_exp = strcmp({data.(groupname)},groupvalues{groupi}{groupvaluei});
      else
        results_exp = find([data.(groupname)] == groupvalues{groupi}(groupvaluei));
      end
      hsearch.results = unique(lineidx(results_exp));
    end
    
    hsearch.resulti = hsearch.resulti + 1;
    if hsearch.resulti > numel(hsearch.results),
      hsearch.resulti = 1;
    end
    
    linei_selected = hsearch.results(hsearch.resulti);
    fprintf('Selecting line %s\n',linenames{linei_selected});
    
    UpdateSelected();
    
  end


end

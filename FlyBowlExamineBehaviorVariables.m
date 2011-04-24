function handles = FlyBowlExamineBehaviorVariables(varargin)

handles = struct;

[analysis_protocol,settingsdir,datalocparamsfilestr,hfig,period,maxdatenum,...
  figpos,datenumnow,sage_params_path,sage_db,username,rootdatadir,leftovers] = ...
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
  'rootdatadir','/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data');

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
  if ~isempty(username),
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
note_metadata_fns = {'notes_behavioral','notes_technical'};
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

queries = leftovers;
queries(end+1:end+2) = {'daterange',daterange};
queries(end+1:end+2) = {'data_type','stats_perframe_*'};
queries(end+1:end+2) = {'flag_aborted',0};
queries(end+1:end+2) = {'automated_pf','P'};
queries(end+1:end+2) = {'experiment_name','FlyBowl_*'};
%data = SAGEGetBowlData(queries{:},'removemissingdata',true,'rootdir',rootdatadir);
load('datacache20110420.mat','data');
if isempty(data),
  uiwait(warndlg(sprintf('No data for date range %s to %s',daterange{:}),'No data found'));
  if didinputweek,
    fprintf('Trying previous week...\n');
    maxdatenum = maxdatenum - period;
    handles = FlyBowlExamineBehaviorVariables(...
      'analysis_protocol',analysis_protocol,...
      'settingsdir',settingsdir,...
      'datalocparamsfilestr',datalocparamsfilestr,...
      'period',period,...
      'maxdatenum',maxdatenum,...
      'datenumnow',datenumnow,...
      'sage_params_path',sage_params_path,...
      'sage_db',sage_db,...
      'username',rc.username);
    return;
  else
    fprintf('Choose a different week.\n');
    handles = FlyBowlExamineBehaviorVariables(...
      'analysis_protocol',analysis_protocol,...
      'settingsdir',settingsdir,...
      'datalocparamsfilestr',datalocparamsfilestr,...
      'period',period,...
      'maxdatenum',[],...
      'datenumnow',datenumnow,...
      'sage_params_path',sage_params_path,...
      'sage_db',sage_db,...
      'username',rc.username);
    return;
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
for i = 1:nexpdirs,
  if ~isfield(data,'file_system_path'),
    data(i).file_system_path = fullfile(rootdatadir,expdir_bases{i});
  end
  if ~isfield(data,'manual_behavior'),
    data(i).manual_behavior = 'U';
  end
end
expdirs = {data.file_system_path};


%% get behavior variables

nstats = numel(examinestats);
% average over lines

[linenames,~,lineidx] = unique({data.line_name});
nlines = numel(linenames);
stat = nan(nlines,nstats);
expstat = nan(nexpdirs,nstats);

for i = 1:nstats,

  % flags are now binary
  % special case: flags
  if ismember(examinestats{i}{1},flag_metadata_fns),
    % make -3, 3
    v = double([data.(examinestats{i}{1})])*2*3-3;
  end
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
  if ismember(examinestats{i}{1},note_metadata_fns),
    
    v = cellfun(@(s) ~isempty(s) && ~strcmpi(s,'None'),{data.(examinestats{i}{1})});
    v = double(v)*2*3-3; % make -3, 3
    % take average over lines
    for j = 1:nlines,
      stat(j,i) = nanmean(v(lineidx==j));
    end
    %stat(:,i) = v;

  % special case: strings
  elseif ismember(examinestats{i}{1},string_metadata_fns),
    [~,~,v] = unique({data.(examinestats{i}{1})});
    v = v - 1;
    v = (v/max(v)*2-1)*3;
    for j = 1:nlines,
      stat(j,i) = nanmean(v(lineidx==j));
    end
    %stat(:,i) = v;

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

statnames = cell(1,nstats);
for i = 1:nstats,
  
  statnames{i} = sprintf('%s_',examinestats{i}{:});
  statnames{i} = statnames{i}(1:end-1);
  statnames{i} = strrep(statnames{i},'stats_perframe_','');
  statnames{i} = strrep(statnames{i},'meanmean_perexp_','');
  statnames{i} = strrep(statnames{i},'flyany_frame','');
  
end


%% z-score stats

z = sig;
z(z == 0) = 1;
normstat = bsxfun(@rdivide,bsxfun(@minus,stat,mu),z);

%% create figure

if ishandle(hfig),
  delete(hfig);
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
  'Callback',@undo_Callback,'Enable','off');
hmenu.save = uimenu(hmenu.file,'Label','Save...',...
  'Callback',@save_Callback);
hmenu.load = uimenu(hmenu.file,'Label','Load Spreadsheet...',...
  'Callback',@load_Callback);
hmenu.open = uimenu(hmenu.file,'Label','View Experiment',...
  'Callback',@open_Callback,'Enable','off');

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

%% print function

  function s = printfun(linei,stati)
    
    expdiris = find(lineidx == linei);
    s = cell(1,2+numel(expdiris));
    s{1} = linenames{linei};
    % special case: notes, strings
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
      
      if strcmp(SelectionType,'open'),
        
        % create an html page with links to all experiments
        filenamecurr = fullfile(tempdir,[linenames{linei_selected},'.html']);
        fid = fopen(filenamecurr,'w');
        fprintf(fid,'<html>\n<title>%s</title>\n<body>\n',linenames{linei_selected});
        fprintf(fid,'<h1>Line %s</h1>\n',linenames{linei_selected});
        fprintf(fid,'<ul>\n');
        expdiris = find(lineidx == linei_selected);
        for expdiri = expdiris(:)',
          fprintf(fid,'  <li><a href="%s">%s</a></li>\n',expdirs{expdiri},expdir_bases{expdiri});
        end
        fprintf(fid,'</ul>\n');
        fprintf(fid,'</body>\n</html>\n');
        fclose(fid);
        tempfilescreated{end+1} = filenamecurr;
        % open this page
        web(filenamecurr,'-browser');
      end
      
    catch ME,

      fprintf('Error evaluating buttondownfcn, disabling:\n');
      set(hax,'ButtonDownFcn','');
      rethrow(ME);
    end
    
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
      'username',rc.username);

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
      'username',rc.username);

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
   notes_behavioral_pos = [.02,.14,.96,.35];
   text_behavioral_pos = [.02,.49,.96,.06];
   notes_technical_pos = [.02,.57,.96,.35];
   text_technical_pos = [.02,.92,.96,.06];
   
   % TODO
   if iscell(data(expdiri_selected).notes_technical),
     notes_technical = data(expdiri_selected).notes_technical;
   else
     notes_technical = regexp(data(expdiri_selected).notes_technical,'\\n','split');
   end
   if iscell(data(expdiri_selected).notes_behavioral),
     notes_behavioral = data(expdiri_selected).notes_behavioral;
   else
     notes_behavioral = regexp(data(expdiri_selected).notes_behavioral,'\\n','split');
   end

   hnotes.pushbutton_done = uicontrol(hnotes.dialog,'Style','pushbutton',...
     'Units','normalized','Position',done_pos,...
     'String','Done','Callback',@notes_done_Callback);
   hnotes.pushbutton_cancel = uicontrol(hnotes.dialog,'Style','pushbutton',...
     'Units','normalized','Position',cancel_pos,...
     'String','Cancel','Callback',@notes_cancel_Callback);
   hnotes.edit_behavioral = uicontrol(hnotes.dialog,'Style','edit',...
     'Units','normalized','Position',notes_behavioral_pos,...
     'Min',0,'Max',10,...
     'String',notes_behavioral,...
     'HorizontalAlignment','left',...
     'BackgroundColor','w');
   hnotes.text_behavioral = uicontrol(hnotes.dialog,'Style','text',...
     'Units','normalized','Position',text_behavioral_pos,...
     'String','Behavior notes:',...
     'HorizontalAlignment','left');
   hnotes.edit_technical = uicontrol(hnotes.dialog,'Style','edit',...
     'Units','normalized','Position',notes_technical_pos,...
     'Min',0,'Max',10,...
     'String',notes_technical,...
     'HorizontalAlignment','left',...
     'BackgroundColor','w');
   hnotes.text_technical = uicontrol(hnotes.dialog,'Style','text',...
     'Units','normalized','Position',text_technical_pos,...
     'String','Technical notes:',...
     'HorizontalAlignment','left');
   uiwait(hnotes.dialog);
   
 end

  function notes_done_Callback(hObject,event) %#ok<INUSD>
    
    addToUndoList();
    
    data(expdiri_selected).notes_technical = get(hnotes.edit_technical,'String');
    data(expdiri_selected).notes_behavioral = get(hnotes.edit_behavioral,'String');
    %SetNotes(data(expdiri_selected),sage_db);
    
    tmpi = find(strcmpi(examinestats,'notes_behavioral'),1);
    if ~isempty(tmpi),
      s = data(expdiri_selected).notes_behavioral;
      v = ~isempty(s) && ~strcmpi(s,'None');
      stat(expdiri_selected,tmpi) = double(v);
      normstat(expdiri_selected,tmpi) = ...
        (stat(linei_selected,tmpi) - mu(tmpi))/z(tmpi);
      set(h(linei_selected),'YData',normstat(linei_selected,:));
      set(hselected,'YData',normstat(linei_selected,:));
      set(hselected1,'YData',normstat(linei_selected,stati_selected));
    end

    tmpi = find(strcmpi(examinestats,'notes_technical'),1);
    if ~isempty(tmpi),
      s = data(expdiri_selected).notes_technical;
      v = ~isempty(s) && ~strcmpi(s,'None');
      stat(linei_selected,tmpi) = double(v);
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
    set_rest_manual_Callback('P');
  end

  function set_rest_manual_d_Callback(hObject,event) %#ok<INUSD>
    set_rest_manual_Callback('F');
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
      data(tmpi).manual_behavior = s;
      %SetManualPF(data(tmpi),sage_db);
    end
    manual_behavior(idx) = s;
    
    if lower(s) == 'p',

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
    
    [savefilename1,savepath1] = uiputfile(savefilename,'Save manual_behavior tsv');
    if ~ischar(savefilename1),
      return;
    end
    savefilename = fullfile(savepath1,savefilename1);
    rc.savepath = savepath1;
    [savepath2,savefilename2] = fileparts(savefilename);
    savefilename2 = fullfile(savepath2,[savefilename2,'_diagnosticinfo.mat']);

    fid = fopen(savefilename,'w');
    fprintf(fid,'#experiment_name\tmanual_behavior\tnotes_behavioral\tnotes_technical\n');
    % TODO
    for tmpi = 1:nexpdirs,
      if iscell(data(tmpi).notes_behavioral),
        notes_behavioral = sprintf('%s\\n',data(tmpi).notes_behavioral{:});
      else
        notes_behavioral = data(tmpi).notes_behavioral;
      end
       if iscell(data(tmpi).notes_technical),
        notes_technical = sprintf('%s\\n',data(tmpi).notes_technical{:});
      else
        notes_technical = data(tmpi).notes_technical;
      end
      fprintf(fid,'%s\t%s\t%s\t%s\n',...
        data(tmpi).experiment_name,...
        data(tmpi).manual_behavior,...
        notes_behavioral,...
        notes_technical);
    end
    fclose(fid);
    save(savefilename2,'info','data','examinestats');
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
    
    save(examine_params.rcfile,'-struct','rc');
    
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
      'notes_behavioral',{{data.notes_behavioral}},...
      'notes_technical',{{data.notes_technical}}));
    set(hmenu.undo,'Enable','on');
  end


%% reset state

  function setManualPFState(s)
  
    manual_behavior = s.manual_behavior;
    for tmpi = 1:nexpdirs,
      data(tmpi).manual_behavior = s.manual_behavior(tmpi);
      data(tmpi).notes_behavioral = s.notes_behavioral{tmpi};
      data(tmpi).notes_technical = s.notes_technical{tmpi};
    end
    idx_manual_n = lower(manual_behavior) == 'p';
    idx_manual_d = lower(manual_behavior) == 'f';
    idx_visible = true(1,nexpdirs);
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
    
    if ~isempty(expdiri_selected),
      manual_behavior_curr = data(expdiri_selected).manual_behavior;
      s = get(hmanual_behavior.popup,'String');
      vcurr = find(strncmpi(manual_behavior_curr,s,1),1);
      if isempty(vcurr),
        error('Unknown manual_behavior %s',manual_behavior_curr);
      end
      set(hmanual_behavior.popup,'Value',vcurr);
    end
    
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
   dialog_pos = get(hinfo.dialog,'Position');
   dialog_pos(1) = dialog_pos(1) - (dialog_w-dialog_pos(3))/2;
   dialog_pos(2) = dialog_pos(2) - (dialog_h-dialog_pos(4))/2;
   dialog_pos(3) = dialog_w;
   dialog_pos(4) = dialog_h;
   set(hinfo.dialog,'Position',dialog_pos);
   hinfo.checkboxes = nan(1,nstats);
   
   for tmpi = 1:nstats,
     [tmpr,tmpc] = ind2sub([nrstats,ncstats],tmpi);
     checkbox_pos = [border+(tmpc-1)*checkbox_w,dialog_h-(border+tmpr*checkbox_h),...
       checkbox_w,checkbox_h];
     hinfo.checkboxes(tmpi) = uicontrol(hinfo.dialog,'Style','checkbox',...
       'String',statnames{tmpi},'Position',checkbox_pos,...
       'Value',info(expdiri_selected,tmpi));
   end
   hinfo.donebutton = uicontrol(hinfo.dialog,'Style','pushbutton',...
     'String','Done','Position',...
     [dialog_w/2-donebutton_w/2,border,donebutton_w,donebutton_h],...
     'Callback',@info_done_Callback);
   uiwait(hinfo.dialog);
   
 end

  function info_done_Callback(hObject,event) %#ok<INUSD>
    
    for tmpi = 1:nstats,
      info(expdiri_selected,tmpi) = get(hinfo.checkboxes(tmpi),'Value') == 1;
    end
    needsave = true;
    delete(hinfo.dialog);
    
  end

  function MotionFcn(hObject,event) %#ok<INUSD>
    
    if isempty(expdiri_selected),
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
    if isempty(expdiri_selected),
      return;
    end
    if ispc,
      winopen(expdirs{expdiri_selected});
    else
      web(expdirs{expdiri_selected},'-browser');
    end
  end

%% load curation spreadsheet callback

  function load_Callback(hObject,event) %#ok<INUSD>
    
    if needsave,
      res = questdlg('Load state from %s? All changes will be lost.');
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
    [loadpath2,loadfilename2] = fileparts(loadfilename);
    loadfilename2 = fullfile(loadpath2,[loadfilename2,'_diagnosticinfo.mat']);
    loadedinfo = false;
    if ~exist(loadfilename2,'file'),
      warndlg(sprintf('diagnostic info file %s does not exist',loadfilename2));
    else
      try
        tmpinfo = load(loadfilename2,'info','data','examinestats');
        tmpinfo.experiment_names = {tmpinfo.data.experiment_name};
        loadedinfo = true;
        if numel(examinestats) ~= numel(tmpinfo.examinestats),
          warning('examined stats different, not loading diagnostic info');
        end
      catch ME,
        warning('Could not load info from %s:\n%s',loadfilename2,getReport(ME));
      end
    end

    state = undolist(end);
    
    fid = fopen(loadfilename,'r');
    experiment_names = {data.experiment_name};
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
      if numel(m) ~= 4,
        warning('Skipping line %s: wrong number of fields',s);
      end
      
      experiment_name = m{1};
      manual_behavior_curr = m{2};
      notes_behavioral = m{3};
      notes_technical = m{4};
      notes_behavioral = regexp(notes_behavioral,'\\n','split');
      notes_technical = regexp(notes_technical,'\\n','split');
      
      tmpi = find(strcmp(experiment_names,experiment_name),1);
      if isempty(tmpi),
        fprintf('experiment %s not currently examined, skipping\n',experiment_name);
        continue;
      end
      
      state.manual_behavior(tmpi) = manual_behavior_curr(1);
      state.notes_behavioral{tmpi} = notes_behavioral;
      state.notes_technical{tmpi} = notes_technical;
    
    end
    fclose(fid);

    if loadedinfo,
      [isintersect,idxcurr] = ismember(tmpinfo.experiment_names,experiment_names);
      info(idxcurr(isintersect),:) = tmpinfo.info(isintersect,:);
    end

    needsave = true;
    setManualPFState(state);
    
  end

end

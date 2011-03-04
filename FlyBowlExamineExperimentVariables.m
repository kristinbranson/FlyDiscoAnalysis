function handles = FlyBowlExamineExperimentVariables(varargin)

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

if isempty(datenumnow),
  datenumnow = now;
  daycurr = weekday(datenumnow);
  daywant = 7;
  datenumnow = floor(datenumnow)+daywant-daycurr;
end

if isempty(maxdatenum),
  maxdatenum = datenumnow;
end
format = 'yyyy-mm-ddTHH:MM:SS';
mindatenum = maxdatenum - period;
mindatestr = datestr(mindatenum,format);
maxdatestr = datestr(maxdatenum,format);
daterange = {mindatestr,maxdatestr};
%maxdaysprev = datenumnow-mindatenum;
mindaysprev = datenumnow-maxdatenum;

%% read parameters

datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
dataloc_params = ReadParams(datalocparamsfile);
examineparamsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.examineexperimentvariablesparamsfilestr);
examine_params = ReadParams(examineparamsfile);
examinefile = fullfile(settingsdir,analysis_protocol,dataloc_params.examineexperimentvariablesfilestr);
examinestats = ReadExamineExperimentVariables(examinefile);
registrationparamsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.registrationparamsfilestr);
registration_params = ReadParams(registrationparamsfile);

if exist(examine_params.rcfile,'file'),
  rc = load(examine_params.rcfile);
else
  rc = struct('username','','savepath','');
end

%% get user name

if isempty(username),
  tmpusername = inputdlg('User name:','Set user name',1,{rc.username});
  if isempty(tmpusername),
    delete(hfig);
    return;
  end
  rc.username = tmpusername{1};
else
  rc.username = username;
end
savefilename = fullfile(rc.savepath,sprintf('DataCuration_FlyBowl_%s_%sto%s.tsv',...
  rc.username,datestr(mindatenum,'yyyymmdd'),datestr(maxdatenum,'yyyymmdd')));
needsave = false;

%% string metadata

flag_metadata_fns = {'flag_redo','flag_review'};
note_metadata_fns = {'notes_behavioral','notes_technical'};
flag_options = {{'None',''},{'Rearing problem','Flies look sick',...
  'See behavioral notes','See technical notes'}};

%% connect to SAGE for storing manualpf
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
data_types = {'ufmf_diagnostics_summary_*','temperature_stream',...
  'ctrax_diagnostics_*','registrationdata_*','sexclassifier_diagnostics_*',...
  'stats_perframe_areasmooth*'};
%data_types = examinestats;

queries = leftovers;
queries(end+1:end+2) = {'daterange',daterange};
queries(end+1:end+2) = {'data_type',data_types};
queries(end+1:end+2) = {'flag_aborted',0};
queries(end+1:end+2) = {'automated_pf','P'};
queries(end+1:end+2) = {'experiment_name','FlyBowl_*'};
%data = SAGEGetBowlData(queries{:});
load('datacache.mat','data');
% sort by date
date = {data.exp_datetime};
[~,order] = sort(date);
data = data(order);
nexpdirs = numel(data);
expdir_bases = {data.experiment_name};
expdir_bases = cellfun(@(s) regexprep(s,'^FlyBowl_',''),expdir_bases,'UniformOutput',false);

% %% which experiments?
% 
% subreadfiles = {};
% for i = 1:numel(examine_params.requiredfiles),
%   subreadfiles{end+1} = dataloc_params.(examine_params.requiredfiles{i}); %#ok<AGROW>
% end
% subreadfiles = unique(subreadfiles);
% 
% exp_params = [{'rootdir',dataloc_params.rootreaddir,'subreadfiles',subreadfiles},leftovers];
% 
% [expdir_bases,expdirs,~,experiments,~,~] = ...
%   getExperimentDirs('settingsdir',settingsdir,...
%   'datalocparamsfilestr',datalocparamsfilestr,...
%   'protocol',analysis_protocol,...
%   exp_params{:});
% nexpdirs = numel(expdirs);
% 
% if nexpdirs == 0,
%   error('No experiments selected');
% end
% 
% % sort by date
% date = {experiments.date};
% [~,order] = sort(date);
% expdirs = expdirs(order);
% %experiments = experiments(order);
% expdir_bases = expdir_bases(order);
% 
% %% load experiment variables for each experiment
% 
% data = [];
% 
% for expdiri = 1:nexpdirs,
%   
%   expdir = expdirs{expdiri};
%   
%   datacurr = struct;
%   
%   % read ufmf diagnostics
%   ufmfdiagnosticsfile = fullfile(expdir,dataloc_params.ufmfdiagnosticsfilestr);
%   datacurr.ufmf_diagnostics = readUFMFDiagnostics(ufmfdiagnosticsfile);
%   
%   % read quickstats
%   quickstatsfile = fullfile(expdir,dataloc_params.quickstatsfilestr);
%   datacurr.quickstats = ReadParams(quickstatsfile);
%   
%   % read ctrax diagnostics
%   ctraxdiagnosticsfile = fullfile(expdir,dataloc_params.ctraxdiagnosticsfilestr);
%   datacurr.ctrax_diagnostics = ReadParams(ctraxdiagnosticsfile);
%   % add in mean_nsplit
%   datacurr.ctrax_diagnostics.mean_nsplit = ...
%     datacurr.ctrax_diagnostics.sum_nsplit / datacurr.ctrax_diagnostics.nlarge_split;
%   % add in nframes_not_tracked
%   datacurr.ctrax_diagnostics.nframes_not_tracked = ...
%     datacurr.ufmf_diagnostics.summary.nFrames - datacurr.ctrax_diagnostics.nframes_analyzed;
%   
%   % read registration data
%   registrationtxtfile = fullfile(expdir,dataloc_params.registrationtxtfilestr);
%   datacurr.registrationData = ReadParams(registrationtxtfile);
% 
%   % read sex classification diagnostics
%   sexclassifierdiagnosticsfile = fullfile(expdir,dataloc_params.sexclassifierdiagnosticsfilestr);
%   datacurr.sexclassifierdiagnostics = ReadParams(sexclassifierdiagnosticsfile);
% 
%   % read temperature file 
%   temperaturefile = fullfile(expdir,dataloc_params.temperaturefilestr);
%   tempdata = importdata(temperaturefile,',');
%   datacurr.temperature.stream = tempdata(:,2);
%   % add in mean, max, std, maxdiff, nreadings
%   datacurr.temperature.mean = nanmean(datacurr.temperature.stream);
%   datacurr.temperature.max = max(datacurr.temperature.stream);
%   datacurr.temperature.maxdiff = datacurr.temperature.max - min(datacurr.temperature.stream);
%   datacurr.temperature.nreadings = numel(datacurr.temperature.stream);
%   
%   data = structappend(data,datacurr);
%   
% end
% 

%% for now, add on extra diagnostics here. we will store these later
for i = 1:nexpdirs,
  % add in mean_nsplit
  data(i).ctrax_diagnostics_mean_nsplit = ...
    data(i).ctrax_diagnostics_sum_nsplit / data(i).ctrax_diagnostics_nlarge_split;
  % add in nframes_not_tracked
  data(i).ctrax_diagnostics_nframes_not_tracked = ...
    data(i).ufmf_diagnostics_summary_nFrames - data(i).ctrax_diagnostics_nframes_analyzed;
  % add in mean, max, std, maxdiff, nreadings
  data(i).temperature_mean = nanmean(data(i).temperature_stream);
  data(i).temperature_max = max(data(i).temperature_stream);
  data(i).temperature_maxdiff = data(i).temperature_max - min(data(i).temperature_stream);
  data(i).temperature_nreadings = numel(data(i).temperature_stream);
  data(i).file_system_path = fullfile(rootdatadir,expdir_bases{i});
  data(i).registrationdata_pxpermm = data(i).registrationdata_circleRadius / registration_params.circleRadius_mm;
end
expdirs = {data.file_system_path};

%% get experiment variables

nstats = numel(examinestats);
stat = nan(nexpdirs,nstats);

for i = 1:nstats,

  % special case: flags
  if ismember(examinestats{i}{1},flag_metadata_fns),
    
    % which set of flags
    v = zeros(1,nexpdirs);
    for j = 1:numel(flag_options),
      v(ismember({data.(examinestats{i}{1})},flag_options{j})) = j;
    end
    stat(:,i) = v;
    
  % special case: notes
  elseif ismember(examinestats{i}{1},note_metadata_fns),
    
    v = cellfun(@(s) ~isempty(s) && ~strcmpi(s,'None'),{data.(examinestats{i}{1})});
    stat(:,i) = double(v);
    
  % numbers
  else
    
    datacurr = {data.(examinestats{i}{1})};
    for k = 2:numel(examinestats{i}),
      datacurr = cellfun(@(s) s.(examinestats{i}{k}),datacurr,'UniformOutput',false);
    end
    badidx = cellfun(@isempty,datacurr);
    for j = find(badidx),
      warning('No data for stat %s experiment %s',sprintf('%s,',examinestats{i}),strtrim(data(j).experiment_name));
    end
    stat(~badidx,i) = cell2mat(datacurr);
    
  end
end
  

% for expdiri = 1:nexpdirs,
%   
%   for i = 1:nstats,
% 
%     pathcurr = examinestats{i};
%   
%     % get examine stats for all expdirs. maybe there is a more efficient way
%     % to do this?
%     statcurr = data(expdiri);
%     for j = 1:numel(pathcurr),
%       statcurr = statcurr.(pathcurr{j});
%     end
%     stat(expdiri,i) = statcurr;
% 
%   end
% 
% end

%% get mean and standard deviation for z-scoring

if isempty(examine_params.examineexperimentvariablesstatsfile),
  mu = nanmean(stat,1);
  sig = nanstd(stat,1,1);  
else
  normstats = load(examine_params.examineexperimentvariablesstatsfile);
  mu = nan(1,nstats);
  sig = nan(1,nstats);
  for i = 1:nstats,

    pathcurr = examinestats{i};
    statcurr = normstats.mu;
    for j = 1:numel(pathcurr),
      statcurr = statcurr.(pathcurr{j});
    end
    mu(i) = statcurr;
    statcurr = normstats.sig;
    for j = 1:numel(pathcurr),
      statcurr = statcurr.(pathcurr{j});
    end
    sig(i) = statcurr;

  end
end

%% abbr names of stats

statnames = cell(1,nstats);
for i = 1:nstats,
  
  statnames{i} = sprintf('%s_',examinestats{i}{:});
  statnames{i} = statnames{i}(1:end-1);
  statnames{i} = strrep(statnames{i},'ufmf_diagnostics_summary','ufmf');
  statnames{i} = strrep(statnames{i},'ufmf_diagnostics_stream','ufmf');
  statnames{i} = strrep(statnames{i},'ctrax_diagnostics','ctrax');
  statnames{i} = strrep(statnames{i},'registrationdata','reg');
  statnames{i} = strrep(statnames{i},'sexclassifier_diagnostics','sex');
  statnames{i} = strrep(statnames{i},'stats_perframe_','');
  statnames{i} = strrep(statnames{i},'flyany_frame','');
  statnames{i} = strrep(statnames{i},'_perexp','');
  
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
colors = jet(nexpdirs)*.7;
drawnow;

%% plot diagnostics

if nexpdirs == 1,
  off = 0;
else
  off = (2*((1:nexpdirs)-(nexpdirs+1)/2)/(nexpdirs-1))*examine_params.offx;
end

h = nan(1,nexpdirs);
x = nan(nexpdirs,nstats);
for expdiri = 1:nexpdirs,
  x(expdiri,:) = (1:nstats)+off(expdiri);
  h(expdiri) = plot(hax,x(expdiri,:),normstat(expdiri,:),'o',...
    'color',colors(expdiri,:),'markerfacecolor',colors(expdiri,:),...
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

%% set manual_pf = p, f marker

%gray_p_color = [.7,.7,.7];
%gray_f_color = [.5,.5,.5];
badidx = cellfun(@isempty,{data.manual_pf});
for j = find(badidx),
  warning('No data for manual_pf, experiment %s',strtrim(data(j).experiment_name));
end
manual_pf = repmat('U',[1,nexpdirs]);
manual_pf(~badidx) = [data.manual_pf];
idx_manual_p = lower(manual_pf) == 'p';
idx_manual_f = lower(manual_pf) == 'f';
idx_visible = true(1,nexpdirs);
set(h(idx_manual_p),'Marker','+');
set(h(idx_manual_f),'Marker','x');

%% selected experiment

hselected = plot(0,0,'o','color','k','Visible','off','HitTest','off','MarkerSize',10,'MarkerFaceColor','k');
expdiri_selected = [];
stati_selected = [];


%% Examine menu

hmenu = struct;
hmenu.file = uimenu('Label','File','Parent',hfig);
hmenu.options = uimenu('Label','Options','Parent',hfig);
hmenu.set = uimenu('Label','Set','Parent',hfig);
hmenu.info = uimenu('Label','Info','Parent',hfig);
hmenu.plot_manual_p = uimenu(hmenu.options,'Label','Plot manual_pf = p',...
  'Checked','on','Callback',@plot_manual_p_Callback);
hmenu.plot_manual_f = uimenu(hmenu.options,'Label','Plot manual_pf = f',...
  'Checked','on','Callback',@plot_manual_f_Callback);
hmenu.set_rest_manual_p = uimenu(hmenu.set,'Label','Set manual_pf == u -> p',...
  'Callback',@set_rest_manual_p_Callback);
hmenu.set_rest_manual_f = uimenu(hmenu.set,'Label','Set manual_pf == u -> f',...
  'Callback',@set_rest_manual_f_Callback);
hmenu.undo = uimenu(hmenu.set,'Label','Undo',...
  'Callback',@undo_Callback,'Enable','off');
hmenu.save = uimenu(hmenu.file,'Label','Save...',...
  'Callback',@save_Callback);
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

%% manual pf

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
manualpf_textpos = [textpos_px(1)+textpos_px(3)+margin,c1-h1/2,w1,h1];
hmanualpf = struct;
hmanualpf.text = uicontrol(hfig,'Style','text','Units','Pixels',...
  'Position',manualpf_textpos,'String','Manual PF:',...
  'BackgroundColor',get(hfig,'Color'),'Visible','off');
manualpf_popuppos = [manualpf_textpos(1)+manualpf_textpos(3)+margin,...
  c1-h2/2,w2,h2];
hmanualpf.popup = uicontrol(hfig,'Style','popupmenu','Units','Pixels',...
  'Position',manualpf_popuppos,'String',{'Pass','Fail','Unknown'},...
  'Value',3,'Visible','off','Callback',@manualpf_popup_Callback);
manualpf_pushbuttonpos = [manualpf_popuppos(1)+manualpf_popuppos(3)+margin,...
  c1-h2/2,w3,h2];
hmanualpf.pushbutton = uicontrol(hfig,'Style','pushbutton','Units','Pixels',...
  'Position',manualpf_pushbuttonpos,'String','Add Note...',...
  'Callback',@manualpf_pushbutton_Callback,...
  'Visible','off');
manualpf_pushbutton_info_pos = [manualpf_pushbuttonpos(1)+manualpf_pushbuttonpos(3)+margin,...
  c1-h2/2,w3,h2];
hmanualpf.pushbutton_info = uicontrol(hfig,'Style','pushbutton','Units','Pixels',...
  'Position',manualpf_pushbutton_info_pos,'String','Add Info...',...
  'Callback',@manualpf_pushbutton_info_Callback,...
  'Visible','off');
set([hmanualpf.text,hmanualpf.popup,hmanualpf.pushbutton,hmanualpf.pushbutton_info],'Units','normalized');

hnotes = struct;
hinfo = struct;
info = false(nexpdirs,nstats);


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
if mindaysprev < 7,
  s = 'this week';
elseif mindaysprev < 14,
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

%% context menu for setting manual_pf
% hcmenu = uicontextmenu;
% hcmenu_manualpf = uimenu(hcmenu,'Label','Set manual_pf...','Callback',@set_manual_pf_Callback);
% set(hax,'uicontextmenu',hcmenu);

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
handles.htext = htext;
handles.hdate = hdate;
% handles.hcmenu_manualpf = hcmenu_manualpf;
% handles.hcmenu = hcmenu;

%% print function

  function s = printfun(expdiri,stati)
    
    s = cell(1,2);
    s{1} = expdir_bases{expdiri};
    % special case: flags
    if ismember(examinestats{stati}{1},flag_metadata_fns) || ...
        ismember(examinestats{stati}{1},note_metadata_fns),
      s{2} = sprintf('%s = %s',statnames{stati},data(expdiri).(examinestats{stati}{1}));
    else
      s{2} = sprintf('%s = %s = %s std',statnames{stati},...
        num2str(stat(expdiri,stati)),num2str(normstat(expdiri,stati)));
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
        stati_selected = [];
        expdiri_selected = [];
        set(htext,'String','');
        set([hmanualpf.text,hmanualpf.popup,hmanualpf.pushbutton,hmanualpf.pushbutton_info],'Visible','off');
        return;
      end
      
      SelectionType = get(hfig,'SelectionType');
      set([hmanualpf.text,hmanualpf.popup,hmanualpf.pushbutton,hmanualpf.pushbutton_info],'Visible','on');

      [expdiri_selected,stati_selected] = ind2sub([tmpn,nstats],closest);
      expdiri_selected = tmpidx(expdiri_selected);

      UpdateStatSelected();
      
      set(hselected,'XData',x(expdiri_selected,:),'YData',normstat(expdiri_selected,:),'Visible','on');
      
      manual_pf_curr = data(expdiri_selected).manual_pf;
      s = get(hmanualpf.popup,'String');
      vcurr = find(strncmpi(manual_pf_curr,s,1),1);
      if isempty(vcurr),
        error('Unknown manual_pf %s',manual_pf_curr);
      end
      set(hmanualpf.popup,'Value',vcurr);
      
      if strcmp(SelectionType,'open'),
        
        % open experiment
        if ispc,
          winopen(expdirs{expdiri_selected});
        else
          web(expdirs{expdiri_selected},'-browser');
        end
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
    s = printfun(expdiri_selected,stati_selected);
    set(htext,'String',s);
      
    set(hx,'Color','k','FontWeight','normal');
    set(hx(stati_selected),'Color','r','FontWeight','bold');
    %set(hfig,'Interruptible','on');

  end

%% manual_pf = p callback

  function plot_manual_p_Callback(hObject,event) %#ok<INUSD>
    
    v = get(hObject,'Checked');
    set(h(idx_manual_p),'Visible',v);
    if strcmpi(v,'on'),
      set(hObject,'Checked','off');
      idx_visible(idx_manual_p) = false;
      set(h(~idx_visible),'Visible','off');
    else
      set(hObject,'Checked','on');
      idx_visible(idx_manual_p) = true;
      set(h(idx_visible),'Visible','on');
    end
    
  end

%% manual_pf = f callback

  function plot_manual_f_Callback(hObject,event) %#ok<INUSD>
    
    v = get(hObject,'Checked');
    set(h(idx_manual_f),'Visible',v);
    if strcmpi(v,'on'),
      set(hObject,'Checked','off');
      idx_visible(idx_manual_f) = false;
      set(h(~idx_visible),'Visible','off');
    else
      set(hObject,'Checked','on');
      idx_visible(idx_manual_f) = true;
      set(h(idx_visible),'Visible','on');
    end
    
  end

%% go to next date range

  function nextdate_Callback(hObject,event)
    
    if needsave,
      answer = questdlg('Save manual_pf and notes?');
      switch lower(answer),
        case 'cancel',
          return;
        case 'yes',
          save_Callback(hObject,event);
      end
    end
    
    maxdatenum = min(datenumnow,maxdatenum + period);
    handles = FlyBowlExamineExperimentVariables(...
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
      answer = questdlg('Save manual_pf and notes?');
      switch lower(answer),
        case 'cancel',
          return;
        case 'yes',
          save_Callback(hObject,event);
      end
    end
   
    maxdatenum = maxdatenum - period;
    handles = FlyBowlExamineExperimentVariables(...
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
  
  function manualpf_popup_Callback(hObject,event) %#ok<INUSD>
    
    if isempty(expdiri_selected),
      warning('No experiment selected');
      return;
    end

    addToUndoList();
    
    s = get(hObject,'String');
    vcurr = get(hObject,'Value');
    manual_pf_full = s{vcurr};
    manual_pf_new = manual_pf_full(1);
    data(expdiri_selected).manual_pf = manual_pf_new;
    %SetManualPF(data(expdiri_selected),sage_db);
    manual_pf(expdiri_selected) = manual_pf_new;

    switch lower(manual_pf_new),
      
      case 'p',
        
        % update indices
        idx_manual_p(expdiri_selected) = true;
        idx_manual_f(expdiri_selected) = false;

        % update marker
        set(h(expdiri_selected),'Marker','+');

        % update visible
        set(h(expdiri_selected),'Visible',get(hmenu.plot_manual_p,'Checked'));

      case 'f',

        % update indices
        idx_manual_p(expdiri_selected) = false;
        idx_manual_f(expdiri_selected) = true;

        % update marker
        set(h(expdiri_selected),'Marker','x');
        
        % update visible
        set(h(expdiri_selected),'Visible',get(hmenu.plot_manual_f,'Checked'));
        
      case 'u',
        
        % update indices
        idx_manual_p(expdiri_selected) = false;
        idx_manual_f(expdiri_selected) = false;
        
        % update marker
        set(h(expdiri_selected),'Marker','o');

        % update visible
        idx_visible(expdiri_selected) = true;
        set(h(expdiri_selected),'Visible','on');
        
      otherwise
        
        error('Unknown manual_pf value %s',manual_pf_new);
        
    end
    needsave = true;

    
  end

%% add note callback
 function manualpf_pushbutton_Callback(hObject,event) %#ok<INUSD>
    
   hnotes.dialog = dialog('Name','Add notes','WindowStyle','Normal','Resize','on');
   done_pos = [.29,.02,.2,.1];
   cancel_pos = [.51,.02,.2,.1];
   notes_behavioral_pos = [.02,.14,.96,.35];
   text_behavioral_pos = [.02,.49,.96,.06];
   notes_technical_pos = [.02,.57,.96,.35];
   text_technical_pos = [.02,.92,.96,.06];

   hnotes.pushbutton_done = uicontrol(hnotes.dialog,'Style','pushbutton',...
     'Units','normalized','Position',done_pos,...
     'String','Done','Callback',@notes_done_Callback);
   hnotes.pushbutton_cancel = uicontrol(hnotes.dialog,'Style','pushbutton',...
     'Units','normalized','Position',cancel_pos,...
     'String','Cancel','Callback',@notes_cancel_Callback);
   hnotes.edit_behavioral = uicontrol(hnotes.dialog,'Style','edit',...
     'Units','normalized','Position',notes_behavioral_pos,...
     'String',data(expdiri_selected).notes_behavioral,...
     'HorizontalAlignment','left',...
     'BackgroundColor','w');
   hnotes.text_behavioral = uicontrol(hnotes.dialog,'Style','text',...
     'Units','normalized','Position',text_behavioral_pos,...
     'String','Behavior notes:',...
     'HorizontalAlignment','left');
   hnotes.edit_technical = uicontrol(hnotes.dialog,'Style','edit',...
     'Units','normalized','Position',notes_technical_pos,...
     'String',data(expdiri_selected).notes_technical,...
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
        (stat(expdiri_selected,tmpi) - mu(tmpi))/z(tmpi);
      set(h(expdiri_selected),'YData',normstat(expdiri_selected,:));
      set(hselected,'YData',normstat(expdiri_selected,:));
    end

    tmpi = find(strcmpi(examinestats,'notes_technical'),1);
    if ~isempty(tmpi),
      s = data(expdiri_selected).notes_technical;
      v = ~isempty(s) && ~strcmpi(s,'None');
      stat(expdiri_selected,tmpi) = double(v);
      normstat(expdiri_selected,tmpi) = ...
        (stat(expdiri_selected,tmpi) - mu(tmpi))/z(tmpi);
      set(h(expdiri_selected),'YData',normstat(expdiri_selected,:));
      set(hselected,'YData',normstat(expdiri_selected,:));
    end
    needsave = true;
    
    close(hnotes.dialog);
    
  end

  function notes_cancel_Callback(hObject,event) %#ok<INUSD>
    
    close(hnotes.dialog);
    
  end

  function set_rest_manual_p_Callback(hObject,event) %#ok<INUSD>
    set_rest_manual_Callback('P');
  end

  function set_rest_manual_f_Callback(hObject,event) %#ok<INUSD>
    set_rest_manual_Callback('F');
  end

  function set_rest_manual_Callback(s)

    addToUndoList();
    
    b = questdlg(sprintf('Set all experiments plotted with manual_pf = U to %s?',s),...
      'Set manual_pf for rest?');
    if ~strcmpi(b,'Yes'),
      return;
    end
    
    idx = find(lower(manual_pf) == 'u');
    
    if isempty(idx),
      return;
    end
    
    needsave = true;
    
    for tmpi = idx,
      data(tmpi).manual_pf = s;
      %SetManualPF(data(tmpi),sage_db);
    end
    manual_pf(idx) = s;
    
    if lower(s) == 'p',

      % update indices
      idx_manual_p(idx) = true;
      idx_manual_f(idx) = false;
      
      % update marker
      set(h(idx),'Marker','+');
      
      % update visible
      set(h(idx),'Visible',get(hmenu.plot_manual_p,'Checked'));
      
    else

      % update indices
      idx_manual_p(idx) = false;
      idx_manual_f(idx) = true;
      
      % update marker
      set(h(idx),'Marker','x');
      
      % update visible
      set(h(idx),'Visible',get(hmenu.plot_manual_f,'Checked'));
      
    end
    
  end

  function save_Callback(hObject,event) %#ok<INUSD>
    
    [savefilename1,savepath1] = uiputfile(savefilename,'Save manual_pf tsv');
    if ~ischar(savefilename1),
      return;
    end
    savefilename = fullfile(savepath1,savefilename1);
    rc.savepath = savepath1;
    [savepath2,savefilename2] = fileparts(savefilename);
    savefilename2 = fullfile(savepath2,[savefilename2,'_diagnosticinfo.mat']);

    fid = fopen(savefilename,'w');
    fprintf(fid,'#experiment_name\tmanual_pf\tnotes_behavioral\tnotes_technical\n');
    for tmpi = 1:nexpdirs,
      fprintf(fid,'%s\t%s\t%s\t%s\n',...
        data(tmpi).experiment_name,...
        data(tmpi).manual_pf,...
        regexprep(data(tmpi).notes_behavioral,'\s+',' '),...
        regexprep(data(tmpi).notes_technical,'\s+',' '));
    end
    fclose(fid);
    save(savefilename2,'info','data','examinestats');
    needsave = false;
    
  end

  function close_fig_Callback(hObject,event)
    
    if needsave,
      answer = questdlg('Save manual_pf and notes?');
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
      struct('manual_pf',{[data.manual_pf]},...
      'notes_behavioral',{{data.notes_behavioral}},...
      'notes_technical',{{data.notes_technical}}));
    set(hmenu.undo,'Enable','on');
  end


%% reset state

  function setManualPFState(s)
  
    manual_pf = s.manual_pf;
    for tmpi = 1:nexpdirs,
      data(tmpi).manual_pf = s.manual_pf(tmpi);
      data(tmpi).notes_behavioral = s.notes_behavioral{tmpi};
      data(tmpi).notes_technical = s.notes_technical{tmpi};
    end
    idx_manual_p = lower(manual_pf) == 'p';
    idx_manual_f = lower(manual_pf) == 'f';
    idx_visible = true(1,nexpdirs);
    if strcmpi(get(hmenu.plot_manual_p,'Checked'),'off'),
      idx_visible(idx_manual_p) = false;
    end
    if strcmpi(get(hmenu.plot_manual_f,'Checked'),'off'),
      idx_visible(idx_manual_f) = false;
    end
    set(h,'Marker','o');
    set(h(idx_manual_p),'Marker','+');
    set(h(idx_manual_f),'Marker','x');
    set(h(idx_visible),'Visible','on');
    set(h(~idx_visible),'Visible','off');
    
    if ~isempty(expdiri_selected),
      manual_pf_curr = data(expdiri_selected).manual_pf;
      s = get(hmanualpf.popup,'String');
      vcurr = find(strncmpi(manual_pf_curr,s,1),1);
      if isempty(vcurr),
        error('Unknown manual_pf %s',manual_pf_curr);
      end
      set(hmanualpf.popup,'Value',vcurr);
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
 function manualpf_pushbutton_info_Callback(hObject,event) %#ok<INUSD>

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

end

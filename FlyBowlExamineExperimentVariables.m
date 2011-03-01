function handles = FlyBowlExamineExperimentVariables(varargin)

[analysis_protocol,settingsdir,datalocparamsfilestr,hfig,period,maxdatenum,figpos,datenumnow,sage_params_path,sage_db,leftovers] = ...
  myparse_nocheck(varargin,...
  'analysis_protocol','current',...
  'settingsdir','/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings',...
  'datalocparamsfilestr','dataloc_params.txt',...
  'hfig',1,...
  'period',7,...
  'maxdatenum',[],...
  'figpos',[],...
  'datenumnow',now,...
  'sage_params_path','',...
  'sage_db',[]);

if isempty(maxdatenum),
  maxdatenum = datenumnow;
end
format = 'yyyy-mm-ddTHH:MM:SS';
mindatenum = maxdatenum - period;
mindatestr = datestr(mindatenum,format);
maxdatestr = datestr(maxdatenum,format);
daterange = {mindatestr,maxdatestr};
maxdaysprev = datenumnow-mindatenum;
mindaysprev = datenumnow-maxdatenum;

%% read parameters

datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
dataloc_params = ReadParams(datalocparamsfile);
examineparamsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.examineexperimentvariablesparamsfilestr);
examine_params = ReadParams(examineparamsfile);
examinefile = fullfile(settingsdir,analysis_protocol,dataloc_params.examineexperimentvariablesfilestr);
examinestats = ReadExamineExperimentVariables(examinefile);

%% string metadata

flag_metadata_fns = {'flag_redo','flag_review'};
note_metadata_fns = {'notes_behavioral','notes_technical'};
flag_options = {{'None',''},{'Rearing problem','Flies look sick',...
  'See behavioral notes','See technical notes'}};

%% connect to SAGE for storing manualpf

if isempty(sage_db) && ~isempty(sage_params_path),
  try
    sage_db = ConnectToSage(sage_params_path);
  catch ME,
    getReport(ME);
    warning('Could not connect to Sage');
    sage_db = [];
  end
else
  sage_db = [];
end

%% get data
data_types = {'ufmf_diagnostics_summary_*','temperature_stream',...
  'ctrax_diagnostics_*','registrationdata_*','sexclassifier_diagnostics_*'};
%data_types = examinestats;

queries = leftovers;
queries(end+1:end+2) = {'daterange',daterange};
queries(end+1:end+2) = {'data_type',data_types};
queries(end+1:end+2) = {'flag_aborted',0};
queries(end+1:end+2) = {'automated_pf','P'};
data = SAGEGetBowlData(queries{:});
%load('datacache.mat','data');
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
  data(i).file_system_path = fullfile(dataloc_params.rootreaddir,expdir_bases{i});
end
expdirs = {data.file_system_path};

%% get experiment variables

nstats = numel(examinestats);
stat = nan(nexpdirs,nstats);

for i = 1:nstats,

  % special case: flags
  if ismember(examinestats{i},flag_metadata_fns),
    
    % which set of flags
    v = zeros(1,nexpdirs);
    for j = 1:numel(flag_options),
      v(ismember({data.(examinestats{i})},flag_options{j})) = j;
    end
    stat(:,i) = v;
    
  % special case: notes
  elseif ismember(examinestats{i},note_metadata_fns),
    
    v = cellfun(@(s) ~isempty(s) && ~strcmpi(s,'None'),{data.(examinestats{i})});
    stat(:,i) = double(v);
    
  % numbers
  else
    
    badidx = cellfun(@isempty,{data.(examinestats{i})});
    for j = find(badidx),
      warning('No data for stat %s, experiment %s',examinestats{i},strtrim(data(j).experiment_name));
    end
    stat(~badidx,i) = [data.(examinestats{i})];
    
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
  
  statnames{i} = examinestats{i};
  statnames{i} = strrep(statnames{i},'ufmf_diagnostics_summary','ufmf');
  statnames{i} = strrep(statnames{i},'ufmf_diagnostics_stream','ufmf');
  statnames{i} = strrep(statnames{i},'ctrax_diagnostics','ctrax');
  statnames{i} = strrep(statnames{i},'registrationdata','reg');
  statnames{i} = strrep(statnames{i},'sexclassifier_diagnostics','sex');
  
end


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
set(hfig,'Units','Pixels','Position',figpos);
hax = axes('Parent',hfig,'Units','Normalized','Position',examine_params.axespos);
% plot 0
plot(hax,[0,nstats+1],[0,0],'k-','HitTest','off');
hold(hax,'on');
colors = jet(nexpdirs)*.7;
drawnow;

%% selected experiment

hselected = plot(0,0,'o','color','k','Visible','off','HitTest','off','MarkerSize',10,'MarkerFaceColor','k');
expdiri_selected = [];
stati_selected = [];

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
dx = nstats+1;
miny = min(normstat(:));
maxy = max(normstat(:));
dy = maxy - miny;
ylim = [miny-.01*dy,maxy+.01*dy];
set(hax,'XLim',xlim,'YLim',ylim,'XTick',1:nstats,'XTickLabel',statnames,'XGrid','on');

ylabel(hax,'Stds from mean');

%% set manual_pf = p marker

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

%% options menu

hmenu = struct;
hmenu.options = uimenu('Label','Options','Parent',hfig);
hmenu.plot_manual_p = uimenu(hmenu.options,'Label','Plot manual_pf = p',...
  'Checked','on','Callback',@plot_manual_p_Callback);
hmenu.plot_manual_f = uimenu(hmenu.options,'Label','Plot manual_pf = f',...
  'Checked','on','Callback',@plot_manual_f_Callback);
hmenu.set_rest_manual_p = uimenu(hmenu.options,'Label','Set manual_pf == u -> p',...
  'Callback',@set_rest_manual_p_Callback);
hmenu.set_rest_manual_f = uimenu(hmenu.options,'Label','Set manual_pf == u -> f',...
  'Callback',@set_rest_manual_f_Callback);
hmenu.daterange = uimenu(hmenu.options,'Label',...
  sprintf('Date range: %s - %s ...',daterange{:}),...
  'Callback',@daterange_Callback);

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
w2 = 120;
w3 = 120;
h1 = 20;
h2 = 30;
c1 = ((figpos(4)-margin) + textpos_px(2))/2;
manualpf_textpos = [textpos_px(1)+textpos_px(3)+margin,c1-h1/2,w1,h1];
hmanualpf_text = uicontrol(hfig,'Style','text','Units','Pixels',...
  'Position',manualpf_textpos,'String','Manual PF:',...
  'BackgroundColor',get(hfig,'Color'),'Visible','off');
manualpf_popuppos = [manualpf_textpos(1)+manualpf_textpos(3)+margin,...
  c1-h2/2,w2,h2];
hmanualpf_popup = uicontrol(hfig,'Style','popupmenu','Units','Pixels',...
  'Position',manualpf_popuppos,'String',{'Pass','Fail','Unknown'},...
  'Value',3,'Visible','off','Callback',@manualpf_popup_Callback);
manualpf_pushbuttonpos = [manualpf_popuppos(1)+manualpf_popuppos(3)+margin,...
  c1-h2/2,w3,h2];
hmanualpf_pushbutton = uicontrol(hfig,'Style','pushbutton','Units','Pixels',...
  'Position',manualpf_pushbuttonpos,'String','Add Note...',...
  'Callback',@manualpf_pushbutton_Callback,...
  'Visible','off');   
hnotes = struct;


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
hcurrdate = uicontrol(hfig,'Style','text','Units','Pixels',...
  'Position',currpos,'String',...
  sprintf('%.1f - %.1f days ago',maxdaysprev,mindaysprev));
hnextdate = uicontrol(hfig,'Style','pushbutton','Units','Pixels',...
  'Position',nextpos,'String','>','Callback',@nextdate_Callback);
hprevdate = uicontrol(hfig,'Style','pushbutton','Units','Pixels',...
  'Position',prevpos,'String','<','Callback',@prevdate_Callback);
if mindaysprev < .0001,
  set(hnextdate,'Enable','off');
else
  set(hnextdate,'Enable','on');
end
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

%% return values

handles = struct;
handles.hx = hx;
handles.hax = hax;
handles.hfig = hfig;
handles.h = h;
handles.hselected = hselected;
handles.htext = htext;
handles.hprevdate = hprevdate;
handles.hcurrdate = hcurrdate;
handles.hnextdate = hnextdate;
% handles.hcmenu_manualpf = hcmenu_manualpf;
% handles.hcmenu = hcmenu;

%% print function

  function s = printfun(expdiri,stati)
    
    s = cell(1,2);
    s{1} = expdir_bases{expdiri};
    % special case: flags
    if ismember(examinestats{stati},flag_metadata_fns) || ...
        ismember(examinestats{stati},note_metadata_fns),
      s{2} = sprintf('%s = %s',statnames{stati},data(expdiri).(examinestats{stati}));
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
      
      if d > examine_params.maxdistclick,
        set(hselected,'Visible','off');
        stati_selected = [];
        expdiri_selected = [];
        set(htext,'String','');
        set([hmanualpf_text,hmanualpf_popup,hmanualpf_pushbutton],'Visible','off');
        return;
      end
      
      SelectionType = get(hfig,'SelectionType');
      set([hmanualpf_text,hmanualpf_popup,hmanualpf_pushbutton],'Visible','on');

      [expdiri_selected,stati_selected] = ind2sub([tmpn,nstats],closest);
      expdiri_selected = tmpidx(expdiri_selected);

      s = printfun(expdiri_selected,stati_selected);
      set(htext,'String',s);
      set(hselected,'XData',x(expdiri_selected,:),'YData',normstat(expdiri_selected,:),'Visible','on');
      if ~isempty(stati_selected),
        set(hx(stati_selected),'Color','k','FontWeight','normal');
      end
      %hxselected = stati_selected;
      set(hx(stati_selected),'Color','r','FontWeight','bold');

      manual_pf_curr = data(expdiri_selected).manual_pf;
      s = get(hmanualpf_popup,'String');
      vcurr = find(strncmpi(manual_pf_curr,s,1),1);
      if isempty(vcurr),
        error('Unknown manual_pf %s',manual_pf_curr);
      end
      set(hmanualpf_popup,'Value',vcurr);
      
      if strcmp(SelectionType,'open'),
        
        % open experiment
        if ispc,
          winopen(expdirs{expdiri});
        else
          web(expdirs{expdiri},'-browser');
        end
      end
      
    catch ME,

      fprintf('Error evaluating buttondownfcn, disabling:\n');
      set(hax,'ButtonDownFcn','');
      rethrow(ME);
    end
    
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

  function nextdate_Callback(hObject,event) %#ok<INUSD>
    
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
      'sage_params_path',sage_params_path);

  end

%% go to previous date range

  function prevdate_Callback(hObject,event) %#ok<INUSD>
    
    maxdatenum = maxdatenum - period;
    handles = FlyBowlExamineExperimentVariables(...
      'analysis_protocol',analysis_protocol,...
      'settingsdir',settingsdir,...
      'datalocparamsfilestr',datalocparamsfilestr,...
      'hfig',hfig,...
      'period',period,...
      'maxdatenum',maxdatenum,...
      'figpos',get(hfig,'Position'),...
      'datenumnow',datenumnow);

  end

%% set manual pf callback
  
  function manualpf_popup_Callback(hObject,event) %#ok<INUSD>
    
    if isempty(expdiri_selected),
      warning('No experiment selected');
      return;
    end

    s = get(hObject,'String');
    vcurr = get(hObject,'Value');
    manual_pf_full = s{vcurr};
    manual_pf_new = manual_pf_full(1);
    data(expdiri_selected).manual_pf = manual_pf_new;
    SetManualPF(data(expdiri_selected),sage_db);
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
    
    data(expdiri_selected).notes_technical = get(hnotes.edit_technical,'String');
    data(expdiri_selected).notes_behavioral = get(hnotes.edit_behavioral,'String');
    SetNotes(data(expdiri_selected),sage_db);
    
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

    b = questdlg(sprintf('Set all experiments plotted with manual_pf = U to %s?',s),...
      'Set manual_pf for rest?');
    if ~strcmpi(b,'Yes'),
      return;
    end
    
    idx = find(lower(manual_pf) == 'u');
    for tmpi = idx,
      data(tmpi).manual_pf = s;
      SetManualPF(data(tmpi),sage_db);
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

end
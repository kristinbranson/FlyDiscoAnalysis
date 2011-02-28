function handles = FlyBowlExamineExperimentVariables(varargin)

[analysis_protocol,settingsdir,datalocparamsfilestr,hfig,period,maxdatenum,figpos,datenumnow,leftovers] = ...
  myparse_nocheck(varargin,...
  'analysis_protocol','current',...
  'settingsdir','/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings',...
  'datalocparamsfilestr','dataloc_params.txt',...
  'hfig',1,...
  'period',7,...
  'maxdatenum',[],...
  'figpos',[],...
  'datenumnow',now);

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

%% get data
data_types = {'ufmf_diagnostics_summary_*','temperature_stream',...
  'ctrax_diagnostics_*','registrationdata_*','sexclassifier_diagnostics_*'};
%data_types = examinestats;

queries = leftovers;
queries(end+1:end+2) = {'daterange',daterange};
queries(end+1:end+2) = {'data_type',data_types};
queries(end+1:end+2) = {'flag_aborted',0};
queries(end+1:end+2) = {'automated_pf','P'};
%data = SAGEGetBowlData(queries{:});
load('datacache.mat','data');
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
    
    v = cellfun(@isempty,{data.(examinestats{i})});
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

% if ishandle(hfig),
%   close(hfig);
% end
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

%% gray out manual_pf = p

gray_p_color = [.7,.7,.7];
gray_f_color = [.5,.5,.5];
badidx = cellfun(@isempty,{data.manual_pf});
for j = find(badidx),
  warning('No data for manual_pf, experiment %s',strtrim(data(j).experiment_name));
end
manual_pf = repmat('U',[1,nexpdirs]);
manual_pf(~badidx) = [data.manual_pf];
idx_manual_p = lower(manual_pf) == 'p';
idx_manual_f = lower(manual_pf) == 'f';
idx_visible = true(1,nexpdirs);
set(h(idx_manual_p),'Color',gray_p_color);
set(h(idx_manual_f),'Color',gray_f_color);

%% options menu

hmenu = struct;
hmenu.options = uimenu('Label','Options','Parent',hfig);
hmenu.plot_manual_p = uimenu(hmenu.options,'Label','Plot manual_pf = p',...
  'Checked','on','Callback',@plot_manual_p_Callback);
hmenu.plot_manual_f = uimenu(hmenu.options,'Label','Plot manual_pf = f',...
  'Checked','on','Callback',@plot_manual_f_Callback);
hmenu.daterange = uimenu(hmenu.options,'Label',...
  sprintf('Date range: %s - %s ...',daterange{:}),...
  'Callback',@daterange_Callback);

%% text box

set(hax,'Units','normalized');
axpos = get(hax,'Position');
textpos = [axpos(1),axpos(2)+axpos(4),axpos(3)*.75,1-axpos(2)-axpos(4)-.01];
htext = annotation('textbox',textpos,'BackgroundColor','k','Color','g',...
  'String','Experiment info','Interpreter','none');

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

hxselected = [];
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
        if ~isempty(hxselected),
          set(hx(hxselected),'Color','k','FontWeight','normal');
        end          
        return;
      end
      
      SelectionType = get(hfig,'SelectionType');
      
      [expdiri_loc,stati_loc] = ind2sub([tmpn,nstats],closest);
      expdiri_loc = tmpidx(expdiri_loc);

      s = printfun(expdiri_loc,stati_loc);
      set(htext,'String',s);
      set(hselected,'XData',x(expdiri_loc,:),'YData',normstat(expdiri_loc,:),'Visible','on');
      if ~isempty(hxselected),
        set(hx(hxselected),'Color','k','FontWeight','normal');
      end
      hxselected = stati_loc;
      set(hx(hxselected),'Color','r','FontWeight','bold');
      
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
    else
      set(hObject,'Checked','on');
      idx_visible(idx_manual_p) = true;
    end
    
  end

%% manual_pf = f callback

  function plot_manual_f_Callback(hObject,event) %#ok<INUSD>
    
    v = get(hObject,'Checked');
    set(h(idx_manual_f),'Visible',v);
    if strcmpi(v,'on'),
      set(hObject,'Checked','off');
      idx_visible(idx_manual_f) = false;
    else
      set(hObject,'Checked','on');
      idx_visible(idx_manual_f) = true;
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
      'datenumnow',datenumnow);

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

end

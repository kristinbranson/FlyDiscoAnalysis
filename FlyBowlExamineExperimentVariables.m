function FlyBowlExamineExperimentVariables(varargin)

[analysis_protocol,settingsdir,datalocparamsfilestr,hfig,leftovers] = ...
  myparse_nocheck(varargin,...
  'analysis_protocol','current',...
  'settingsdir','/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings',...
  'datalocparamsfilestr','dataloc_params.txt',...
  'hfig',1);

%% read parameters

datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
dataloc_params = ReadParams(datalocparamsfile);
examineparamsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.examineexperimentvariablesparamsfilestr);
examine_params = ReadParams(examineparamsfile);
examinefile = fullfile(settingsdir,analysis_protocol,dataloc_params.examineexperimentvariablesfilestr);
examinestats = ReadExamineExperimentVariables(examinefile);

%% get data
data_types = {'ufmf_diagnostics_summary_*','temperature_stream',...
  'ctrax_diagnostics_*','registrationData_*','sexclassifier_daignostics_*'};
%data_types = examinestats;

queries = leftovers;
queries(end+1:end+2) = {'data_type',data_types};
data = SAGEGetBowlData(queries{:});

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
%% get experiment variables

nstats = numel(examinestats);
stat = nan(nexpdirs,nstats);

for expdiri = 1:nexpdirs,
  
  for i = 1:nstats,

    pathcurr = examinestats{i};
  
    % get examine stats for all expdirs. maybe there is a more efficient way
    % to do this?
    statcurr = data(expdiri);
    for j = 1:numel(pathcurr),
      statcurr = statcurr.(pathcurr{j});
    end
    stat(expdiri,i) = statcurr;

  end

end

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
  statnames{i} = strrep(statnames{i},'registrationData','reg');
  statnames{i} = strrep(statnames{i},'sexclassifierdiagnostics','sex');
  
end


%% z-score stats

z = sig;
z(z == 0) = 1;
normstat = bsxfun(@rdivide,bsxfun(@minus,stat,mu),z);

%% create figure

figure(hfig);
clf;
set(hfig,'Units','Pixels','Position',examine_params.figpos);
hax = axes('Parent',hfig,'Units','Normalized','Position',examine_params.axespos);
hold(hax,'on');
colors = jet(nexpdirs)*.7;

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

%% text box

axpos = get(hax,'Position');
textpos = [axpos(1),axpos(2)+axpos(4),axpos(3)/4,1-axpos(2)-axpos(4)-.01];
htext = annotation('textbox',textpos,'BackgroundColor','k','Color','g',...
  'String','Experiment info','Interpreter','none');

%% set buttondownfcn

set(hax,'ButtonDownFcn',@ButtonDownFcn);

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

%% print function

  function s = printfun(expdiri,stati)
    
    s = cell(1,2);
    s{1} = expdir_bases{expdiri};
    s{2} = sprintf('%s = %s = %s std',statnames{stati},...
      num2str(stat(expdiri,stati)),num2str(normstat(expdiri,stati)));
    
  end

%% button down function

  function ButtonDownFcn(varargin)
        
    try
    
      % get current point
      tmp = get(hax,'CurrentPoint');
      xclicked = tmp(1,1);
      yclicked = tmp(1,2);
      [d,closest] = min( ((x(:)-xclicked)/dx).^2+((normstat(:)-yclicked)/dy).^2 );
      if d > examine_params.maxdistclick,
        set(hselected,'Visible','off');
        return;
      end
      
      SelectionType = get(hfig,'SelectionType');
      
      [expdiri_loc,stati_loc] = ind2sub([nexpdirs,nstats],closest);
      s = printfun(expdiri_loc,stati_loc);
      set(htext,'String',s);
      set(hselected,'XData',x(expdiri_loc,:),'YData',normstat(expdiri_loc,:),'Visible','on');
    
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

end

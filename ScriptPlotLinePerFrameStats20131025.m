
%% parameters

analysis_protocol = '20130909';
settingsdir = 'settings';
datalocparamsfile = fullfile(settingsdir,analysis_protocol,'dataloc_params.txt');
dataloc_params = ReadParams(datalocparamsfile);
outdir = '/groups/branson/bransonlab/projects/olympiad/LineResults';
statsplotparamsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.statsplotparamsfilestr);
stats_plotparams = ReadStatsPlotParams(statsplotparamsfile);
groups = fieldnames(stats_plotparams);

%% load in data

load CollectedPrimaryPerFrameStats20131024;

%% plot stats

for linei = nlines:-1:1,
  line_name = line_names{linei};
  fprintf('Plotting stats for line %d / %d %s\n',linei,nlines,line_names{linei});
  linedir = fullfile(outdir,line_name);
  if ~exist(linedir,'dir'),
    mkdir(linedir);
  end

  for i = 1:numel(groups),
    group = groups{i};
    name = sprintf('%s, %s',line_name,group);
    if false,
      visible = 'on';
    else
      visible = 'off';
    end
    stathandles = PlotLinePerFrameStats(stats_plotparams.(group),linestats,setstats,allstatsnorm,allnframestotal,metadata,line_name,normcontrolmean,normcontrolstd,name,'visible',visible);
    if false,
      drawnow;
    end
    
    savename = sprintf('stats_%s.png',group);
    savename = fullfile(linedir,savename);
    if exist(savename,'file'),
      delete(savename);
    end
    set(stathandles.hfig,'Units','pixels','Position',stathandles.position);
    savefig(savename,stathandles.hfig,'png');

    savename = sprintf('stats_%s.pdf',group);
    savename = fullfile(linedir,savename);
    if exist(savename,'file'),
      delete(savename);
    end
    set(stathandles.hfig,'Units','pixels','Position',stathandles.position);
    savefig(savename,stathandles.hfig,'pdf');
    
  end
end

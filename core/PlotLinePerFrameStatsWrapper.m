function PlotLinePerFrameStatsWrapper(line_name,matfilename)

fprintf('Plotting stats for line %s...\n',line_name);

fprintf('Loading data from file %s...\n',matfilename);
load(matfilename);

linedir = fullfile(outdir,line_name);
if ~exist(linedir,'dir'),
  mkdir(linedir);
end

for i = 1:numel(groups), %#ok<USENS>
  
  group = groups{i};
  fprintf('Plotting group %s...\n',group);

  name = sprintf('%s, %s',line_name,group);
  stathandles = PlotLinePerFrameStats(stats_plotparams.(group),linestats,setstats,allstatsnorm,allnframestotal,metadata,line_name,normcontrolmean,normcontrolstd,name);
  savename = sprintf('stats_%s.png',group);
  savename = fullfile(linedir,savename);
  if exist(savename,'file'),
    delete(savename);
  end
  set(stathandles.hfig,'Units','pixels','Position',stathandles.position);
  savefig(savename,stathandles.hfig,'png');

  %stathandles = PlotLinePerFrameStats(stats_plotparams.(group),linestats,setstats,allstatsnorm,allnframestotal,metadata,line_name,normcontrolmean,normcontrolstd,name);
  savename = sprintf('stats_%s.pdf',group);
  savename = fullfile(linedir,savename);
  if exist(savename,'file'),
    delete(savename);
  end
  set(stathandles.hfig,'Units','pixels','Position',stathandles.position);%,'Renderer','painters');
  savefig(savename,stathandles.hfig,'pdf');
  
end

if isdeployed,
  delete(findall(0,'type','figure'));
end
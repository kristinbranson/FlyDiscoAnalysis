function newexpdir = CopyStatsPlots(expdir,rootdir)

[~,expname] = fileparts(expdir);
newexpdir = fullfile(rootdir,expname);
if ~exist(newexpdir,'dir'),
  mkdir(newexpdir);
end
copyfile(fullfile(expdir,'stats.html'),fullfile(newexpdir,'stats.html'),'f');
unix(sprintf('cp -rf %s/analysis_plots %s/analysis_plots',expdir,newexpdir));
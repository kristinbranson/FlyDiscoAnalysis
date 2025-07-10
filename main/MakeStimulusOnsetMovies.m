function [outfiles,ions,fliesplot] = MakeStimulusOnsetMovies(moviefile,trx,outfilebase,varargin)

[nfliesplot,fliesplot,ions,nperiodsplot,stimsets,visible,force,leftovers] = ...
  myparse_nocheck(varargin,'nfliesplot',[],...
  'fliesplot',[],'ions',[],'nperiodsplot',[],...
  'stimsets',cell(0,2),...
  'visible','on','force',false);

if ~isempty(fliesplot),
  nfliesplot = numel(fliesplot);
end
if isempty(nfliesplot),
  nfliesplot = trx.nflies;
end

ind = trx.getIndicatorLED(1);
non = numel(ind.starton);

if ~ischar(moviefile),
  readframe = moviefile;
  fid = [];
else
  [readframe,~,fid] = get_readframe_fcn(moviefile);
end

hfig = figure('Visible',visible);
hax = axes('Position',[0,0,1,1],'Parent',hfig);

% if isempty(fliesplot),
%   if nfliesplot < trx.nflies,
%     fliesplot = round(linspace(1,trx.nflies,nfliesplot));
%   else
%     fliesplot = 1:trx.nflies;
%   end
% end
if isempty(ions),
  if isempty(nperiodsplot) || nperiodsplot < 0,
    nperiodsplot = non;
  end
  if isempty(stimsets),
    ions = unique(round(linspace(1,non,nperiodsplot)));
  else
    ions = [];
    for seti = 1:size(stimsets,1),
      idxcurr = unique(round(linspace(1,numel(stimsets{seti,2}),nperiodsplot)));
      ions = [ions,stimsets{seti,2}(idxcurr)]; %#ok<AGROW>
    end
  end
end
nperiodsplot = numel(ions);

if isempty(fliesplot),
  if nfliesplot < trx.nflies,
    fliesplot = ChooseFliesPlot(trx,ind,ions,nfliesplot);
    fliesplot = sort(fliesplot);
    nfliesplot = numel(fliesplot);
  else
    fliesplot = 1:trx.nflies;
  end
end
outfiles = cell(numel(ions),numel(fliesplot));
for ioni = 1:numel(ions),
  ion = ions(ioni);
  for fliesploti = 1:numel(fliesplot),
    fly = fliesplot(fliesploti);
    outfile = sprintf('%s_period%02d_fly%02d.gif',outfilebase,ion,fly);
    outfiles{ioni,fliesploti} = outfile;
    if ~force && exist(outfile,'file'),
      continue;
    end
    MakeStimulusOnsetMovie(readframe,trx,outfile,ion,fly,ind,'hax',hax,leftovers{:});
  end
end

close(hfig);
if ~isempty(fid) && fid > 1,
  fclose(fid);
end
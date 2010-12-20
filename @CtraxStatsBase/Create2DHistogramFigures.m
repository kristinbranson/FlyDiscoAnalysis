% [hfig,hax,nfigs,figi_main,figi_perexp,figi_perfly,...
% nax,axi_hist,axi_pluserror,axi_minuserror,axi_perexp,axi_perfly] = ...
% CtraxStatsBase.Create2DHistogramFigures(hax,hfig,nflies,nexpdirs,...
% doplotperexp,doplotperfly,doploterrorbars,figpos)
% Returns handles to figures, axes required for plotting 2D histograms.
% It uses axes in HAX, HFIG as possible, otherwise it creates new axes,
% figures based on number required by NFLIES, NEXPDIRS, DOPLOTPEREXP,
% DOPLOTPERFLY, DOPLOTERRORBARS.

function [hfig,hax,nfigs,figi_main,figi_perexp,figi_perfly,...
  nax,axi_hist,axi_pluserror,axi_minuserror,axi_perexp,axi_perfly] = ...  
  Create2DHistogramFigures(hax,hfig,nflies,nexpdirs,...
  doplotperexp,doplotperfly,doploterrorbars,figpos,docla)

%% initialize in case not set
figi_perexp = [];
figi_perfly = [];
axi_pluserror =[];
axi_minuserror = [];
axi_perexp = [];
axi_perfly = [];

%% index within hfig
nfigs = 1;
figi_main = nfigs;
if doplotperexp,
  nfigs = nfigs + 1;
  figi_perexp = nfigs;
end
if doplotperfly,
  nfigs = nfigs + 1;
  figi_perfly = nfigs;
end
  
%% index within hax
nax = 1;
axi_hist = nax;
if doploterrorbars,
  nax = nax+1;
  axi_pluserror = nax;
  nax = nax+1;
  axi_minuserror = nax;
end
if doploterrorbars,
  nax_main = 3;
  axi_main = [axi_hist,axi_pluserror,axi_minuserror];
else
  nax_main = 1;
  axi_main = axi_hist;
end
if doplotperexp,
  axi_perexp = nax + (1:nexpdirs);
  nax = nax + nexpdirs;
end
if doplotperfly,
  axi_perfly = nax + (1:nflies);
  nax = nax + nflies;
end

%% get an axis for the main histogram, error hists

[hax,hfig] = get_axes(hax,hfig,...
  'figpos',figpos,...
  'axparams',{1,nax_main,.05},...
  'axi_create',axi_main,...
  'figi',figi_main);
if nax_main > 1,
  hax([1,2,3]) = hax([2,3,1]);
end

%% get axes for exps
if doplotperexp,
  
  figpos_perexp = get(hfig(figi_main),'Position');
  figpos_perexp(4) = figpos_perexp(4)*nax_main/nexpdirs;

  [hax,hfig] = get_axes(hax,hfig,...
    'figpos',figpos_perexp,...
    'axparams',{1,nexpdirs,[[.01,.01];[.1,.1]]},...
    'axi_create',axi_perexp,...
    'figi',figi_perexp);
end

%% get axes for flies
if doplotperfly,
  figpos_perfly = get(hfig(figi_main),'Position');
  figpos_perfly(4) = figpos_perfly(4)*nax_main/nflies;
  
  [hax,hfig] = get_axes(hax,hfig,...
    'figpos',figpos_perfly,...
    'axparams',{1,nflies,[[.01,.01];[.1,.1]]},...
    'axi_create',axi_perfly,...
    'figi',figi_perfly);
end

%% clear the axes
if docla,
  for axi = 1:nax,
    cla(hax(axi));
  end
end

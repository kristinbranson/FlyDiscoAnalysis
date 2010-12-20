% [edges,nbins,centers] = obj.SelectHistEdges(fn,flies,conditions,edges,nbins,...
% lim,lim_prctile,outputfun,binmode)
% Selects the edges based on all the information given in the input
% arguments.
% If EDGES is not empty, then EDGES is used.
% Otherwise, it sets the lowest edge to LIM(1) if LIM(1) is not NaN.
% Otherwise, it sets the lowest edge to the LIM_PRCTILE(1) of all data
% corresponding to field FN, flies FLIES, conditions CONDITIONS, and
% output function OUTPUTFUN.
% Similarly, the highest edge is set to LIM(1) if LIM(1) is not NaN.
% Otherwise, it sets the highest edge to the LIM_PRCTILE(2) of all data
% corresponding to field FN, flies FLIES, conditions CONDITIONS, and
% output function OUTPUTFUN.
% It selects edges between these bounds according to BINMODE, which is
% either 'linear', 'log' (log binning), or 'logabs' (log binning in
% both directions from zero).

function [edges,nbins,centers] = SelectHistEdges(obj,...
  fn,edges,flies,conditions,nbins,lim,lim_prctile,outputfun,binmode)


% set edges if not input
if ~isempty(edges),
  nbins = length(edges)-1;
  centers = (edges(1:end-1)+edges(2:end))/2;
  return;
end
if any(isnan(lim)),
  minx = inf;
  maxx = -inf;
  % for now, don't worry about memory, store all the data
  useprctile = any(~isnan(lim_prctile));
  if useprctile,
    xall = [];
  end
  for my_fly = flies,
    if ~isempty(conditions),
      my_x = obj.trx(my_fly).(fn)(conditions(obj.trx(my_fly)));
    else
      my_x = obj.trx(my_fly).(fn);
    end
    if ~isempty(outputfun),
      my_x = outputfun(my_x);
    end
    badidx = isinf(my_x) | isnan(my_x);
    minx = min(minx,min(my_x(~badidx)));
    maxx = max(maxx,max(my_x(~badidx)));
    if useprctile,
      xall = [xall,my_x]; %#ok<AGROW>
    end
  end
  if isnan(lim(1)),
    if isnan(lim_prctile(1)),
      lim(1) = minx;
    else
      lim(1) = prctile(xall,lim_prctile(1));
    end
  end
  if isnan(lim(2)),
    if isnan(lim_prctile(2)),
      lim(2) = maxx;
    else
      lim(2) = prctile(xall,lim_prctile(2));
    end
  end
end
switch lower(binmode),
  case 'linear',
    edges = linspace(lim(1),lim(2),nbins+1);
  case 'log',
    % make sure > 0

      if lim(1) > 0,
        tmplim = lim;
      else
        off = (lim(2)-lim(1))/nbins;
        tmplim = lim - lim(1) + off;
      end
      edges = exp(linspace(log(tmplim(1)),log(tmplim(2)),nbins+1));
      if lim(1) < 0,
        edges = edges + lim(1) - off;
      end
      
  case 'logabs',
    % from lim(1) to 0, we do log spacing on neg value
    % from 0 to lim(2), we do log spacing on pos value
    
    if lim(1) > 0,
      % corner case: both are positive
      tmplim = lim;
      edges = exp(linspace(log(tmplim(1)),log(tmplim(2)),nbins+1));
      
    elseif lim(2) < 0,
      % corner case: both are negative
      tmplim = -lim;
      edges = exp(linspace(log(tmplim(1)),log(tmplim(2)),nbins+1));
      edges = fliplr(-edges);
      
    else
      % how much of data is below 0
      fracneg = -lim(1) / (lim(2)-lim(1));
      fracpos = 1 - fracneg;
      % how many bins will we have on one side of 0
      if fracneg > .5,
        nbinsneg = floor(fracneg*nbins);
        nbinspos = nbins - nbinsneg;
      else
        nbinspos = floor(fracpos*nbins);
        nbinsneg = nbins - nbinspos;
      end
      % positive edges
      off = (lim(2)-lim(1))/nbins;
      tmplim = [0,lim(2)] + off;
      edgespos = exp(linspace(log(tmplim(1)),log(tmplim(2)),nbinspos+1));
      edgespos = edgespos - off;
      % negative edges
      tmplim = [0,-lim(1)];
      tmplim = tmplim + off;
      edgesneg = exp(linspace(log(tmplim(1)),log(tmplim(2)),nbinsneg+1));
      edgesneg = fliplr(-(edgesneg - off));
      edges = [edgesneg(1:end-1),edgespos];
    end
    
  otherwise,
    error('Unknown binmode %s',binmode);
end

centers = (edges(1:end-1)+edges(2:end))/2;
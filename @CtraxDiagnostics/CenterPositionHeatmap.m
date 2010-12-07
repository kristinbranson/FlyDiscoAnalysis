% heatmap = obj.CenterPositionHeatmap(varargin)

function heatmap = CenterPositionHeatmap(obj,varargin)

% parse inputs
nbins = 100;
[edges_x,edges_y,nbins_x,nbins_y,nbins,xlim,ylim,hfig,doplot,conditions,...
  flies,expdirs,jackknife] = ...
  myparse(varargin,'edges_x',[],'edges_y',[],'nbins_x',[],'nbins_y',[],'nbins',nbins,...
  'xlim',nan(1,2),'ylim',nan(1,2),'hfig',[],'doplot',false,'conditions',[],...
  'flies',1:obj.nflies,'expdirs',obj.expdir_bases,'jacknife','perfly');

%% set edges if not input
if isempty(edges_x),
  if isempty(nbins_x),
    nbins_x = nbins;
  end
  % left-most point in the arena
  if isnan(xlim(1)),
    xlim(1) = obj.arena_center(1)-obj.arena_radius;
  end
  % right-most point in the arena
  if isnan(xlim(2)),
    xlim(2) = obj.arena_center(1)+obj.arena_radius;
  end
  edges_x = linspace(xlim(1),xlim(2),nbins_x+1);
else
  nbins_x = length(edges_x) - 1;  
end
if isempty(edges_y),
  if isempty(nbins_y),
    nbins_y = nbins;
  end
  % top-most point in the arena
  if isnan(ylim(1)),
    ylim(1) = obj.arena_center(2)-obj.arena_radius;
  end
  % bottom-most point in the arena
  if isnan(ylim(2)),
    ylim(2) = obj.arena_center+obj.arena_radius;
  end
  edges_y = linspace(ylim(1),ylim(2),nbins_y+1);
else
  nbins_y = length(edges_y) - 1;  
end
centers_x = (edges_x(1:end-1)+edges_x(2:end))/2;
centers_y = (edges_y(1:end-1)+edges_y(2:end))/2;

%% take intersection of flies and expdir_bases
[didfind,ns] = ismember(expdirs,obj.expdir_bases);
if ~all(didfind),
  warning(['The following expdirs are not loaded: ',sprintf('%s ',expdirs{~didfind})]);
end
ns = ns(didfind);

% which flies does this correspond to?
allflies_perexp = [obj.movie2flies{ns}];
flies = intersect(flies,allflies_perexp);
ns = unique(obj.fly2movie(flies));
nflies = length(flies);
% save for returning
heatmap.ns = ns;
heatmap.nflies = nflies;

% check that there is at least one experiment
nexpdirs = length(ns);
if nexpdirs == 0,
  error('No experiments selected.');
end


%% set jackknife method if not set
if isempty(jackknife),
  if nexpdirs > 1,
    jackknife = 'perexp';
  elseif nflies == 1,
    jackknife = 'none';
  else
    jackknife = 'perfly';
  end
end

% to split over experiments, we need more than one experiment
if nexpdirs == 1 && strcmpi(jackknife,'perexp'),
  warning('Only one experiment selected, but jackknife = ''perexp''. Splitting over flies instead.');
  jackknife = 'perfly';
end
% to split over flies, we need more than one fly
if nflies == 1 && strcmpi(jackknife,'perfly'),
  warning('Only one fly selected, but jackknife = ''perfly''. Splitting over frames not implemented. Not jackknifing.');
  jackknife = 'none';
end
res.jackknife = jackknife;


%% histogram for each fly
heatmap.countsperfly = nan([nbins_y,nbins_x,nflies]);
for i = 1:nflies,

  fly = flies(i);

  % filter out frames where conditions apply
  if ~isempty(conditions),
    idx = conditions(obj.trx(fly));
  else
    idx = true(size(obj.trx(fly).x_mm));
  end
  x = obj.trx(fly).x_mm(idx);
  y = obj.trx(fly).y_mm(idx);
  
  % histogram positions for the current fly
  countscurr = hist3([y(:),x(:)]',{edges_y,edges_x});
  countscurr(:,end-1) = countscurr(:,end-1)+countscurr(:,end);
  countscurr(end-1,:) = countscurr(end-1,:)+countscurr(end,:);
  countscurr = countscurr(1:end-1,1:end-1);
  
end

%% set up to index into data structures only computed for selected flies, selected exps

% flyidx i corresponds to fly flies(i)
fly2idx = sparse(ones(1,nflies),flies,1:nflies,1,obj.nflies);
% expidx i corresponds to experiment ns(i)
n2idx = sparse(ones(1,nexpdirs),ns,1:nexpdirs,1,obj.nexpdirs);

expdirs = obj.expdir_bases(ns);
movie2flies = cell(1,nexpdirs);
nfliespermovie = zeros(1,nexpdirs);
for i = 1:length(ns),
  n = ns(i);
  movie2flies{i} = full(fly2idx(obj.movie2flies{n}));
  movie2flies{i}(movie2flies{i} == 0) = [];
  nfliespermovie(i) = length(movie2flies{i});
end
fly2movie = full(n2idx(obj.fly2movie(flies)));

%% per-fly fractions

% these won't all be used for all averaging conditions, but it should be
% quick and will speed up jackknifing, so compute them
Zperfly = sum(sum(countsperfly,1),2);
fracperfly = bsxfun(@rdivide,countsperfly,Zperfly);

%% compute frac using all the data

% compute histogram using all the data
frac = compute_frac(1:nflies,1:nexpdirs);


% compute mean, std over flies
heatmap.meanfracperfly = nanmean(heatmap.fracperfly,3);
heatmap.stdfracperfly = nanstd(heatmap.fracperfly,1,3);

% normalize counts over all flies
heatmap.fracallflies = heatmap.counts / sum(heatmap.n);

% store the centers and edges
heatmap.centers_x = centers_x;
heatmap.centers_y = centers_y;
heatmap.edges_x = edges_x;
heatmap.edges_y = edges_y;

if doplot,
  doresize = isempty(hfig) || ~ishandle(hfig);
  if isempty(hfig),
    hfig = figure;
  else
    figure(hfig);
    clf;
  end
  if doresize,
    set(hfig,'Position',[1,1,1000,800]);
  end
  heatmap.hfig = hfig;
  heatmap.hax = createsubplots(2,2,.05);
  axes(heatmap.hax(1)); %#ok<*MAXES>
  imagesc(heatmap.fracallflies);
  axis image;
  title('Fraction of time spent in each bin for all flies combined');
  colorbar;
  axes(heatmap.hax(2));
  imagesc(heatmap.meanfracperfly);
  axis image;
  title('Mean fraction of time spent in each bin per fly');
  colorbar;
  axes(heatmap.hax(3));
  imagesc(log(heatmap.fracallflies));
  axis image;
  title('Log-fraction of time spent in each bin for all flies combined');
  colorbar;
  axes(heatmap.hax(4));
  imagesc(log(heatmap.meanfracperfly));
  axis image;
  title('Log-mean fraction of time spent in each bin per fly');
  colorbar;

  linkaxes(heatmap.hax);
  
end

function my_frac = compute_frac(my_flies,my_ns)

  % take the intersection of the specified flies and experiments
  my_allflies_perexp = [movie2flies{my_ns}];
  my_flies = intersect(my_flies,my_allflies_perexp);
  my_ns = unique(fly2movie(my_flies));
  
  switch lower(averaging),
    case 'allexps_allflies',
      
      % treat all frames of data the same: just add up all the counts
      my_counts = sum(countsperfly(:,:,my_flies),3);
      my_frac = my_counts / sum(my_counts(:));
      
    case 'allexps_perfly',
      
      % get per-fly fracs, but treat all flies in all experiments the same
      my_frac = nanmean(fracperfly(:,:,my_flies),3);
      
    case 'perexp_allflies',
      
      % get per-exp fracs, but treat all frames within an experiment the
      % same
      
      % loop over experiments
      my_frac = zeros(nbins_y,nbins_x);
      my_nexps = zeros(nbins_y,nbins_x);
      for my_n = my_ns(:)',
        % which flies in this experiment
        my_flies_curr = intersect(my_flies,movie2flies{my_n});
        % frac for these flies
        my_frac_curr = sum(countsperfly(:,:,my_flies_curr),3);
        my_frac_curr = my_frac_curr / sum(my_frac_curr(:));
        % only include good data
        my_isdata = ~isnan(my_frac_curr);
        my_frac(my_isdata) = my_frac(my_isdata) + my_frac_curr(my_isdata);
        my_nexps(my_isdata) = my_nexps(my_isdata) + 1;
      end
      % normalize
      my_frac = my_frac ./ my_nexps;
      
    case 'perexp_perfly',
      
      % loop over experiments
      my_frac = zeros(1,nbins);
      my_nexps = zeros(1,nbins);
      for my_n = my_ns(:)',
        
        % which flies in this experiment
        my_flies_curr = intersect(my_flies,movie2flies{my_n});
        
        % average histogram for this experiment
        my_frac_curr = nanmean(fracperfly(:,:,my_flies_curr),3);
        
        % only include good data
        my_isdata = ~isnan(my_frac_curr);
        my_frac(my_isdata) = my_frac(my_isdata) + my_frac_curr(my_isdata);
        my_nexps(my_isdata) = my_nexps(my_isdata) + 1;
        
      end
      % normalize
      my_frac = my_frac ./ my_nexps;
      
    otherwise
      error('Unknown averaging method %s',averaging);
  end
  
end

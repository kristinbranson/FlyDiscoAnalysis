%% ScriptBehaviorNormalization

%% set up path


addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
addpath /groups/branson/bransonlab/projects/olympiad/SAGE/MATLABInterface/Trunk;
settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
rootdatadir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data';
analysis_protocol = 'current';

collecteddatamatfile = '/groups/branson/home/bransonk/behavioranalysis/code/Jdetect/Jdetect/experiments/gal4screen/ExperimentsFracTimePerSet20120906T160000.mat';

%% load in data

% fractime is nexps x nbehaviors 
% metadata is 1 x nexps
% mean_fractime_perset is nsets x nbehaviors
% setmetadata is 1 nsets
% TODO: more behavior measures
load(collecteddatamatfile);
nbehaviors = size(fractime,2);

%% parameters

normalize_datatype = 'control';
normalize_unit = 'experiment';
nbins_pairwise_plot = 50;
timestamp = datestr(now,'yyyymmdd');

savefilestr = sprintf('BehaviorNormalization_%s_%s_%s',normalize_datatype,normalize_unit,timestamp);

%% get relevant data, depending on type of normalization

% normalize by control data
if strcmpi(normalize_datatype,'control'),
  
  if strcmpi(normalize_unit,'set'),
    
    idx = strcmp(linename_perset,'pBDPGAL4U');
    X = mean_fractime_perset(idx,:);
    
  else
    
    idx = strcmp({metadata.line_name},'pBDPGAL4U');
    X = fractime(idx,:);
    
  end
  
% normalize by all GMRs
elseif strcmpi(normalize_datatype,'gmr'),
  
  if strcmpi(normalize_unit,'line'),
    
    [linenames,~,lineidx] = unique(linename_perset);
    nlines = numel(linenames);
    mean_fractime_perline = nan(nlines,nbehaviors);
    for linei = 1:nlines,
      mean_fractime_perline(linei,:) = nanmean(mean_fractime_perset(lineidx==linei,:),1);
    end
    X = mean_fractime_perline;
    
  elseif strcmpi(normalize_unit,'set'),
    
    X = mean_fractime_perset;
    
  elseif strcmpi(normalize_unit,'experiment'),
    
    X = fractime;
    
  end
  
end

%% all pairs 

% normalize yi by xi
% pairwise_norm_fcns(xi,yi) is how to normalize measure yi wrt xi
pairwise_norm_fcns = struct('fcn',cell(nbehaviors,nbehaviors),...
  'coeffs',cell(nbehaviors,nbehaviors));
for xi = 1:nbehaviors,
  for yi = 1:nbehaviors,
    if xi == yi,
      continue;
    end
    [fcn_curr,coeffs_curr] = GetBehaviorNormalizationFunction(X(:,xi), X(:,yi));
    pairwise_norm_fcns(xi,yi).fcn = fcn_curr;
    pairwise_norm_fcns(xi,yi).coeffs = coeffs_curr;
  end
end

%% plot 

lims = cell(1,nbehaviors);
edges = cell(1,nbehaviors);
centers = cell(1,nbehaviors);
max_lim_orders = 2;
for i = 1:nbehaviors,
  lims{i} = [min(X(:,i)),prctile(X(:,i),99)];
  limcurr = lims{i};
  limcurr(1) = max(limcurr(1),limcurr(2)*10^(-max_lim_orders));
  [edges{i},centers{i}] = SelectHistEdges(nbins_pairwise_plot,limcurr,'log');
%   edges{i} = linspace(lims{i}(1),lims{i}(2),nbins_pairwise_plot+1);
%   centers{i} = (edges{i}(1:end-1)+edges{i}(2:end))/2;
  %edges{i}(end) = edges{i}(end) + diff(lims{i})*.01;
  edges{i}(1) = min(X(:,i));
  %edges{i}(end) = max(X(:,i));
end  

hfig = 1;
figure(hfig);
clf;
hax = createsubplots(nbehaviors,nbehaviors,[.02,.01;.02,.01]);


maxcounts = 0;
%him = [];
for xi = 1:nbehaviors,
  xlims = edges{xi}([1,nbins_pairwise_plot+1]);
  %xlims = [-1,1]*diff(xlims)*.025 + xlims;
  for yi = 1:nbehaviors,
    axi = sub2ind([nbehaviors,nbehaviors],yi,xi);
    
    ylims = edges{yi}([1,nbins_pairwise_plot+1]);
    %ylims = [-1,1]*diff(ylims)*.025 + ylims;

    if xi == yi,
      text(mean(xlims),mean(xlims),scoresfilestr{xi,1},...
        'HorizontalAlignment','center','Interpreter','none',...
        'Parent',hax(axi));
      axis(hax(axi),'off');
    else
      counts = hist3(X(:,[yi,xi]),'Edges',{edges{yi},edges{xi}});
      counts(end-1,:) = counts(end-1,:)+counts(end,:);
      counts(:,end-1) = counts(:,end-1)+counts(:,end);
      counts = counts(1:end-1,1:end-1) / sum(sum(counts(1:end-1,1:end-1)));
      maxcounts = max(maxcounts,max(counts(:)));
      [xgrid00,ygrid00] = meshgrid(edges{xi}(1:end-1),edges{yi}(1:end-1));
      [xgrid10,ygrid10] = meshgrid(edges{xi}(2:end),edges{yi}(1:end-1));
      [xgrid01,ygrid01] = meshgrid(edges{xi}(1:end-1),edges{yi}(2:end));
      [xgrid11,ygrid11] = meshgrid(edges{xi}(2:end),edges{yi}(2:end));
      patch([xgrid00(:)';xgrid10(:)';xgrid11(:)';xgrid01(:)';xgrid00(:)'],...
        [ygrid00(:)';ygrid10(:)';ygrid11(:)';ygrid01(:)';ygrid00(:)'],counts(:)',...
        'Parent',hax(axi),'EdgeColor','none');
      hold(hax(axi),'on');
      ycurr = [ones(2,1),xlims(:)]*pairwise_norm_fcns(xi,yi).coeffs;
      plot(hax(axi),xlims,ycurr','w-');

%      imagesc(xlims,ylims,counts,'Parent',hax(axi));
%       plot(hax(axi),X(:,xi),X(:,yi),'k.');
    end
    axis(hax(axi),[xlims,ylims],'xy');
    if yi < nbehaviors && yi ~= xi-1,
      set(hax(axi),'XTickLabel',{});
    end
    if xi > 1 && xi ~= yi+1,
      set(hax(axi),'YTickLabel',{});
    end
  end
end

cminterp = logscale_colormap(jet(256),[0,1]);
colormap(cminterp);

savefig([savefilestr,'.pdf'],hfig,'pdf');

%% 
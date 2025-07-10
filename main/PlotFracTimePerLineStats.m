function hfigs = PlotFracTimePerLineStats(stats,varargin)

nlines = numel(stats.line_names);
nbehaviors = numel(stats.behaviornames);

% parameters
npadx = 1;
prctiles_plot_control = [1,5,25,50,75,95,99];
orderby = 'Chase';
figoffset = 1000;

[orderby,plotexpdata] = myparse(varargin,'orderby',orderby,'plotexpdata',true);

% compute control prctiles
idxcontrol = strcmp(stats.linename_perset,'pBDPGAL4U');
mucontrol = mean(stats.mean_fractime_perset(idxcontrol,:),1); %#ok<NASGU>
stdcontrol = std(stats.mean_fractime_perset(idxcontrol,:),1,1); %#ok<NASGU>
prctiles_control = prctile(stats.mean_fractime_perset(idxcontrol,:),prctiles_plot_control,1);
medi = find(prctiles_plot_control == 50); %#ok<NASGU>

npregions = (numel(prctiles_plot_control)-1)/2;

% choose colors 
linecolors = jet(nlines)*.7;
linecolors(strcmpi(stats.line_names,'pBDPGAL4U'),:) = 0;

% timestamps
%set_exp_datetime = regexp(stats.set_names,'\d{8}T\d{6}$','match','once');
% is this old data?
%set_isprelimdata = datenum(set_exp_datetime,'yyyymmddTHHMMSS')' <= datenum(mindatemodern,'yyyymmdd');

% order all lines by fraction of time performing a behavior


behaviori = find(strcmpi(orderby,stats.behaviornames),1);
if ~isempty(behaviori),
  [~,order] = sort(stats.mean_fractime_perline(:,behaviori));
  [~,reorder] = sort(order);
else
  order = 1:nlines;
  reorder = 1:nlines;
end
[~,lineidx] = ismember(stats.linename_perset,stats.line_names);
[~,explineidx] = ismember(stats.linename_perexp,stats.line_names);

xticklabel = stats.line_names(order);
xticklabel = regexprep(xticklabel,'GMR_','');
xticklabel = regexprep(xticklabel,'^MB','');
xticklabel = regexprep(xticklabel,'_AE_01$','');
xticklabel = regexprep(xticklabel,'_01$','');

% loop through behaviors
hfigs = nan(1,nbehaviors);

hleg = [];
sleg = {};
isexpleg = false;
issetleg = false;
islineleg = false;

for behaviori = 1:nbehaviors,

  % create figure
  hfigs(behaviori) = figoffset+behaviori;
  figure(hfigs(behaviori));
  clf;
  set(hfigs(behaviori),'Units','pixels','Position',[110,1090,1600,700]);
  hax = axes('Position',[.05,.1,.9,.85]); %#ok<LAXES,NASGU>

  % draw control distribution
  for j = npregions:-1:1
    ycurr = prctiles_control([npregions+1-j,npregions+1+j],behaviori);
    patch([0,nlines+1,nlines+npadx,1-npadx,1-npadx],ycurr([1,1,2,2,1]),diff(5+prctiles_plot_control([npregions+1-j,npregions+1+j]))/105*[1,1,1],'LineStyle','none');
    hold on;
  end
  plot([1-npadx,nlines+npadx],prctiles_control(npregions+1,behaviori)+[0,0],'k-');
  for j = 1:numel(prctiles_plot_control),
    text(1,prctiles_control(j,behaviori),sprintf('%d%%ile',prctiles_plot_control(j)),'HorizontalAlignment','left');
  end
    
  for j = 1:nlines,
    idxcurr = lineidx == j;
    expidxcurr = explineidx == j;
    plot(reorder([j,j]),stats.mean_fractime_perline(j,behaviori)+[-1,1]*stats.std_fractime_perline(j,behaviori),'-','color',(linecolors(j,:)+1)/2);
    if strcmpi('pBDPGAL4U',stats.line_names{j}),
      setcolor = (linecolors(j,:)+3)/4;
    else
      setcolor = (linecolors(j,:)+1)/2;
    end
    if plotexpdata,
      hcurr = plot(reorder(j)+zeros(size(stats.mean_fractime_perexp(expidxcurr,behaviori))),stats.mean_fractime_perexp(expidxcurr,behaviori),'.','color',(linecolors(j,:)+3)/4);
      if ~isexpleg,
        hleg(end+1) = hcurr;
        sleg{end+1} = 'Experiment';
        isexpleg = true;
      end
    end
    if any(idxcurr),
      hcurr = plot(reorder(j)+zeros(size(stats.mean_fractime_perset(idxcurr,behaviori))),stats.mean_fractime_perset(idxcurr,behaviori),'o','color',setcolor,'markerfacecolor',setcolor,'markersize',4);
      if ~issetleg,
        hleg(end+1) = hcurr;
        sleg{end+1} = 'Set';
        issetleg = true;
      end
    end
    hcurr = plot(reorder(j),stats.mean_fractime_perline(j,behaviori),'s','color',linecolors(j,:),'markerfacecolor',linecolors(j,:));
    if ~islineleg,
      hleg(end+1) = hcurr;
      sleg{end+1} = 'Line';
      islineleg = true;
    end

  end  
  
  set(gca,'XLim',[1-npadx,nlines+npadx],'YLim',[0,max(stats.mean_fractime_perset(:,behaviori))*1.01]);

  %xlabel('Line');
  ylabel(sprintf('Frac time %s',stats.behaviornames{behaviori}),'Interpreter','none');

%   isleft = true;
%   
%   for jj = nlines-9:nlines,
%     j = order(jj);
%     shortname = regexprep(stats.line_names{j},'GMR_','');
%     shortname = regexprep(shortname,'_AE_01$','');
%     shortname = regexprep(shortname,'_01$','');
%     if isleft,
%       text(nlines+npadx,stats.mean_fractime_perline(j,behaviori),shortname,'color',linecolors(j,:),'HorizontalAlignment','left','Interpreter','none');
%     else
%       text(jj,stats.mean_fractime_perline(j,behaviori),shortname,'color',linecolors(j,:),'HorizontalAlignment','right','Interpreter','none');
%     end
%     isleft = ~isleft;
%   end

  set(gca,'XTick',1:nlines,'XTickLabel',xticklabel);

  rotateticklabel(gca,90);
  
  if behaviori == 1,
    legend(hleg,sleg,'Location','NorthWest');
  end
      
end
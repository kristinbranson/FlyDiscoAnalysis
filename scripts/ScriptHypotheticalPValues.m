pas = [.005,.01,.02,.05,.1];
pbgivenas = .1:.1:1;
pagivenbs = .1:.05:1;

pval = nan(numel(pas),numel(pbgivenas),numel(pagivenbs));
for pai = 1:numel(pas),
  pa = pas(pai);
  for pbgivenai = 1:numel(pbgivenas),
    pbgivena = pbgivenas(pbgivenai);
    for pagivenbi = 1:numel(pagivenbs),
      pagivenb = pagivenbs(pagivenbi);
      
      pb = pbgivena*pa/pagivenb;
      nb = round(pb*nlines);
      na = round(nb*pagivenb);
        
      if na <= 0,
        pval(pai,pbgivenai,pagivenbi) = 1;
      else
        pval(pai,pbgivenai,pagivenbi) = 1-binocdf(na-1,nb,pa);
      end
      
    end
  end
end
%%

hfig = 1435;
figure(hfig);
clf;
hold on;
colors = jet(numel(pbgivenas))*.7;
nc = ceil(sqrt(numel(pas)));
nr = ceil(numel(pas)/nc);
hax = createsubplots(nr,nc,[.05,.1]);
hax = reshape(hax,[nr,nc])';
if nr*nc > numel(pas),
  delete(hax(numel(pas)+1:end));
end

xlim = [pagivenbs(1),pagivenbs(end)];
xlim = xlim + diff(xlim)*.05*[-1,1];

for pai = 1:numel(pas),
  pa = pas(pai);
  hold(hax(pai),'on');
  h = nan(1,numel(pbgivenas));
%   for pv = [.01,.001,.0001],
%     plot(hax(pai),xlim,pv+[0,0],'--','Color',[.5,.5,.5]);
%   end
  for pbgivenai = 1:numel(pbgivenas),
    pbgivena = pbgivenas(pbgivenai);
    idxcurr = 1:numel(pagivenbs);
    h(pbgivenai) = plot(hax(pai),pagivenbs(idxcurr),max(eps,squeeze(pval(pai,pbgivenai,idxcurr)))','.-','Color',colors(pbgivenai,:));
  end
  title(hax(pai),sprintf('P(expression) = %f',pa));
  if pai == 1,
    s = [{sprintf('P(behavior|expression) = %s',num2str(pbgivenas(1)))}
      cellstr(num2str(pbgivenas(2:end)'))];
    legend(h,s);
  end
  set(hax(pai),'YScale','log');
  box(hax(pai),'off');
end
set(hax(ishandle(hax)),'YTick',10.^[-15,-10,-5:0]);
set(hax(:,1),'XTickLabel',{});
xlabel(hax(1,2),'P(expression|behavior)');
ylabel(hax(1,2),'p-value');
linkaxes(hax(ishandle(hax)));
set(hax(ishandle(hax)),'XLim',xlim,'YLim',[eps/10,10]);

%%

SaveFigLotsOfWays(hfig,'HypotheticalPValues20130927');

% 
% nhits = round(bsxfun(@times,pas,pbgivenas')*nlines);
% 
% 
% pv = nan(100,10,15);
% for h = 1:100,
%   for b = 1:10,
%     if b > h,
%       continue;
%     end    
%     if round(h/b) ~= h/b,
%       continue;
%     end
%     x = h/b;
%     for p = 1:15,
%       pv(h,b,p) = 1-binocdf(x-1,h,p/100);
%     end
%   end
% end
% 
% figure(1434);
% clf;
% hold on;
% colors = jet(10)*.7;
% psplot = [5];
% markersplot = {'*-','x--','s:'};
% for b = 1:10,
%   for pi = 1:numel(psplot),
%     p = psplot(pi);
%     idxcurr = find(~isnan(pv(:,b,p)));
%     plot(idxcurr,pv(idxcurr,b,p),markersplot{pi},'Color',colors(b,:));
%   end
% end

%%

maxfrac = max(fraclineswithexpression(:));
pas = (2:ceil(nlines*maxfrac))./nlines;
nbs = 2:8;

pval = nan(numel(pas),numel(nbs));
for i = 1:numel(nbs),
  nb = nbs(i);
  pval(:,i) = 1-binocdf(nb-1,nb,pas);
end

hfig = 348;
figure(hfig);
clf;
colors = jet(numel(nbs))*.7;
for i = 1:numel(nbs),
  plot(pas,pval(:,i),'.-','Color',colors(i,:));
  hold on;
end
set(gca,'YScale','log','YTick',10.^(-16:2:0),'YGrid','on','YLim',[10^-16,1])
ylim = get(gca,'YLim');
for i = 1:numel(nbs),
  j = find(pval(:,i)>10^-6,1);
  if isempty(j),
    continue;
  end
  plot(pas(j)+[0,0],ylim,':','Color','k');
end
s = cell(size(nbs));
for i = 1:numel(nbs),
  s{i} = sprintf('N. lines = %d',nbs(i));
end
legend(s);
xlabel('Fraction of lines with expression');
ylabel('p-value');
set(gca,'XLim',[0,max(pas)]);
box off;

SaveFigLotsOfWays(hfig,'HypotheticalPValues20131007');

paedges = nan(1,numel(nbs));
for i = 1:numel(nbs),
  pa = pas(find(pval(:,i)>10^-6,1));
  if isempty(pa),
    continue;
  end
  paedges(i) = pa;
  fprintf('%d: %f\n',nbs(i),pa);
end

%edges = linspace(2/nlines,maxfrac+eps,100);
idxdata = find(~isnan(paedges));
edges = [2/nlines,paedges(idxdata),maxfrac+eps];
centers = (edges(1:end-1)+edges(2:end))/2;

counts = histc(fraclineswithexpression(:),edges);
counts = counts(1:end-1);

hfig = 349;
figure(hfig);
clf;
x = [edges(1:end-1);edges(2:end)];
y = [zeros(1,numel(centers));counts'/sum(counts)];
for i = 1:numel(edges)-1,
  if i > numel(idxdata),
    color = 'k';
  else
    color = colors(idxdata(i),:);
  end
  patch(x([1,1,2,2],i),y([1,2,2,1],i),color);  
end
%plot(centers,counts/sum(counts),'.-');
box off;
xlabel('Fraction of lines with expression');
ylabel('Fraction of voxels');
set(gca,'XLim',[edges(1),edges(end)],'XTick',paedges(idxdata));
rotateticklabel(gca);

SaveFigLotsOfWays(hfig,'HypotheticalPValuesPart2_20131007');

% 2: 0.001357
% 3: 0.010403
% 4: 0.031660
% 5: 0.063320
% 6: 0.100407

%%

maxfrac = max(fraclineswithexpression(:));
pas = (2:ceil(nlines*maxfrac))./nlines;
nbs = 2:2:30;

pval = nan(numel(pas),numel(nbs));
for i = 1:numel(nbs),
  nb = nbs(i);
  pval(:,i) = 1-binocdf(nb/2-1,nb,pas);
end

hfig = 350;
figure(hfig);
clf;
colors = jet(numel(nbs))*.7;
for i = 1:numel(nbs),
  plot(pas,pval(:,i),'.-','Color',colors(i,:));
  hold on;
end
set(gca,'YScale','log','YTick',10.^(-16:2:0),'YGrid','on','Ylim',[10^-16,1])
ylim = get(gca,'YLim');
for i = 1:numel(nbs),
  j = find(pval(:,i)>10^-6,1);
  if isempty(j),
    continue;
  end
  plot(pas(j)+[0,0],ylim,':','Color','k');
end
s = cell(size(nbs));
for i = 1:numel(nbs),
  s{i} = sprintf('N. lines = %d',nbs(i));
end
legend(s,'Location','EastOutside');
xlabel('Fraction of lines with expression');
ylabel('p-value');
set(gca,'XLim',[0,max(pas)]);
box off;

SaveFigLotsOfWays(hfig,'HypotheticalPValuesNoise20131007');

paedges = nan(1,numel(nbs));
for i = 1:numel(nbs),
  pa = pas(find(pval(:,i)>10^-6,1));
  if isempty(pa),
    continue;
  end
  paedges(i) = pa;
  fprintf('%d: %f\n',nbs(i),pa);
end

%edges = linspace(2/nlines,maxfrac+eps,100);
idxdata = find(~isnan(paedges)&[paedges(1:end-1)~=paedges(2:end),true]);
edges = [2/nlines,paedges(idxdata),maxfrac+eps];
centers = (edges(1:end-1)+edges(2:end))/2;

counts = histc(fraclineswithexpression(:),edges);
counts = counts(1:end-1);

hfig = 349;
figure(hfig);
clf;
x = [edges(1:end-1);edges(2:end)];
y = [zeros(1,numel(centers));counts'/sum(counts)];
for i = 1:numel(edges)-1,
  if i > numel(idxdata),
    color = 'k';
  else
    color = colors(idxdata(i),:);
  end
  patch(x([1,1,2,2],i),y([1,2,2,1],i),color);  
end
%plot(centers,counts/sum(counts),'.-');
box off;
xlabel('Fraction of lines with expression');
ylabel('Fraction of voxels');
set(gca,'XLim',[edges(1),edges(end)],'XTick',paedges(idxdata));
rotateticklabel(gca);

SaveFigLotsOfWays(hfig,'HypotheticalPValuesNoisePart2_20131007');
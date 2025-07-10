function h = DrawBehaviorVisualization(stats,colors,statlims,x1,x2,y1,y2,varargin)

[hax,leftovers] = myparse_nocheck(varargin,'hax',[]);

if isempty(hax),
  hax = gca;
end
[nstats,nlines] = size(stats);
x1 = x1(:)';
x2 = x2(:)';
y1 = y1(:)';
y2 = y2(:)';

fracperstat = bsxfun(@rdivide,abs(stats),sum(abs(stats),1));
xs = bsxfun(@times,bsxfun(@plus,x1,cumsum([zeros(1,nlines);fracperstat])),x2-x1);
ispos = stats >= 0;
zstats = nan(size(stats));
for stati = 1:nstats,
  zstats(stati,ispos(stati,:)) = stats(stati,ispos(stati,:)) / statlims(stati,1);
  zstats(stati,~ispos(stati,:)) = stats(stati,~ispos(stati,:)) / statlims(stati,2);
end
zstats = min(1,max(-1,zstats));
ys = bsxfun(@times,bsxfun(@plus,y1,zstats,x2-x1));

h = nan(nstats+1,1);
for stati = 1:nstats,
  
  xscurr = xs(stati:stati+1,:);
  yscurr = [zeros(1,nlines);ys(stati,:)];
  
  h(stati) = patch([xscurr([1,1,2,2,1],:);nans(1,nlines)],...
    [yscurr([1,2,2,1,1],:);nan(1,nlines)],colors(stati,:),...
    'LineStyle','none',leftovers{:},'Parent',hax);
  
end

h(nstats+1) = plot([x1;x2;nans(1,nlines)],...
  repmat((y1+y2)/2,[3,1]),'--','Color',[.7,.7,.7]);
function h = semilogy_with_zeros(hax,x,y,xlim,ylim,varargin)

% get the minimum positive value for y
minvplus = min(y(y>0));
if isempty(minvplus), minvplus = .001; end
% what is the label directly below this?
logminvlabel = floor(log10(minvplus));
minvlabel = 10^logminvlabel;
% set zero to be min label / 2
zerov = minvlabel / 2;
% set axis lim to be zero / 2
limv = zerov / 2;
ylim = [limv,ylim];
% reset y = 0 to zerov
y(y==0) = zerov;

lastisinf = isinf(xlim(2));
if lastisinf,
  xlim(2) = x(end)+.5*(x(end)-x(end-1));
end
firstisinf = isinf(xlim(1));
if firstisinf,
  xlim(1) = x(1)+.5*(x(1)-x(2));
end

% plot
h = semilogy(hax,x,y,varargin{:});

% adjust axis
axis(hax,[xlim,ylim]);

% choose y-ticks
ytick = [zerov,10.^(logminvlabel:0)];

% labels corresponding to these
yticklabel = cell(1,1-logminvlabel);
for i = 1:1-logminvlabel,
  yticklabel{i} = sprintf('10^%d',i-1+logminvlabel);
end
% zero for bottom tick
yticklabel = [{'0'},yticklabel];

% apply
set(hax,'ytick',ytick,'YTickLabel',yticklabel);

% set xticklabel if inf
if lastisinf || firstisinf,
  xtick = get(hax,'xtick');
  xticklabel = cellstr(get(hax,'xticklabel'));
  if lastisinf,
    i = find(xtick == x(end),1,'last');
    if isempty(i), i = numel(xtick); end
    xtick(i) = x(end);
    xticklabel{i} = ['>=',num2str(x(end)-.5*(x(end)-x(end-1)))];
  end
  if firstisinf,
    i = find(xtick == x(1),1,'first');
    if isempty(i), i = 1; end
    xtick(i) = x(1);
    xticklabel{i} = ['<=',num2str(x(end)-.5*(x(end)-x(end-1)))];
  end
  set(hax,'xtick',xtick,'xticklabel',xticklabel);
end
function barhist(ax,dat,grp,Nbin,grpleg,varargin)
% barhist(ax,dat,grp,Nbin,grpleg)
% stacked histogram
%
% ax: axis handle.
% dat: col vec. data to be histed.
% grp: col vec. grouping vec of ints, same size as dat.
% Nbin: number of hist bins
% grpleg: cellstr, number of els equals max(grp)

assert(numel(dat)==numel(grp));

Ngrp = max(grp);
bins = linspace(min(dat),max(dat)+eps(max(dat)),Nbin);
[~,bin] = histc(dat,bins);
assert(numel(bin)==numel(dat));
Zs = hist(grp,1:Ngrp);

if ~isempty(varargin),
  if ismember('stacked',varargin),
    bartype = 'stacked';
  elseif ismember('histc',varargin),
    bartype = 'histc';
  elseif ismember('lineplot',varargin),
    bartype = 'lineplot';
  else
    bartype = 'stacked';
  end
  if ismember('normalize',varargin),
    donormalize = true;
  else
    donormalize = false;
  end
end

stackdata = zeros(Nbin,Ngrp);
for binidx = 1:Nbin
    tfThisBin = bin==binidx;
    for grpidx = 1:Ngrp
        tfThisGrp = grp==grpidx;
        stackdata(binidx,grpidx) = nnz(tfThisBin&tfThisGrp);
    end
end
if donormalize,
  stackdata = bsxfun(@rdivide,stackdata,Zs(:)');
end
  
ctrs = [(bins(1:end-1)+bins(2:end))/2 bins(end)];
if strcmpi(bartype,'lineplot'),
  plot(ax,ctrs,stackdata);
else
  bar(ax,ctrs,stackdata,1,bartype);
end
if ~isempty(grpleg)
    for i = 1:Ngrp,
      grpleg{i} = [grpleg{i},sprintf(', N=%d',Zs(i))];
    end
    legend(ax,grpleg,'Location','best');
end
set(ax,'XTick',ctrs(1:2:end));
% xticklbls = [(bins(1:end-1)+bins(2:end))/2 bins(end)];
% xticklbltxt = cell(numel(xticklbls),1);
% for c = 1:numel(xticklbls)
%     xticklbltxt{c} = sprintf('%.3g',xticklbls(c));
% end
% set(ax,'XTick',1:2:Nbin);
% set(ax,'XTickLabel',xticklbltxt(1:2:Nbin));

end
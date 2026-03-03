function fliesplot = ChooseFliesPlot(trx,ind,ions,nfliesplot)

nalive = zeros(1,trx.nflies);
firstframe = trx.firstframe;
endframe = trx.endframe;
nions = numel(ions);
for ioni = 1:nions,
  ion = ions(ioni);
  f0 = ind.starton(ion);
  nalive = nalive + double(f0 >= firstframe & f0 <= endframe);
end
nfliesplot = min(nfliesplot,nnz(nalive>0));
unique_nalive = fliplr(unique(nalive));
sortby = nan(1,trx.nflies);
for i = 1:trx.nflies,
  sortby(i) = nanmean(trx(i).velmag_ctr);
end
nselected = 0;
fliesplot = nan(1,nfliesplot);
for i = 1:numel(unique_nalive),
  nalivecurr = unique_nalive(i);
  idxcurr = find(nalive==nalivecurr);
  if numel(idxcurr) <= nfliesplot-nselected,
    fliesplot(nselected+1:nselected+numel(idxcurr)) = idxcurr;
    nselected = nselected + numel(idxcurr);
    if nselected >= nfliesplot,
      break;
    end
  else
    % all flies with nalive == nalivecurr
    sortbycurr = sortby(idxcurr);
    [~,~,idxselect] = furthestfirst(sortbycurr',nfliesplot-nselected,'start','nearmean');
    fliesplot(nselected+1:end) = sort(idxcurr(idxselect));
    break;
  end
end

function [errsummary,err] = ComputeDayError(date,linestats,setstats,allstats,metadata,varargin)

[minnexps,leftovers] = myparse_nocheck(varargin,'minnexps',2);

alldates = cellfun(@(x) x(1:8),{setstats.metadata.exp_datetime},'UniformOutput',false);
idx = find(strcmp(alldates,date) & setstats.nexps.velmag_ctr_flyany_frameany >= minnexps);

for i = 1:numel(idx),
  set_name = setstats.metadata(idx(i)).set;
  errcurr = ComputeSetError(set_name,linestats,setstats,allstats,metadata,leftovers{:});
  if i == 1,
    nfns = numel(errcurr);
    err = nan(numel(idx),nfns);
  end
  err(i,:) = errcurr;
end

meanerr = nanmean(err,2);
errsummary = struct;
errsummary.meanerr = nanmean(meanerr);
errsummary.maxerr = max(meanerr);
idxgood = find(~isnan(meanerr));
[errsummary.sortederr,order] = sort(meanerr(idxgood),'descend');
errsummary.set_names = {setstats.metadata(idx(idxgood(order))).set};
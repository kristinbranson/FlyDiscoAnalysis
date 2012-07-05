function hists = CombinePerExpHists(histsperexp)

% which exps we will look at
goodidx = ~isnan(histsperexp.Z) & histsperexp.Z > 0;

hists = struct;
hists.Z = sum(histsperexp.Z(goodidx));
hists.meanZ = nanmean(histsperexp.Z(goodidx));
hists.meanfrac_linear = nanmean(histsperexp.meanfrac_linear(goodidx,:),1);
hists.stdfrac_linear = nanstd(histsperexp.meanfrac_linear(goodidx,:),1,1);
hists.meanfrac_log = nanmean(histsperexp.meanfrac_log(goodidx,:),1);
hists.stdfrac_log = nanstd(histsperexp.meanfrac_log(goodidx,:),1,1);
hists.meanfracless_linear = nanmean(histsperexp.meanfracless_linear(goodidx,:),1);
hists.stdfracless_linear = nanstd(histsperexp.meanfracless_linear(goodidx,:),1,1);
hists.meanfracless_log = nanmean(histsperexp.meanfracless_log(goodidx,:),1);
hists.stdfracless_log = nanstd(histsperexp.meanfracless_log(goodidx,:),1,1);
hists.meanfracmore_linear = nanmean(histsperexp.meanfracmore_linear(goodidx,:),1);
hists.stdfracmore_linear = nanstd(histsperexp.meanfracmore_linear(goodidx,:),1,1);
hists.meanfracmore_log = nanmean(histsperexp.meanfracmore_log(goodidx,:),1);
hists.stdfracmore_log = nanstd(histsperexp.meanfracmore_log(goodidx,:),1,1);
hists.n = nnz(goodidx);
function hists = CombinePerExpHists(histsperexp)

% which exps we will look at
goodidx = ~isnan(histsperexp.Z);

hists = struct;
hists.Z = sum(histsperexp.Z(goodidx));
hists.meanZ = mean(histsperexp.Z(goodidx));
hists.meanfrac_linear = mean(histsperexp.meanfrac_linear(goodidx,:),1);
hists.stdfrac_linear = std(histsperexp.meanfrac_linear(goodidx,:),1,1);
hists.meanfrac_log = mean(histsperexp.meanfrac_log(goodidx,:),1);
hists.stdfrac_log = std(histsperexp.meanfrac_log(goodidx,:),1,1);
hists.meanfracless_linear = mean(histsperexp.meanfracless_linear(goodidx,:),1);
hists.stdfracless_linear = std(histsperexp.meanfracless_linear(goodidx,:),1,1);
hists.meanfracless_log = mean(histsperexp.meanfracless_log(goodidx,:),1);
hists.stdfracless_log = std(histsperexp.meanfracless_log(goodidx,:),1,1);
hists.meanfracmore_linear = mean(histsperexp.meanfracmore_linear(goodidx,:),1);
hists.stdfracmore_linear = std(histsperexp.meanfracmore_linear(goodidx,:),1,1);
hists.meanfracmore_log = mean(histsperexp.meanfracmore_log(goodidx,:),1);
hists.stdfracmore_log = std(histsperexp.meanfracmore_log(goodidx,:),1,1);
hists.n = nnz(goodidx);
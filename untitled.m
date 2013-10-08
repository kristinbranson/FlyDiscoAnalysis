all_jump_lines = importdata('/groups/branson/home/bransonk/Downloads/jumpline_hits.txt');

[meanim,minim,pvalue,maxpvaluesig,nlinesread] = ComputeAverageAnatomyAndPValue(all_jump_lines);

stati_jump = find(strcmp(statfnscurr,'fractime_flyany_framejump'));
stati_dcenter = find(strcmp(statfnscurr,'dcenter_flyany_frameany'));

clf;
plot(zdatacluster_norm(:,stati_jump),zdatacluster_norm(:,stati_dcenter),'k.');
hold on;
idx_alljump = ismember(line_names,all_jump_lines);
plot(zdatacluster_norm(idx_alljump,stati_jump),zdatacluster_norm(idx_alljump,stati_dcenter),'o','markerfacecolor',[0,.7,.7],'Color',[0,.7,.7]);
idx_alljump = ismember(line_names,linenames_jump_avoid0);
plot(zdatacluster_norm(idx_alljump,stati_jump),zdatacluster_norm(idx_alljump,stati_dcenter),'o','markerfacecolor',[.7,0,.7],'Color',[.7,0,.7]);


%%

pvaluetmp = min(1,min(pvalue_bigger(:,:),pvalue_smaller(:,:))*2);
idx = ~isnan(pvaluetmp);
[issig1,crit_p1,tmp] = fdr_bh(pvaluetmp(idx),.05,'dep','yes');
issig = nan(size(pvaluetmp));
issig(idx) = issig1;
crit_p1 = nan(size(pvaluetmp));
crit_p1(idx) = crit_p1;
nsig = sum(issig,2);
adj_pvalue = nan(size(pvaluetmp));
adj_pvalue(idx) = tmp;
% issig2 = nan(nlines,nstatscurr);
% issig2(1:end-1,idxfew) = issig;
pvalue = min(fracbigger_upperbound,fracsmaller_upperbound);
statfn = 'fractime_flyany_framewingextension';
stati = find(strcmp(statfns,statfn));

q = .01;
[issig,crit_p,adj_pvalue_we] = fdr_bh(fracbigger_upperbound(:,stati),q,'pdep','yes');
adj_pvalue_we = min(adj_pvalue_we,1);

figure(1);
clf;
plot(linestats.normmeans.(statfn)+linestats.means.(statfn)(end),adj_pvalue_we,'k.');
set(gca,'YScale','log','XScale','log');
xlabel('Fraction of time performing wing extension');
ylabel('Adjusted p-value');
box off;

figure(2);
v = linestats.normmeans.(statfn)+linestats.means.(statfn)(end);
[sortedv,order] = sort(v);
[~,reorder] = sort(order);
clf;
h = scatter(1:nlines,sortedv,[],log10(adj_pvalue_we(order)),'.');
linei = order(end);
lineis = order(end:-1:end-20);

colormap(flipud(jet(256))*.75);

xticks = [reorder(end),reorder(linei)];
set(gca,'XTick',xticks,'XTickLabel',{'pBDPGAL4U',linestats.line_names{linei}});
hcb = colorbar;
ylabel('Fraction of time performing wing extension');
xlabel('Line');
set(gca,'CLim',[log10(min(adj_pvalue_wa)),0]);
set(hcb,'YTick',[-2,log10(.05)],'YTickLabel',{'.01','.05'});
axisalmosttight;

%%

statfn = 'max_wing_angle_flyany_frameany';
stati = find(strcmp(statfns,statfn));

q = .01;
[issig,crit_p,adj_pvalue_wa] = fdr_bh(fracbigger_upperbound(:,stati),q,'pdep','yes');
adj_pvalue_wa = min(adj_pvalue_wa,1);

figure(3);
clf;
plot(linestats.normmeans.(statfn)+linestats.means.(statfn)(end),adj_pvalue_wa,'k.');
set(gca,'YScale','log','XScale','log');
xlabel('Maximum wing angle');
ylabel('Adjusted p-value');
box off;

figure(4);
v = linestats.normmeans.(statfn)+linestats.means.(statfn)(end);
[sortedv,order] = sort(v);
[~,reorder] = sort(order);
clf;
h = scatter(1:nlines,sortedv,[],log10(adj_pvalue_wa(order)),'.');
set(h,'SizeData',24);

xticks = [reorder(end),reorder(linei)];
xticklabels = {'pBDPGAL4U',linestats.line_names{linei}};
[xticks,tmp] = sort(xticks);
xticklabels = xticklabels(tmp);
set(gca,'XTick',xticks,'XTickLabel',xticklabels);
colormap(flipud(jet(256))*.75);
hcb = colorbar;
ylabel('Maximum wing angle');
xlabel('Line');
set(hcb,'YTick',[-2,log10(.05)],'YTickLabel',{'.01','.05'});
set(gca,'CLim',[log10(min(adj_pvalue_wa)),0]);
set(hcb,'YTick',[-2,log10(.05)],'YTickLabel',{'.01','.05'});
axisalmosttight;

%%

pvalue = min(fracbigger_upperbound,fracsmaller_upperbound);
statfn = 'fractime_flyany_framebackup';
stati = find(strcmp(statfns,statfn));

q = .01;
[issig,crit_p,adj_pvalue_b] = fdr_bh(fracbigger_upperbound(:,stati),q,'pdep','yes');
adj_pvalue_b = min(adj_pvalue_b,1);

figure(1);
clf;
plot(linestats.normmeans.(statfn)+linestats.means.(statfn)(end),adj_pvalue_b,'k.');
set(gca,'YScale','log','XScale','log');
xlabel('Fraction of time performing wing extension');
ylabel('Adjusted p-value');
box off;

figure(2);
v = linestats.normmeans.(statfn)+linestats.means.(statfn)(end);
[sortedv,order] = sort(v);
[~,reorder] = sort(order);
clf;
h = scatter(1:nlines,sortedv,[],log10(adj_pvalue_b(order)),'o','filled');
set(h,'SizeData',24);
linei = order(end);
lineis = order(end:-1:end-20);

colormap(flipud(jet(256))*.75);

xticks = [reorder(end),reorder(linei)];
set(gca,'XTick',xticks,'XTickLabel',{'pBDPGAL4U',linestats.line_names{linei}});
hcb = colorbar;
ylabel('Fraction of time backing up');
xlabel('Line');
set(gca,'CLim',[log10(min(adj_pvalue_du)),0]);
set(hcb,'YTick',[-2,log10(.05)],'YTickLabel',{'.01','.05'});
axisalmosttight;

%%

statfn = 'du_ctr_flyany_frameany';
stati = find(strcmp(statfns,statfn));

q = .01;
[issig,crit_p,adj_pvalue_du] = fdr_bh(fracsmaller_upperbound(:,stati),q,'pdep','yes');
adj_pvalue_du = min(adj_pvalue_du,1);

figure(3);
clf;
plot(linestats.normmeans.(statfn)+linestats.means.(statfn)(end),adj_pvalue_du,'k.');
set(gca,'YScale','log','XScale','log');
xlabel('Maximum wing angle');
ylabel('Adjusted p-value');
box off;

figure(4);
v = linestats.normmeans.(statfn)+linestats.means.(statfn)(end);
[sortedv,order] = sort(v,'descend');
[~,reorder] = sort(order);
clf;
h = scatter(1:nlines,sortedv,[],log10(adj_pvalue_du(order)),'o','filled');
set(h,'SizeData',24);

xticks = [reorder(end),reorder(linei)];
xticklabels = {'pBDPGAL4U',linestats.line_names{linei}};
[xticks,tmp] = sort(xticks);
xticklabels = xticklabels(tmp);
set(gca,'XTick',xticks,'XTickLabel',xticklabels);
colormap(flipud(jet(256))*.75);
hcb = colorbar;
ylabel('Forward velocity');
xlabel('Line');
set(hcb,'YTick',[-2,log10(.05)],'YTickLabel',{'.01','.05'});
set(gca,'CLim',[log10(min(adj_pvalue_du)),0]);
set(hcb,'YTick',[-2,log10(.05)],'YTickLabel',{'.01','.05'});
axisalmosttight;

%%

pvalue = min(fracbigger_upperbound,fracsmaller_upperbound);
statfn = 'fractime_flyany_framejump';
stati = find(strcmp(statfns,statfn));

q = .01;
[issig,crit_p,adj_pvalue_j] = fdr_bh(fracbigger_upperbound(:,stati),q,'pdep','yes');
adj_pvalue_j = min(adj_pvalue_j,1);

figure(1);
clf;
plot(linestats.normmeans.(statfn)+linestats.means.(statfn)(end),adj_pvalue_j,'k.');
set(gca,'YScale','log','XScale','log');
xlabel('Fraction of time jumping');
ylabel('Adjusted p-value');
box off;

figure(2);
v = linestats.normmeans.(statfn)+linestats.means.(statfn)(end);
[sortedv,order] = sort(v);
[~,reorder] = sort(order);
clf;
h = scatter(1:nlines,sortedv,[],log10(adj_pvalue_j(order)),'o','filled');
set(h,'SizeData',24);
linei = order(end);
lineis = order(end:-1:end-20);

colormap(flipud(jet(256))*.75);

xticks = [reorder(end)];
set(gca,'XTick',xticks,'XTickLabel',{'pBDPGAL4U'});
hcb = colorbar;
ylabel('Fraction of time jumping');
xlabel('Line');
set(hcb,'YTick',[-2,log10(.05)],'YTickLabel',{'.01','.05'});
axisalmosttight;
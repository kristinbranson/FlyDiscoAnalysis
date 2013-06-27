%% load in data

load CollectedPerFrameStatsMB20130317.mat;
statfns = fieldnames(allstats);
nstats = numel(statfns);

metadata = metadata(idxanalyze);
for i = 1:numel(statfns),
  allstats.(statfns{i})(~idxanalyze) = [];
end
nlines = numel(linestats.line_names);
exp2set = setidx;

linecontrolidx = find(strcmp(linestats.line_names,'pBDPGAL4U'));
expiscontrol = strcmp({metadata.line_name},'pBDPGAL4U');
expidxcontrol = find(expiscontrol);
setiscontrol = set2lineidx == linecontrolidx;
setidxcontrol = find(setiscontrol);

%% parameters

nsamples = 100000;
epsilon_compare = .000001;

statfns_plot = {
  'fractime_flyany_framestop'
  'fractime_flyany_framewalk'
  'fractime_flyany_framejump'
  'fractime_flyany_framerighting'
  'fractime_flyany_framepivottail'
  'fractime_flyany_framepivotcenter'
  'fractime_flyany_framebodyturn'
  'fractime_flyany_framebackup'
  'fractime_flyany_framecrabwalkall'
  'fractime_flyany_framecrabwalkextreme'
  'fractime_flyany_framechase'
  'fractime_flyany_frameattemptedcopulation'};

[ism,idx] = ismember(statfns,statfns_plot);
statis = find(ism);
[~,order] = sort(idx(ism));
statis = statis(order);

%% load in sigmahat_set and sigmahat_exp weights

load VarianceEstimates20130314.mat;

%% compute means for each line

[~,exp2lineidx] = ismember({metadata.line_name},linestats.line_names);
muhat_weighted_perstat = nan(nlines,nstats);
sigmahat_weighted_perstat = nan(nlines,nstats);

for stati = 1:nstats,
  
  statfn = statfns{stati};
  
  x = allstats.(statfn);
  badexpidx = isnan(x);
  nexpsperset = setnexps(:,stati);
  n = nexpsperset(exp2set);
  
  sigmahat_set = sigmahat_set_perstat(stati);
  sigmahat_exp = sigmahat_exp_perstat(stati);
  
  expweights = 1./(n.*sigmahat_set^2 + sigmahat_exp^2);
  expweights = expweights / sum(expweights);
  
  sum_perline = accumarray(exp2lineidx(~badexpidx)',x(~badexpidx)'.*expweights(~badexpidx));
  sum2_perline = accumarray(exp2lineidx(~badexpidx)',x(~badexpidx).^2'.*expweights(~badexpidx));
  z_perline = accumarray(exp2lineidx(~badexpidx)',expweights(~badexpidx));
  muhat_perline = sum_perline ./ z_perline;
  sigmahat_perline = sqrt(sum2_perline ./ z_perline - muhat_perline.^2);
  
  muhat_weighted_perstat(:,stati) = muhat_perline;
  sigmahat_weighted_perstat(:,stati) = sigmahat_perline;
  
end

% %% if we just weight each experiment equally
% 
% muhat_set_perstat = nan(nlines,nstats);
% sigmahat_set_perstat = nan(nlines,nstats);
% 
% for stati = 1:nstats,
%   
%   statfn = statfns{stati};
%    
%   x = setstats.means.(statfn);
%   badsetidx = isnan(x);
%   
%   sum_perline = accumarray(set2lineidx(~badsetidx)',x(~badsetidx)');
%   sum2_perline = accumarray(set2lineidx(~badsetidx)',x(~badsetidx).^2');
%   z_perline = accumarray(set2lineidx(~badsetidx)',ones(nnz(~badsetidx),1));
%   muhat_perline = sum_perline ./ z_perline;
%   sigmahat_perline = sqrt(sum2_perline ./ z_perline - muhat_perline.^2);
%   
%   muhat_set_perstat(:,stati) = muhat_perline;
%   sigmahat_set_perstat(:,stati) = sigmahat_perline;
%   
% end
% 
% %% plot 
% 
% hfig = 123;
% figure(hfig);
% clf;
% 
% nc = ceil(sqrt(numel(statis)));
% nr = ceil(numel(statis)/nc);
% 
% hax = createsubplots(nr,nc,.05);
% 
% for statii = 1:numel(statis),
%   
%   stati = statis(statii);
% 
%   statfn = statfns{stati};
%   okidx = find(~isnan(linestats.nsets.(statfn)));
%   [ns,~,nidx] = unique(linestats.nsets.(statfn)(okidx));
%   colors = jet(numel(ns))*.7;
%   
%   axes(hax(statii));
%   maxv = max(max(muhat_weighted_perstat(:,stati)),max(muhat_set_perstat(:,stati)))*1.05;
%   plot([0,maxv],[0,maxv],'k-');
%   hold on;
%   for i = 1:numel(ns),    
%     plot(muhat_weighted_perstat(okidx(nidx==i),stati),muhat_set_perstat(okidx(nidx==i),stati),'.','Color',colors(i,:));    
%   end
%   axis equal tight;
%   title(statfn,'Interpreter','none');
%   xlabel('weighted mean');
%   ylabel('unweighted mean');
%   
% end
% 
% hfig = 124;
% figure(hfig);
% clf;
% hax = createsubplots(nr,nc,.05);
% 
% 
% for statii = 1:numel(statis),
%   
%   stati = statis(statii);
% 
%   statfn = statfns{stati};
%   okidx = setdiff(find(~isnan(linestats.nsets.(statfn))),linecontrolidx);
%   [ns,~,nidx] = unique(linestats.nsets.(statfn)(okidx));
%   
%   err = abs(muhat_weighted_perstat(okidx,stati)-muhat_set_perstat(okidx,stati));
%   
%   axes(hax(statii));
%   plot(linestats.nsets.(statfn)(okidx),err,'.');
%   axisalmosttight;
%   title(statfn,'Interpreter','none');
%   xlabel('N. sets');
%   ylabel('Difference in mean estimates');
%   
% end
% 
% %% what does the distribution of set variances look like?
% 
% hfig = 125;
% figure(hfig);
% clf;
% 
% hax = createsubplots(nr,nc,.05);
% 
% 
% for statii = 1:numel(statis),
%   
%   stati = statis(statii);
%   statfn = statfns{stati};
%   x = setstats.means.(statfn);
%   err = nan(size(x));
%   
%   for seti = 1:numel(x),
%     
%     linei = set2lineidx(seti);
%     err(seti) = abs(x(seti) - muhat_set_perstat(linei,stati));
%     
%   end
%   
%   [counts,centers] = hist(err);
%   
%   axes(hax(statii));
%   scatter(muhat_set_perstat(set2lineidx,stati),err,[],set2lineidx,'.');
%   %axisalmosttight;
%   title(statfn,'Interpreter','none');
%   xlabel('Line mean');
%   ylabel('Difference in set and line means');
%   
% end
% 
% %% distribution of difference in set and line means
% 
% 
% hfig = 126;
% figure(hfig);
% clf;
% 
% hax = createsubplots(nr,nc,.05);
% 
% 
% for statii = 1:numel(statis),
%   
%   stati = statis(statii);
%   statfn = statfns{stati};
%   x = setstats.means.(statfn);
%   err = nan(size(x));
%   
%   for seti = 1:numel(x),
%     
%     linei = set2lineidx(seti);
%     err(seti) = abs(x(seti) - muhat_set_perstat(linei,stati));
%     
%   end
%   
%   [counts,centers] = hist(err);
%   frac = counts / sum(counts);
%   axes(hax(statii));
%   bar(centers,counts);
%   %axisalmosttight;
%   title(statfn,'Interpreter','none');
%   xlabel('Difference in set and line means');
%   ylabel('Fraction of sets');
%   
% end

%% randomly construct similar-sized data sets

nexpsperset_control = setstats.nexps.(statfn)(setiscontrol);
[~,exp2set_control] = ismember(setidx(expiscontrol),setidxcontrol);
set2exp_control = cell(1,numel(setidxcontrol));
for i = 1:numel(setidxcontrol),
  set2exp_control{i} = find(exp2set_control == i);
end
fracsmallersamples = nan(nlines,nstats);
fracbiggersamples = nan(nlines,nstats);

idx = exp2lineidx==linecontrolidx;
xexp_control = nan(numel(statis),nnz(idx));

for statii = 1:numel(statis),
  
  stati = statis(statii);
  statfn = statfns{stati};
  xexp_control(statii,:) = allstats.(statfn)(idx);
  
end

sigmahat_set = sigmahat_set_perstat(statis);
sigmahat_exp = sigmahat_exp_perstat(statis);

for linei = 1:nlines,
%for linei = linecontrolidx,
  
  muline = muhat_weighted_perstat(linei,statis);
  exp2lineidx_curr = exp2lineidx == linei;
  set2lineidx_curr = set2lineidx == linei;
  nexpsperset_line = nexpsperset(set2lineidx_curr);
  [fracsmallersamples(linei,statis),fracbiggersamples(linei,statis)] = ...
    ComputePValueRandomPermutations(muline,nexpsperset_line,...
    sigmahat_set,sigmahat_exp,nsamples,...
    xexp_control,nexpsperset_control,set2exp_control);

  for stati = statis(:)',
    statfn = statfns{stati};
    fprintf('%s, %s:\n',statfn,linestats.line_names{linei});
    fprintf('Percentile of %s mean among sample means: %.1f%%\n',...
      linestats.line_names{linei},(fracsmallersamples(linei,stati)+(1-fracbiggersamples(linei,stati)))/2*100);
  end
end
  
%% choose the most significant test, and correct for this

pvalue = min(1,min(1-fracsmallersamples,1-fracbiggersamples)*2);
min_pvalue_observable = 2 / nsamples;

%% make a table

% sort by mean for statfns_plot{1};
stati = find(strcmp(statfns,statfns_plot{1}));
[~,lineorder] = sort(muhat_weighted_perstat(:,stati),'descend');

% normalize strength
ismaybesig = pvalue < .05;
dcontrol = bsxfun(@minus,muhat_weighted_perstat,muhat_weighted_perstat(linecontrolidx,:));
zscore = bsxfun(@rdivide,bsxfun(@minus,muhat_weighted_perstat,muhat_weighted_perstat(linecontrolidx,:)),sigmahat_weighted_perstat(linecontrolidx,:));
colorby = dcontrol;
% minv = prctile(colorby(ismaybesig),1);
% maxv = prctile(colorby(ismaybesig),99);
tmp = colorby;
tmp(~ismaybesig) = nan;
minv = min(tmp,[],1);
maxv = max(tmp,[],1);
minv = min(minv,0);
maxv = max(maxv,0);

log10pvalue = log10(max(min_pvalue_observable,pvalue));

ncolors = 256;
cm = redblue(2*ncolors);
% cm = zeros(2*ncolors,3);
% cm(:,3) = linspace(1,0,2*ncolors);
% cm(:,1) = linspace(0,1,2*ncolors);
% cm(1:ncolors,2) = linspace(0,1,ncolors);
% cm(ncolors+1:end,2) = linspace(1,0,ncolors);
cm = bsxfun(@rdivide,cm,sqrt(sum(cm.^2,2)));

isneg = colorby < 0;
coloridx = nan(size(colorby));
z = -minv;
z(minv==0) = 1;
tmp = max(1,floor(bsxfun(@rdivide,bsxfun(@minus,colorby,minv),z)*ncolors)+1);
coloridx(isneg) = tmp(isneg);
z = maxv;
z(maxv==0) = 1;
tmp = min(ncolors,floor(bsxfun(@rdivide,colorby,z)*ncolors)+1)+ncolors;
coloridx(~isneg) = tmp(~isneg);
coloridx = min(2*ncolors,max(1,coloridx));

hue = reshape(cm(coloridx,:),[nlines,nstats,3]);

maxp = max(log10pvalue(:));
minp = min(log10pvalue(:));
intensity = 1-(log10pvalue-minp)/(maxp-minp);
color = 1-bsxfun(@times,1-hue,intensity);

color_ordered = color(lineorder,:,:);

[ism,idx] = ismember(statfns,statfns_plot);
statis = find(ism);
[~,order] = sort(idx(ism));
statis = statis(order);

color_ordered = color_ordered(:,statis,:);

hfig = 1;
figure(hfig);
clf;
set(hfig,'Units','pixels','Position',[10,10,1200,1500]);

image(color_ordered);

behaviornames = regexprep(statfns_plot,'^fractime_flyany_frame','');
set(gca,'XTick',1:numel(statis),'XTickLabel',behaviornames);

short_line_names = regexprep(linestats.line_names(lineorder),'^GMR_','');
short_line_names = regexprep(short_line_names,'_AE_01$','');
short_line_names = regexprep(short_line_names,'_01$','');
short_line_names = regexprep(short_line_names,'_BX3G_1500154$','');

set(gca,'YTick',1:nlines,'YTickLabel',short_line_names,'box','off');

htext = nan(nlines,numel(statis));
for lineii = 1:nlines,
  linei = lineorder(lineii);
  for statii = 1:numel(statis),
    stati = statis(statii);
    if pvalue(linei,stati) >= .05,
      stars = '   ';
    elseif pvalue(linei,stati) >= .01,
      stars = '*  ';
    elseif pvalue(linei,stati) == 0,
      stars = '***';
    else
      nstars = min(3,floor(-log10(pvalue(linei,stati))));
      stars = [repmat('*',[1,nstars]),repmat(' ',[1,3-nstars])];
    end
    htext(lineii,statii) = text(statii,lineii,...
      sprintf('%.4f%s',muhat_weighted_perstat(linei,stati),stars),...
      'HorizontalAlignment','center','VerticalAlignment','middle');
    if linei == linecontrolidx,
      set(htext(lineii,statii),'FontWeight','bold');
    end
    if color_ordered(lineii,statii,3)/color_ordered(lineii,statii,1) > 1.5,
      set(htext(lineii,statii),'Color',[.99,.99,.99]);
    end
  end
end
title('Fraction of time performing each behavior');
set(gca,'Position',[.04,.02,.905,.95]);
set(gca,'TickLength',[0,0]);

% plot lines every 5 lines so that we can line things up
hline = [];
hold on;
for i = 5:5:nlines,
  hline(end+1) = plot([.5,numel(statis)+.5],[i,i]+.5,'k-');
end
i = find(lineorder==linecontrolidx);
hline(end+1) = plot([.5,numel(statis)+.5],[i,i]+.5,'k-');
hline(end+1) = plot([.5,numel(statis)+.5],[i,i]-.5,'k-');

hcolorbar_ax = axes('Position',[.95,.625,.045,.25]);
image(permute(flipud(cm),[1,3,2]));
axis off;
text(1,size(cm,1),'Less','Color',[.99,.99,.99],...
  'HorizontalAlignment','center','VerticalAlignment','bottom','FontName','times');
text(1,1,'More','Color','k',...
  'HorizontalAlignment','center','VerticalAlignment','top','FontName','times');
text(1,ncolors+.5,'Same','Color','k',...
  'HorizontalAlignment','center','VerticalAlignment','middle','FontName','times');

hcolorbar_ax2 = axes('Position',[.95,.125,.045,.25]);
image(permute(gray(ncolors*2),[1,3,2]));
axis off;
text(1,2*ncolors,'p\geq.05','Color','k',...
  'HorizontalAlignment','center','VerticalAlignment','bottom','FontName','times');
text(1,1,{'***','p<.001'},'Color',[.99,.99,.99],...
  'HorizontalAlignment','center','VerticalAlignment','top','FontName','times');
text(1,2*ncolors*(log10(.01)-minp)/(maxp-minp),{'**','p<.01'},'Color',[.99,.99,.99],...
  'HorizontalAlignment','center','VerticalAlignment','middle','FontName','times');
text(1,2*ncolors*(log10(.05)-minp)/(maxp-minp),{'*','p<.05'},'Color',[.99,.99,.99],...
  'HorizontalAlignment','center','VerticalAlignment','middle','FontName','times');

SaveFigLotsOfWays(hfig,'MBConfidenceAndStrength20130321');

%% make a table showing lines in the set order

% sort by mean for statfns_plot{1};
stati = find(strcmp(statfns,statfns_plot{1}));
line_name_order = importdata('MBLineNameOrder.txt');
short_line_names = regexprep(linestats.line_names,'^GMR_','');
short_line_names = regexprep(short_line_names,'_AE_01$','');
short_line_names = regexprep(short_line_names,'_01$','');
short_line_names = regexprep(short_line_names,'_BX3G_1500154$','');

[lineisrun,lineorder] = ismember(line_name_order,short_line_names);

% normalize strength
ismaybesig = pvalue < .05;
dcontrol = bsxfun(@minus,muhat_weighted_perstat,muhat_weighted_perstat(linecontrolidx,:));
tmp = dcontrol;
tmp(~ismaybesig) = nan;
minv = min(tmp,[],1);
minv = min(minv,0);
maxv = max(tmp,[],1);
maxv = max(maxv,0);

log10pvalue = log10(max(min_pvalue_observable,pvalue));

ncolors = 256;
cm = redblue(2*ncolors);
% cm = zeros(2*ncolors,3);
% cm(:,3) = linspace(1,0,2*ncolors);
% cm(:,1) = linspace(0,1,2*ncolors);
% cm(1:ncolors,2) = linspace(0,1,ncolors);
% cm(ncolors+1:end,2) = linspace(1,0,ncolors);
cm = bsxfun(@rdivide,cm,sqrt(sum(cm.^2,2)));

isneg = dcontrol < 0;
coloridx = nan(size(dcontrol));
z = -minv;
z(minv==0) = 1;
tmp = max(1,floor(bsxfun(@rdivide,bsxfun(@minus,dcontrol,minv),z)*ncolors)+1);
coloridx(isneg) = tmp(isneg);
z = maxv;
z(maxv==0) = 1;
tmp = min(ncolors,floor(bsxfun(@rdivide,dcontrol,z)*ncolors)+1)+ncolors;
coloridx(~isneg) = tmp(~isneg);
coloridx = min(2*ncolors,max(1,coloridx));

hue = reshape(cm(coloridx,:),[nlines,nstats,3]);

maxp = max(log10pvalue(:));
minp = min(log10pvalue(:));
intensity = 1-(log10pvalue-minp)/(maxp-minp);
color = 1-bsxfun(@times,1-hue,intensity);

color_ordered = ones(numel(lineorder),nstats,3);
color_ordered(lineisrun,:,:) = color(lineorder(lineisrun),:,:);

[ism,idx] = ismember(statfns,statfns_plot);
statis = find(ism);
[~,order] = sort(idx(ism));
statis = statis(order);

color_ordered = color_ordered(:,statis,:);

hfig = 2;
figure(hfig);
clf;
set(hfig,'Units','pixels','Position',[10,10,1200,1500]);

image(color_ordered);

behaviornames = regexprep(statfns_plot,'^fractime_flyany_frame','');
set(gca,'XTick',1:numel(statis),'XTickLabel',behaviornames);
set(gca,'YTick',1:nlines,'YTickLabel',line_name_order,'box','off');

htext = nan(nlines,numel(statis));
for lineii = 1:numel(lineorder),
  linei = lineorder(lineii);
  if linei == 0,
    continue;
  end
  for statii = 1:numel(statis),
    stati = statis(statii);
    if pvalue(linei,stati) >= .05,
      stars = '   ';
    elseif pvalue(linei,stati) >= .01,
      stars = '*  ';
    elseif pvalue(linei,stati) == 0,
      stars = '***';
    else
      nstars = min(3,floor(-log10(pvalue(linei,stati))));
      stars = [repmat('*',[1,nstars]),repmat(' ',[1,3-nstars])];
    end
    htext(lineii,statii) = text(statii,lineii,...
      sprintf('%.4f%s',muhat_weighted_perstat(linei,stati),stars),...
      'HorizontalAlignment','center','VerticalAlignment','middle');
    if linei == linecontrolidx,
      set(htext(lineii,statii),'FontWeight','bold');
    end
    if color_ordered(lineii,statii,3)/color_ordered(lineii,statii,1) > 1.5,
      set(htext(lineii,statii),'Color',[.99,.99,.99]);
    end
  end
end
title('Fraction of time performing each behavior');
set(gca,'Position',[.04,.02,.905,.95]);
set(gca,'TickLength',[0,0]);

% plot lines every 5 lines so that we can line things up
hline = [];
hold on;
for i = 5:5:nlines,
  hline(end+1) = plot([.5,numel(statis)+.5],[i,i]+.5,'k-');
end

hcolorbar_ax = axes('Position',[.95,.625,.045,.25]);
image(permute(flipud(cm),[1,3,2]));
axis off;
text(1,size(cm,1),'Less','Color',[.99,.99,.99],...
  'HorizontalAlignment','center','VerticalAlignment','bottom','FontName','times');
text(1,1,'More','Color','k',...
  'HorizontalAlignment','center','VerticalAlignment','top','FontName','times');
text(1,ncolors+.5,'Same','Color','k',...
  'HorizontalAlignment','center','VerticalAlignment','middle','FontName','times');

hcolorbar_ax2 = axes('Position',[.95,.125,.045,.25]);
image(permute(gray(ncolors*2),[1,3,2]));
axis off;
text(1,2*ncolors,'p\geq.05','Color','k',...
  'HorizontalAlignment','center','VerticalAlignment','bottom','FontName','times');
text(1,1,{'***','p<.001'},'Color',[.99,.99,.99],...
  'HorizontalAlignment','center','VerticalAlignment','top','FontName','times');
text(1,2*ncolors*(log10(.01)-minp)/(maxp-minp),{'**','p<.01'},'Color',[.99,.99,.99],...
  'HorizontalAlignment','center','VerticalAlignment','middle','FontName','times');
text(1,2*ncolors*(log10(.05)-minp)/(maxp-minp),{'*','p<.05'},'Color',[.99,.99,.99],...
  'HorizontalAlignment','center','VerticalAlignment','middle','FontName','times');

SaveFigLotsOfWays(hfig,'MBFlyBowlResultsFigure20130321');

%% output p-values to csv file

fid = fopen('MBFlyBowlResultsPValues20130317.csv','w');
fprintf(fid,'line');
fprintf(fid,',%s',statfns{statis});
fprintf(fid,'\n');
for i = 1:numel(lineorder),
  
  ii = lineorder(i);
  fprintf(fid,'%s',line_name_order{i});  
  if i ~= 0,
    fprintf(fid,',%f',sign(dcontrol(i,statis)).*pvalue(i,statis));
  end
  fprintf(fid,'\n');
end
fclose(fid);

%% benjamini correction

m = numel(statis)*nlines;
pvalue1 = max(pvalue,min_pvalue_observable);
q = .01;
[issig,crit_p,tmp] = fdr_bh(pvalue1(:,statis),q,'dep','yes');
adj_pvalue = nan(size(pvalue));
adj_pvalue(:,statis) = tmp;

%% make a table

% sort by mean for statfns_plot{1};
stati = find(strcmp(statfns,statfns_plot{1}));
[~,lineorder] = sort(muhat_weighted_perstat(:,stati),'descend');

% normalize strength
ismaybesig = adj_pvalue < .05;
dcontrol = bsxfun(@minus,muhat_weighted_perstat,muhat_weighted_perstat(linecontrolidx,:));
zscore = bsxfun(@rdivide,bsxfun(@minus,muhat_weighted_perstat,muhat_weighted_perstat(linecontrolidx,:)),sigmahat_weighted_perstat(linecontrolidx,:));
colorby = dcontrol;
% minv = prctile(colorby(ismaybesig),1);
% maxv = prctile(colorby(ismaybesig),99);
tmp = colorby;
tmp(~ismaybesig) = nan;
minv = min(tmp,[],1);
maxv = max(tmp,[],1);
minv = min(minv,0);
maxv = max(maxv,0);

log10pvalue = log10(adj_pvalue);

ncolors = 256;
cm = redblue(2*ncolors);
% cm = zeros(2*ncolors,3);
% cm(:,3) = linspace(1,0,2*ncolors);
% cm(:,1) = linspace(0,1,2*ncolors);
% cm(1:ncolors,2) = linspace(0,1,ncolors);
% cm(ncolors+1:end,2) = linspace(1,0,ncolors);
cm = bsxfun(@rdivide,cm,sqrt(sum(cm.^2,2)));

isneg = colorby < 0;
coloridx = nan(size(colorby));
z = -minv;
z(minv==0) = 1;
tmp = max(1,floor(bsxfun(@rdivide,bsxfun(@minus,colorby,minv),z)*ncolors)+1);
coloridx(isneg) = tmp(isneg);
z = maxv;
z(maxv==0) = 1;
tmp = min(ncolors,floor(bsxfun(@rdivide,colorby,z)*ncolors)+1)+ncolors;
coloridx(~isneg) = tmp(~isneg);
coloridx = min(2*ncolors,max(1,coloridx));

hue = reshape(cm(coloridx,:),[nlines,nstats,3]);

maxp = max(log10pvalue(:));
minp = min(log10pvalue(:));
intensity = 1-(log10pvalue-minp)/(maxp-minp);
color = 1-bsxfun(@times,1-hue,intensity);

color_ordered = color(lineorder,:,:);

[ism,idx] = ismember(statfns,statfns_plot);
statis = find(ism);
[~,order] = sort(idx(ism));
statis = statis(order);

color_ordered = color_ordered(:,statis,:);

hfig = 3;
figure(hfig);
clf;
set(hfig,'Units','pixels','Position',[10,10,1200,1500]);

image(color_ordered);

behaviornames = regexprep(statfns_plot,'^fractime_flyany_frame','');
set(gca,'XTick',1:numel(statis),'XTickLabel',behaviornames);

short_line_names = regexprep(linestats.line_names(lineorder),'^GMR_','');
short_line_names = regexprep(short_line_names,'_AE_01$','');
short_line_names = regexprep(short_line_names,'_01$','');
short_line_names = regexprep(short_line_names,'_BX3G_1500154$','');

set(gca,'YTick',1:nlines,'YTickLabel',short_line_names,'box','off');

htext = nan(nlines,numel(statis));
for lineii = 1:nlines,
  linei = lineorder(lineii);
  for statii = 1:numel(statis),
    stati = statis(statii);
    if adj_pvalue(linei,stati) >= .05,
      stars = '   ';
    elseif adj_pvalue(linei,stati) >= .01,
      stars = '*  ';
    elseif adj_pvalue(linei,stati) == 0,
      stars = '***';
    else
      nstars = min(3,floor(-log10(adj_pvalue(linei,stati))));
      stars = [repmat('*',[1,nstars]),repmat(' ',[1,3-nstars])];
    end
    htext(lineii,statii) = text(statii,lineii,...
      sprintf('%.4f%s',muhat_weighted_perstat(linei,stati),stars),...
      'HorizontalAlignment','center','VerticalAlignment','middle');
    if linei == linecontrolidx,
      set(htext(lineii,statii),'FontWeight','bold');
    end
    if color_ordered(lineii,statii,3)/color_ordered(lineii,statii,1) > 1.5,
      set(htext(lineii,statii),'Color',[.99,.99,.99]);
    end
  end
end
title('Fraction of time performing each behavior');
set(gca,'Position',[.04,.02,.905,.95]);
set(gca,'TickLength',[0,0]);

% plot lines every 5 lines so that we can line things up
hline = [];
hold on;
for i = 5:5:nlines,
  hline(end+1) = plot([.5,numel(statis)+.5],[i,i]+.5,'k-');
end
i = find(lineorder==linecontrolidx);
hline(end+1) = plot([.5,numel(statis)+.5],[i,i]+.5,'k-');
hline(end+1) = plot([.5,numel(statis)+.5],[i,i]-.5,'k-');

hcolorbar_ax = axes('Position',[.95,.625,.045,.25]);
image(permute(flipud(cm),[1,3,2]));
axis off;
text(1,size(cm,1),'Less','Color',[.99,.99,.99],...
  'HorizontalAlignment','center','VerticalAlignment','bottom','FontName','times');
text(1,1,'More','Color','k',...
  'HorizontalAlignment','center','VerticalAlignment','top','FontName','times');
text(1,ncolors+.5,'Same','Color','k',...
  'HorizontalAlignment','center','VerticalAlignment','middle','FontName','times');

hcolorbar_ax2 = axes('Position',[.95,.125,.045,.25]);
image(permute(gray(ncolors*2),[1,3,2]));
axis off;
text(1,2*ncolors,'p\geq.05','Color','k',...
  'HorizontalAlignment','center','VerticalAlignment','bottom','FontName','times');
text(1,1,{'***','p<.001'},'Color',[.99,.99,.99],...
  'HorizontalAlignment','center','VerticalAlignment','top','FontName','times');
text(1,2*ncolors*(log10(.01)-minp)/(maxp-minp),{'**','p<.01'},'Color',[.99,.99,.99],...
  'HorizontalAlignment','center','VerticalAlignment','middle','FontName','times');
text(1,2*ncolors*(log10(.05)-minp)/(maxp-minp),{'*','p<.05'},'Color',[.99,.99,.99],...
  'HorizontalAlignment','center','VerticalAlignment','middle','FontName','times');

SaveFigLotsOfWays(hfig,'MBConfidenceAndStrength20130321Corrected');

%% make a table showing lines in the set order with corrected p-values

% sort by mean for statfns_plot{1};
stati = find(strcmp(statfns,statfns_plot{1}));
line_name_order = importdata('MBLineNameOrder.txt');
short_line_names = regexprep(linestats.line_names,'^GMR_','');
short_line_names = regexprep(short_line_names,'_AE_01$','');
short_line_names = regexprep(short_line_names,'_01$','');
short_line_names = regexprep(short_line_names,'_BX3G_1500154$','');

[lineisrun,lineorder] = ismember(line_name_order,short_line_names);

% normalize strength
ismaybesig = adj_pvalue < .05;
dcontrol = bsxfun(@minus,muhat_weighted_perstat,muhat_weighted_perstat(linecontrolidx,:));
tmp = dcontrol;
tmp(~ismaybesig) = nan;
minv = min(tmp,[],1);
minv = min(minv,0);
maxv = max(tmp,[],1);
maxv = max(maxv,0);

log10pvalue = log10(adj_pvalue);

ncolors = 256;
cm = redblue(2*ncolors);
% cm = zeros(2*ncolors,3);
% cm(:,3) = linspace(1,0,2*ncolors);
% cm(:,1) = linspace(0,1,2*ncolors);
% cm(1:ncolors,2) = linspace(0,1,ncolors);
% cm(ncolors+1:end,2) = linspace(1,0,ncolors);
cm = bsxfun(@rdivide,cm,sqrt(sum(cm.^2,2)));

isneg = dcontrol < 0;
coloridx = nan(size(dcontrol));
z = -minv;
z(minv==0) = 1;
tmp = max(1,floor(bsxfun(@rdivide,bsxfun(@minus,dcontrol,minv),z)*ncolors)+1);
coloridx(isneg) = tmp(isneg);
z = maxv;
z(maxv==0) = 1;
tmp = min(ncolors,floor(bsxfun(@rdivide,dcontrol,z)*ncolors)+1)+ncolors;
coloridx(~isneg) = tmp(~isneg);
coloridx = min(2*ncolors,max(1,coloridx));

hue = reshape(cm(coloridx,:),[nlines,nstats,3]);

maxp = max(log10pvalue(:));
minp = min(log10pvalue(:));
intensity = 1-(log10pvalue-minp)/(maxp-minp);
color = 1-bsxfun(@times,1-hue,intensity);

color_ordered = ones(numel(lineorder),nstats,3);
color_ordered(lineisrun,:,:) = color(lineorder(lineisrun),:,:);

[ism,idx] = ismember(statfns,statfns_plot);
statis = find(ism);
[~,order] = sort(idx(ism));
statis = statis(order);

color_ordered = color_ordered(:,statis,:);

hfig = 2;
figure(hfig);
clf;
set(hfig,'Units','pixels','Position',[10,10,1200,1500]);

image(color_ordered);

behaviornames = regexprep(statfns_plot,'^fractime_flyany_frame','');
set(gca,'XTick',1:numel(statis),'XTickLabel',behaviornames);
set(gca,'YTick',1:nlines,'YTickLabel',line_name_order,'box','off');

htext = nan(nlines,numel(statis));
for lineii = 1:numel(lineorder),
  linei = lineorder(lineii);
  if linei == 0,
    continue;
  end
  for statii = 1:numel(statis),
    stati = statis(statii);
    if adj_pvalue(linei,stati) >= .05,
      stars = '   ';
    elseif adj_pvalue(linei,stati) >= .01,
      stars = '*  ';
    elseif adj_pvalue(linei,stati) == 0,
      stars = '***';
    else
      nstars = min(3,floor(-log10(adj_pvalue(linei,stati))));
      stars = [repmat('*',[1,nstars]),repmat(' ',[1,3-nstars])];
    end
    htext(lineii,statii) = text(statii,lineii,...
      sprintf('%.4f%s',muhat_weighted_perstat(linei,stati),stars),...
      'HorizontalAlignment','center','VerticalAlignment','middle');
    if linei == linecontrolidx,
      set(htext(lineii,statii),'FontWeight','bold');
    end
    if color_ordered(lineii,statii,3)/color_ordered(lineii,statii,1) > 1.5,
      set(htext(lineii,statii),'Color',[.99,.99,.99]);
    end
  end
end
title('Fraction of time performing each behavior');
set(gca,'Position',[.04,.02,.905,.95]);
set(gca,'TickLength',[0,0]);

% plot lines every 5 lines so that we can line things up
hline = [];
hold on;
for i = 5:5:nlines,
  hline(end+1) = plot([.5,numel(statis)+.5],[i,i]+.5,'k-');
end

hcolorbar_ax = axes('Position',[.95,.625,.045,.25]);
image(permute(flipud(cm),[1,3,2]));
axis off;
text(1,size(cm,1),'Less','Color',[.99,.99,.99],...
  'HorizontalAlignment','center','VerticalAlignment','bottom','FontName','times');
text(1,1,'More','Color','k',...
  'HorizontalAlignment','center','VerticalAlignment','top','FontName','times');
text(1,ncolors+.5,'Same','Color','k',...
  'HorizontalAlignment','center','VerticalAlignment','middle','FontName','times');

hcolorbar_ax2 = axes('Position',[.95,.125,.045,.25]);
image(permute(gray(ncolors*2),[1,3,2]));
axis off;
text(1,2*ncolors,'p\geq.05','Color','k',...
  'HorizontalAlignment','center','VerticalAlignment','bottom','FontName','times');
% text(1,1,{'***','p<.001'},'Color',[.99,.99,.99],...
%   'HorizontalAlignment','center','VerticalAlignment','top','FontName','times');
text(1,2*ncolors*(log10(.01)-minp)/(maxp-minp),{'**','p<.01'},'Color',[.99,.99,.99],...
  'HorizontalAlignment','center','VerticalAlignment','middle','FontName','times');
text(1,2*ncolors*(log10(.05)-minp)/(maxp-minp),{'*','p<.05'},'Color',[.99,.99,.99],...
  'HorizontalAlignment','center','VerticalAlignment','middle','FontName','times');

SaveFigLotsOfWays(hfig,'MBFlyBowlResultsFigure20130321Corrected');

%% output p-values to csv file

fid = fopen('MBFlyBowlResultsPValues20130321Corrected.csv','w');
fprintf(fid,'line');
fprintf(fid,',%s',statfns{statis});
fprintf(fid,'\n');
for i = 1:numel(lineorder),
  
  ii = lineorder(i);
  fprintf(fid,'%s',line_name_order{i});  
  if i ~= 0,
    fprintf(fid,',%f',sign(dcontrol(i,statis)).*adj_pvalue(i,statis));
  end
  fprintf(fid,'\n');
end
fclose(fid);

%% save

save MBConfidenceAndStrengthData20130321.mat;
rootdatadir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlData';
expdirs = importdata('/groups/branson/home/bransonk/behavioranalysis/code/Jdetect/Jdetect/experiments/datamanagement/expdirs_primary20130306.txt');
mindatenum = datenum('20130306','yyyymmdd');

allstats = struct;
statfns = {};
for i = 1:numel(expdirs),
  if mod(i,100) == 0,
    fprintf('%d / %d\n',i,numel(expdirs));
  end
  expdir = expdirs{i};
  statsfile = fullfile(expdir,'stats_perframe.mat');
  tmp = dir(statsfile);
  if isempty(tmp),
    fprintf('%s does not exist\n',statsfile);
    continue;
  elseif tmp(1).datenum < mindatenum,
    fprintf('%s is too old (%s)\n',statsfile,tmp(1).date);
    continue;
  end
  statsdata = load(statsfile,'statsperexp');
  fnscurr = fieldnames(statsdata.statsperexp);
  newfns = setdiff(fnscurr,statfns);
  for j = 1:numel(newfns),
    fn = newfns{j};
    allstats.(fn) = nan(1,numel(expdirs));
    fprintf('exp %d: Adding %s\n',i,fn);
  end
  statfns = union(statfns,newfns);
  for j = 1:numel(fnscurr),
    fn = fnscurr{j};
    allstats.(fn)(i) = statsdata.statsperexp.(fn).meanmean;
  end
  
end

% read metadata for these experiments
metadata = SAGEListBowlExperiments('checkflags',false,'rootdir',rootdatadir,'data_type',{'QuickStats_BackSubStats_meanNConnComps','ctrax_diagnostics_nframes_analyzed'},'removemissingdata',false);
[ism,idx] = ismember(expdirs,{metadata.file_system_path});
clear metadata_ordered;
for i = find(ism)',
  metadata_ordered(i) = metadata(idx(i));
end
% ism is always true
badidx = ~ism;

% correct some metadata
checkdata = load('PrimaryDataCheck20130306.mat');
for i = 1:numel(checkdata.expdirs),
  [~,name] = fileparts(checkdata.expdirs{i});
  checkdata.expdirs{i} = fullfile(rootdatadir,name);
end
[ism2,idx2] = ismember(expdirs,checkdata.expdirs);
for i = find(ism2)',
  if checkdata.isok(idx2(i)),
    metadata_ordered(i).automated_pf = 'P';
  end
end

metadata = metadata_ordered(~badidx);
if any(badidx),
  for i = 1:numel(statfns),
    allstats.(statfns{i})(badidx) = [];
  end
end


idxanalyze = strcmp({metadata.screen_type},'primary') & ...
  strcmpi({metadata.automated_pf},'P') & ...
  ~strcmpi({metadata.manual_pf},'F');

%% sanity check that everything still lines up
i = 1234;
expdir = metadata(i).file_system_path;
statsfile = fullfile(expdir,'stats_perframe.mat');
statsdata = load(statsfile,'statsperexp');
for j = 1:numel(statfns),
  fn = statfns{j};
  v1 = statsdata.statsperexp.(fn).meanmean;
  v2 = allstats.(fn)(i);
  if (~isnan(v1) || ~isnan(v2)) && v1 ~= v2,
    fprintf('%d: %s\n',j,fn);
  end
end


%% save

save CollectedPerFrameStats20130306.mat allstats metadata idxanalyze;

%% remove bad data

metadata = metadata(idxanalyze);
for i = 1:numel(statfns),
  allstats.(statfns{i})(~idxanalyze) = [];
end

%% collate by set

[set_names,idx,setidx] = unique({metadata.set});
setstats = struct;
setstats.means = struct;
setstats.stds = struct;
setstats.nexps = struct;
setstats.metadata = metadata(idx);

for i = 1:numel(statfns),
  setstats.means.(statfns{i}) = nan(1,numel(set_names));
  setstats.stds.(statfns{i}) = nan(1,numel(set_names));
  setstats.nexps.(statfns{i}) = nan(1,numel(set_names));
end
for i = 1:numel(set_names),
  idx = setidx == i;
  for j = 1:numel(statfns),
    fn = statfns{j};
    datacurr = allstats.(fn)(idx);
    datacurr = datacurr(~isnan(datacurr));
    if isempty(datacurr), continue; end
    setstats.means.(fn)(i) = mean(datacurr);
    setstats.stds.(fn)(i) = std(datacurr,1);
    setstats.nexps.(fn)(i) = numel(datacurr);
  end
end

%% collate sets by line

minnexpsperset = 2;
[line_names,~,set2lineidx] = unique({setstats.metadata.line_name});
linestats = struct;
linestats.means = struct;
linestats.stds = struct;
linestats.nsets = struct;
linestats.line_names = line_names;

for i = 1:numel(statfns),
  linestats.means.(statfns{i}) = nan(1,numel(line_names));
  linestats.stds.(statfns{i}) = nan(1,numel(line_names));
  linestats.nsets.(statfns{i}) = nan(1,numel(line_names));
end
for i = 1:numel(line_names),
  idx = find(set2lineidx == i);
  for j = 1:numel(statfns),
    fn = statfns{j};
    idxcurr = idx;
    idxcurr(setstats.nexps.(fn)(idx) < minnexpsperset) = [];
    datacurr = setstats.means.(fn)(idxcurr);
    datacurr = datacurr(~isnan(datacurr));
    if isempty(datacurr), continue; end
    linestats.means.(fn)(i) = mean(datacurr);
    linestats.stds.(fn)(i) = std(datacurr,1);
    linestats.nsets.(fn)(i) = numel(datacurr);
  end
end

%% reorganize

nstats = numel(statfns);
nlines = numel(line_names);
nsets = numel(setstats.metadata);

linemeanstats = nan(nlines,nstats);
linestdstats = nan(nlines,nstats);
linensets = nan(nlines,nstats);
setmeanstats = nan(nsets,nstats);
setstdstats = nan(nsets,nstats);
setnexps = nan(nsets,nstats);

for i = 1:nstats,
  
  fn = statfns{i};
  
  linemeanstats(:,i) = linestats.means.(fn);
  linestdstats(:,i) = linestats.stds.(fn);
  linensets(:,i) = linestats.nsets.(fn);
  
  setmeanstats(:,i) = setstats.means.(fn);
  setstdstats(:,i) = setstats.stds.(fn);
  setnexps(:,i) = setstats.nexps.(fn);
  
end

%% save

save -append CollectedPerFrameStats20130306.mat setstats linestats setidx set2lineidx linemeanstats linestdstats linensets setmeanstats setstdstats setnexps;

%% plot each statistic for all lines

%lines_show = line_names(1:20:end)';

linename_perset = {setstats.metadata.line_name};
p = [0,1,5,25,50,75,95,99,100];
setidxcontrol = strcmp(linename_perset,'pBDPGAL4U');
lineidxcontrol = find(strcmp(linestats.line_names,'pBDPGAL4U'));

setmucontrol = setmeanstats(setidxcontrol,:);
setstdcontrol = setstdstats(setidxcontrol,:);
prctiles_control = prctile(setmucontrol,p,1);
medi = find(p == 50);
thresh_low = prctile(setmucontrol,1,1);
thresh_high = prctile(setmucontrol,99,1);

npregions = (numel(p)-1)/2;
linecolors = jet(nlines)*.7;
linecolors(lineidxcontrol,:) = 0;

short_line_names = regexprep(line_names,'^GMR_','');
short_line_names = regexprep(short_line_names,'_AE_01$','');
short_line_names = regexprep(short_line_names,'_01$','');

xtick0 = 0:500:nlines;
xticklabel0 = cellstr(num2str(xtick0'));

outdir = 'PerFrameStatFigs20130307';
mkdir(outdir);

for stati = 1:nstats,

  hfig = 1000;
  if ~ishandle(hfig),
    figure(hfig);
  else
    set(0,'CurrentFigure',hfig);
  end
  clf;
  %set(hfig,'Units','pixels','Position',[110,1090,802,368]);
  set(hfig,'Units','pixels','Position',[10,10,2525,1050])
  
  
  for j = npregions:-1:1
    ycurr = prctiles_control([npregions+1-j,npregions+1+j],stati);
    patch([-4,nlines+5,nlines+5,-4,-4],ycurr([1,1,2,2,1]),diff(5+p([npregions+1-j,npregions+1+j]))/105*[1,1,1],'LineStyle','none');
    hold on;
  end
  plot([-4,nlines+5],prctiles_control(npregions+1,stati)+[0,0],'k-');
  for j = 1:numel(p),
    text(1,prctiles_control(j,stati),sprintf('%d%%ile',p(j)),'HorizontalAlignment','left');
  end
  
  [v,order] = sort(linemeanstats(:,stati));
  [~,reorder] = sort(order);
  for j = 1:nlines,
    idxcurr = set2lineidx == j;
    plot(reorder([j,j]),linemeanstats(j,stati)+[-1,1]*linestdstats(j,stati),'-','color',(linecolors(j,:)+1)/2);
    plot(reorder(j),setmeanstats(idxcurr,stati),'.','color',(linecolors(j,:)+3)/4);
    plot(reorder(j),linemeanstats(j,stati),'*','color',linecolors(j,:));
  end  
  
  minv = min(setmeanstats(:,stati));
  maxv = max(setmeanstats(:,stati));
  if minv >= 0 && minv/(maxv-minv) <= .001,
    ylim = [0,1.01*maxv];
  else
    dv = maxv-minv;
    ylim = [minv-dv*.01,maxv+dv*.01];
  end
  
  set(gca,'XLim',[-4,nlines+5],'YLim',ylim);

  ylabel(statfns{stati},'Interpreter','None');
  
  idxsig = find(linemeanstats(:,stati) > thresh_high(stati) | ...
    linemeanstats(:,stati) < thresh_low(stati));

%   set(gca,'XTick',1:nlines,'XTickLabel',short_line_names(order));
%   hxtick = rotateticklabel(gca,90);

  
  hsig = nan(1,numel(idxsig));
  for jj = 1:numel(hsig),
    j = idxsig(jj);
    shortname = short_line_names{j};
    hsig(jj) = text(reorder(j),linemeanstats(j,stati),[' ',shortname],'color',linecolors(j,:),'HorizontalAlignment','left','VerticalAlignment','middle','Interpreter','none');
  end

  if ~ismember(reorder(lineidxcontrol),xtick0),
    xticki = find(reorder(lineidxcontrol)>xtick0,1,'last');
    xtick = [xtick0(1:xticki),reorder(lineidxcontrol),xtick0(xticki+1:end)];
    xticklabel = [xticklabel0(1:xticki);{'pBDPGAL4U'};xticklabel0(xticki+1:end)];
  else
    xtick = xtick0;
    xticklabel = xticklabel0;
  end
  
  set(gca,'XTick',xtick,'XTickLabel',xticklabel);
      
  drawnow;
  
  SaveFigLotsOfWays(hfig,fullfile(outdir,sprintf('SortedLines_%s_AllGAL4s',statfns{stati})));
  
end


%%

for i = 1:nstats,
  fprintf('%d: %s\n',i,statfns{i});
end
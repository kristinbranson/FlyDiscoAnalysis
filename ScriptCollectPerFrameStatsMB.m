expdirs = importdata('/groups/branson/home/bransonk/behavioranalysis/code/Jdetect/Jdetect/experiments/mushroombody/MBData_primary_explist_20130317.txt');
mindatenum = datenum('20130307','yyyymmdd');

allstats = struct;
statfns = {};
metadata = [];
for i = 1:numel(expdirs),
  if mod(i,100) == 0,
    fprintf('%d / %d\n',i,numel(expdirs));
  end
  expdir = expdirs{i};
  
  metadatacurr = ReadMetadataFile(fullfile(expdirs{i},'Metadata.xml'));
  metadatacurr.file_system_path = expdirs{i};
  metadatacurr.automated_pf = 'F';
  metadatacurr.line_name = metadatacurr.line;

  metadata = structappend(metadata,metadatacurr);
  
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
  
  metadata(end).automated_pf = 'P';
  
end

idxanalyze = strcmp({metadata.screen_type},'non_olympiad_mushroombody') & ...
  strcmpi({metadata.automated_pf},'P') & ...
  ([metadata.flag_redo] ~= 1) & ...
  ([metadata.flag_aborted] ~= 1);


%% add "set" -- super-experiment

MAX_SET_TIMERANGE = 10/(24*60);
MAX_EXPS_PER_SET = 4;

[line_names,~,lineidx] = unique({metadata.line_name});
sets = nan(1,numel(metadata));
seti = 0;
for linei = 1:numel(line_names),
  expidx1 = find(lineidx==linei);
  [rigs,~,rigidx] = unique([metadata(expidx1).rig]);
  for rigi = 1:numel(rigs),
    expidx2 = expidx1(rigidx==rigi);

    % sort by datetime
    [exp_datenum,order] = sort(datenum({metadata(expidx2).exp_datetime},'yyyymmddTHHMMSS'));
    expidx2 = expidx2(order);
    min_set_time = 0;
    nperset = 0;
  
    for i = 1:numel(expidx2),
      
      % start a new set?
      if exp_datenum(i)-min_set_time > MAX_SET_TIMERANGE || nperset >= MAX_EXPS_PER_SET,
        min_set_time = inf;
        seti = seti+1;
        nperset = 0;
      end
      sets(expidx2(i)) = seti;
      nperset = nperset+1;
      min_set_time = min(min_set_time,exp_datenum(i));
      
    end
  end
end

for seti = 1:max(sets),
  expidx = find(sets==seti);
  %fprintf('Experiments in set %d:\n',seti);
  %fprintf('%s\n',metadata(expidx).experiment_name);
  [~,order] = sort({metadata(expidx).exp_datetime});
  min_datetime = metadata(expidx(order(1))).exp_datetime;
  set_name = sprintf('%s__Rig%d__%s',metadata(expidx(1)).line_name,...
    metadata(expidx(1)).rig,...
    min_datetime);
  for i = expidx(:)',
    metadata(i).set = set_name;
  end
end

%% save

save CollectedPerFrameStatsMB20130317.mat allstats metadata idxanalyze;

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

save -append CollectedPerFrameStatsMB20130317.mat setstats linestats setidx set2lineidx linemeanstats linestdstats linensets setmeanstats setstdstats setnexps;

%% plot each statistic for all lines

%lines_show = line_names(1:20:end)';

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
  'fractime_flyany_frameattemptedcopulation'
  'velmag_ctr_flyany_framewalk'};



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
short_line_names = regexprep(short_line_names,'_BX3G_1500154$','');

xtick0 = 0:20:nlines;
xticklabel0 = cellstr(num2str(xtick0'));

outdir = 'MBPerFrameStatFigs20130317';
mkdir(outdir);

% statis = find(~cellfun(@isempty,regexp(statfns,'^fractime','once')) | ...
%   ~cellfun(@isempty,regexp(statfns,'^velmag','once')));
statis = find(ismember(statfns,statfns_plot));

[~,lineidx] = ismember({metadata.line_name},line_names);

npad = 1;

for stati = statis',

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
    patch([-npad+1,nlines+npad,nlines+npad,-npad+1,-npad+1],ycurr([1,1,2,2,1]),diff(5+p([npregions+1-j,npregions+1+j]))/105*[1,1,1],'LineStyle','none');
    hold on;
  end
  plot([-npad+1,nlines+npad],prctiles_control(npregions+1,stati)+[0,0],'k-');
  for j = 1:numel(p),
    text(1,prctiles_control(j,stati),sprintf('%d%%ile',p(j)),'HorizontalAlignment','left');
  end
  
  isnodata = isnan(linemeanstats(:,stati));
  nbadlines = nnz(isnodata);
  [v,order] = sort(linemeanstats(:,stati));
  [~,reorder] = sort(order);
  for j = 1:nlines,
    if isnodata(j),
      continue;
    end
    idxcurr = find(set2lineidx == j);
%     for k = 1:numel(idxcurr),
%       idxexp = setidx == idxcurr(k);
%       if numel(idxcurr) == 1,
%         off = 0;
%       else
%         off = (k-1)/(numel(idxcurr)-1)*.75-.75/2;
%       end
%       plot(reorder(j)+off,allstats.(statfns{stati})(idxexp),'x','color',(linecolors(j,:)+1)/2);
%     end
    plot(reorder([j,j]),linemeanstats(j,stati)+[-1,1]*linestdstats(j,stati),'-','color',(linecolors(j,:)+1)/2);
    plot(reorder(j),setmeanstats(idxcurr,stati),'.','color',(linecolors(j,:)+1)/2);
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
   
  set(gca,'XLim',[-npad+1,nlines-nbadlines+npad],'YLim',ylim);

  ylabel(statfns{stati},'Interpreter','None');
  
  idxsig = find(linemeanstats(:,stati) > thresh_high(stati) | ...
    linemeanstats(:,stati) < thresh_low(stati));

%   set(gca,'XTick',1:nlines,'XTickLabel',short_line_names(order));
%   hxtick = rotateticklabel(gca,90);

  
%   hsig = nan(1,numel(idxsig));
%   for jj = 1:numel(hsig),
%     j = idxsig(jj);
%     shortname = short_line_names{j};
%     hsig(jj) = text(reorder(j),linemeanstats(j,stati),[' ',shortname],'color',linecolors(j,:),'HorizontalAlignment','left','VerticalAlignment','middle','Interpreter','none');
%   end
  
  set(gca,'XTick',1:nlines-nbadlines,'XTickLabel',short_line_names(order(1:nlines-nbadlines)));
  htick = rotateticklabel(gca);    
  for i = 1:nlines-nbadlines,
    set(htick(i),'Color',linecolors(order(i),:));
  end
  
  drawnow;
  
  SaveFigLotsOfWays(hfig,fullfile(outdir,sprintf('SortedLines_%s_MB',statfns{stati})));
  
end


%%

for i = 1:nstats,
  fprintf('%d: %s\n',i,statfns{i});
end
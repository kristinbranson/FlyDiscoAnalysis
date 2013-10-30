addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;

GTDataDir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/GTFlyBowlPaper';

%% parameters

behaviors = {'Jump','WingExtension'};
nbehaviors = numel(behaviors);

orderedtests = {
  'behavior hits 0th','0','Percentile of hits'
  'FP 5th','5','Percentile of hits'
  'FP 15th','15','Percentile of hits'
  'FP 50th','50','Percentile of hits'
  'FP 100th','100','Percentile of hits'
  'control','Ctrl',''
  'cluster 1','1','Cluster representative'
  'cluster 2','2','Cluster representative'
  'cluster 3','3','Cluster representative'
  'cluster 4','4','Cluster representative'
  'cluster 5','5','Cluster representative'
  'cluster 6','6','Cluster representative'
  'cluster 7','7','Cluster representative'
  'cluster 8','8','Cluster representative'
  'cluster 9','9','Cluster representative'
  'cluster 10','10','Cluster representative'
  'not any behavior 0th','0','%ile no-behavior hits'
  'NAB 5th','5','%ile no-behavior hits'
  'NAB 15th','15','%ile no-behavior hits'
  'NAB 50th','50','%ile no-behavior hits'
  'NAB 100th','100','%ile no-behavior hits'
  'uncertain scores 1','1','Low confidence'
  'unceratin scores 2','2','Low confidence'
};
ntests = size(orderedtests,1);

%% load in results

gtdata = struct;
for i = 1:nbehaviors,
  filename = fullfile(GTDataDir,sprintf('GTResults_%s.csv',behaviors{i}));
  fid = fopen(filename,'r');
  gtdatacurr = {};
  while true,
    s = fgetl(fid);
    if ~ischar(s),
      break;
    end
    if isempty(s),
      continue;
    end;
    ss = regexp(s,',','split');
    ss = strrep(ss,'"','');
    gtdatacurr(end+1,1:numel(ss)) = ss;
  end
  fclose(fid);
  
  switch behaviors{i},
    case 'Jump',
      gtfns = struct('line_name','lines',...
        'test','test',...
        'jabfile','jab file name',...
        'ntruepos_imp','j imp = j',...
        'nfalseneg_imp','J imp - N',...
        'nfalsepos_imp', 'N imp = J',...
        'ntrueneg_imp','N imp = N',...
        'ntruepos_all','j = j',...
        'nfalseneg_all','j = n',...
        'nfalsepos_all', 'n = j',...
        'ntrueneg_all','n = n',...
        'boutsize','bout size');
    case 'WingExtension',
      gtfns = struct('line_name','lines',...
        'test','test',...
        'jabfile','',...
        'ntruepos_imp','b imp = b',...
        'nfalseneg_imp','b imp = n',...
        'nfalsepos_imp', 'N imp = b',...
        'ntrueneg_imp','N imp = N',...
        'ntruepos_all','b = b',...
        'nfalseneg_all','b = n',...
        'nfalsepos_all', 'n = b',...
        'ntrueneg_all','n = n',...
        'boutsize','bout size');
  end
  
  fns = fieldnames(gtfns);
  for j = 1:numel(fns),
    fn = gtfns.(fns{j});
    k = find(strcmpi(fn,gtdatacurr(1,:)));
    if isempty(k),
      error('Could not find column header %s',fn);
    end
    s = str2double(gtdatacurr(2:end,k));
    if ~any(isnan(s)),
      gtdata(i).(fns{j}) = s;
    else
      gtdata(i).(fns{j}) = gtdatacurr(2:end,k);
    end
  end
  
  gtdata(i).falseposrate_imp = gtdata(i).nfalsepos_imp ./ ...
    max(1,(gtdata(i).nfalsepos_imp + gtdata(i).ntrueneg_imp));
  gtdata(i).falsenegrate_imp = gtdata(i).nfalseneg_imp ./ ...
    max(1,(gtdata(i).nfalseneg_imp + gtdata(i).ntruepos_imp));

  gtdata(i).falseposrate_all = gtdata(i).nfalsepos_all ./ ...
    max(1,(gtdata(i).nfalsepos_all + gtdata(i).ntrueneg_all));
  gtdata(i).falsenegrate_all = gtdata(i).nfalseneg_all ./ ...
    max(1,(gtdata(i).nfalseneg_all + gtdata(i).ntruepos_all));
  
end

fns = {
  'ntruepos_imp'
  'nfalseneg_imp'
  'nfalsepos_imp'
  'ntrueneg_imp'
  'ntruepos_all'
  'nfalseneg_all'
  'nfalsepos_all'
  'ntrueneg_all'
  };

% compute the total
gttotal = struct;
for j = 1:nbehaviors,
  for i = 1:numel(fns),
    gttotal(j).(fns{i}) = sum(gtdata(j).(fns{i}));
  end
  gttotal(j).falseposrate_imp = gttotal(j).nfalsepos_imp ./ ...
    max(1,(gttotal(j).nfalsepos_imp + gttotal(j).ntrueneg_imp));
  gttotal(j).falsenegrate_imp = gttotal(j).nfalseneg_imp ./ ...
    max(1,(gttotal(j).nfalseneg_imp + gttotal(j).ntruepos_imp));
  
  gttotal(j).falseposrate_all = gttotal(j).nfalsepos_all ./ ...
    max(1,(gttotal(j).nfalsepos_all + gttotal(j).ntrueneg_all));
  gttotal(j).falsenegrate_all = gttotal(j).nfalseneg_all ./ ...
    max(1,(gttotal(j).nfalseneg_all + gttotal(j).ntruepos_all));
end

nlinesgt = numel(gtdata(1).test);

%% plot on one axes

skipsome = [false;~strcmp(orderedtests(1:end-1,3),orderedtests(2:end,3))];

hfig = 2;
figure(hfig);
clf;
set(hfig,'Units','pixels','Position',[10, 10 1250 450]);
hax = createsubplots(1,2,[.075,.03;.2,.015]);
hax = reshape(hax,[1,nbehaviors]);

ntrueneg = [gtdata.nfalsepos_imp] + [gtdata.ntrueneg_imp];
ntruepos = [gtdata.nfalseneg_imp] + [gtdata.ntruepos_imp];
falseposrate = [gtdata.falseposrate_imp];
falsenegrate = [gtdata.falsenegrate_imp];

minntrueneg = min(ntrueneg(:));
maxntrueneg = max(ntrueneg(:));

minntruepos = min(ntruepos(:));
maxntruepos = max(ntruepos(:));
maxn = max(maxntruepos,maxntrueneg);

minfalseposrate = min(falseposrate(:));
maxfalseposrate = max(falseposrate(:));

minfalsenegrate = min(falsenegrate(:));
maxfalsenegrate = max(falsenegrate(:));

x0 = cumsum(skipsome*.5+1);
xpos = [x0'-.4;x0'-.025];  
xneg = [x0'+.025;x0'+.4];

ncolors = 256;
w = linspace(0,1,ncolors)';
red = [.7,0,0];
blue = [0,0,.7];
redmap = bsxfun(@times,[1,1,1],1-w) + bsxfun(@times,red,w);
bluemap = bsxfun(@times,[1,1,1],1-w) + bsxfun(@times,blue,w);

h = [];

for i = 1:nbehaviors,
  [ism,order] = ismember(gtdata(i).test,orderedtests(:,1));
  if ~all(ism),
    error('Could not find all tests');
  end
  
  axes(hax(i));
  y = [zeros(ntests,1),...
    gtdata(i).falseposrate_imp(order)]';
  n = gtdata(i).nfalsepos_imp(order)' + ...
    gtdata(i).ntrueneg_imp(order)';
  s = cell(size(n));
  for j = 1:numel(n),
    s{j} = sprintf(' %d',n(j));
  end
  ni = min(ncolors,floor(n(:)'/maxn*ncolors)+1);
  c = bluemap(ni,:);
  for j = 1:size(xpos,2),
    h(1,j) = patch(xpos([1,1,2,2,1],j),y([1,2,2,1,1],j),c(j,:),'EdgeColor',blue);
  end
  %text(mean(xpos,1),gtdata(i).falseposrate_imp(order),s,...
  %  'HorizontalAlignment','left','VerticalAlignment','middle','Rotation',90);
  
  y = [zeros(ntests,1),...
    gtdata(i).falsenegrate_imp(order)]';
  n = gtdata(i).nfalseneg_imp(order)' + ...
    gtdata(i).ntruepos_imp(order)';
  ni = min(ncolors,floor(n(:)'/maxn*ncolors)+1);
  c = redmap(ni,:);
  for j = 1:size(xneg,2),
    h(2,j) = patch(xneg([1,1,2,2,1],j),y([1,2,2,1,1],j),c(j,:),'EdgeColor',red);
  end
  s = cell(size(n));
  for j = 1:numel(n),
    s{j} = sprintf(' %d',n(j));
  end
  %text(mean(xneg,1),gtdata(i).falsenegrate_imp(order),s,...
  %  'HorizontalAlignment','left','VerticalAlignment','middle','Rotation',90);

  y = [0;gttotal(i).falseposrate_imp];
  patch(xpos([1,1,2,2,1],end)+2,y([1,2,2,1,1]),[0,0,.3],'CDataMapping','direct');
  y = [0;gttotal(i).falsenegrate_imp];
  patch(xneg([1,1,2,2,1],end)+2,y([1,2,2,1,1]),[.3,0,0],'CDataMapping','direct');

  
  if i == 1,
    ylabel('Error rate');
    hleg = legend(h(1:2,1),{'False positive rate','False negative rate'});
    set(hleg,'Box','off');
  end
  title(behaviors{i});
    
  %datacurr = [gtdata(i).nfalseneg_imp(order),...
  %  gtdata(i).ntruepos_imp(order)];
  %bar(1:ntests,datacurr);
  
end

maxerrorrate = max(maxfalseposrate,maxfalsenegrate);
set(hax,'Ylim',[-.01*maxerrorrate,maxerrorrate*1.05],'XLim',[0,max(x0)+3]);

set(hax,'XTick',[x0;x0(end)+2],'XTickLabel',[orderedtests(:,2)',{'Total'}]);
set(hax(2:end),'YTickLabel',{});

[orderedtesttypes,~,testtypeidx] = unique(orderedtests(:,3));
harrow = nan(numel(orderedtesttypes),nbehaviors);
htext = nan(numel(orderedtesttypes),nbehaviors);
for i = 1:numel(orderedtesttypes),
  idxcurr = find(testtypeidx == i);
  if numel(idxcurr) <= 1,
    continue;
  end
  x1 = xpos(1,idxcurr(1));
  x2 = xneg(2,idxcurr(end));
  for j = 1:nbehaviors,
    xfig = dsxy2figxy(hax(j),[x1,x2],[0,0]);
    pos = get(hax(j),'Position');
    yfig1 = pos(2)*.65 + [0,0];
    yfig2 = pos(2)*[.4,.65];

    harrow(i,j) = annotation('doublearrow',xfig,yfig1,...
      'Head1Style','rectangle','Head2Style','rectangle',...
      'Head1Length',1,'Head2Length',1);
    htext(i,j) = annotation('textbox',[xfig(1),yfig2(1),diff(xfig),diff(yfig2)],...
      'String',orderedtesttypes{i},'LineStyle','none','HorizontalAlignment','center');
      
  end
end

hfigcb = 3;
figure(hfigcb);
clf;
set(hfigcb,'Units','pixels','Position',[10, 10 1250 450]);
hax = createsubplots(1,1,[.075,.03;.2,.015]);
imagesc(0);
set(gca,'CLim',[0,maxn]);
colormap(redmap);
hcb = colorbar;
ytick = get(hcb,'YTick');
set(hcb,'YTick',[0,1000,2000]);

hfigcb2 = 4;
figure(hfigcb2);
clf;
set(hfigcb2,'Units','pixels','Position',[10, 10 1250 450]);
hax = createsubplots(1,1,[.075,.03;.2,.015]);
imagesc(0);
set(gca,'CLim',[0,maxn]);
colormap(bluemap);
hcb = colorbar;
ytick = get(hcb,'YTick');
set(hcb,'YTick',[0,1000,2000]);

for i = 1:nbehaviors,
  beh = behaviors{i};
  fprintf('%s\n',beh);
  fprintf('Minimum false positive rate: %f\n',min(gtdata(i).falseposrate_imp));
  fprintf('Maximum false positive rate: %f\n',max(gtdata(i).falseposrate_imp));
  fprintf('Median false positive rate: %f\n',median(gtdata(i).falseposrate_imp));
  fprintf('Mean false positive rate: %f\n',mean(gtdata(i).falseposrate_imp));
  fprintf('Minimum false negative rate: %f\n',min(gtdata(i).falsenegrate_imp));
  fprintf('Maximum false negative rate: %f\n',max(gtdata(i).falsenegrate_imp));
  fprintf('Median false negative rate: %f\n',median(gtdata(i).falsenegrate_imp));
  fprintf('Mean false negative rate: %f\n',mean(gtdata(i).falsenegrate_imp));
  fprintf('Minimum error rate: %f\n',min(gtdata(i).falseposrate_imp+gtdata(i).falsenegrate_imp)/2);
  fprintf('Maximum error rate: %f\n',max(gtdata(i).falseposrate_imp+gtdata(i).falsenegrate_imp)/2);
  fprintf('Median error rate: %f\n',median(gtdata(i).falseposrate_imp+gtdata(i).falsenegrate_imp)/2);
  fprintf('Mean error rate: %f\n',mean(gtdata(i).falseposrate_imp+gtdata(i).falsenegrate_imp)/2);
end

% Jump
% Minimum false positive rate: 0.000000
% Maximum false positive rate: 0.049370
% Median false positive rate: 0.007246
% Mean false positive rate: 0.012426
% Minimum false negative rate: 0.000000
% Maximum false negative rate: 0.047619
% Median false negative rate: 0.000000
% Mean false negative rate: 0.008426
% Minimum error rate: 0.000901
% Maximum error rate: 0.037520
% Median error rate: 0.005533
% Mean error rate: 0.010426
% WingExtension
% Minimum false positive rate: 0.000000
% Maximum false positive rate: 0.111903
% Median false positive rate: 0.024639
% Mean false positive rate: 0.028013
% Minimum false negative rate: 0.000000
% Maximum false negative rate: 0.322581
% Median false negative rate: 0.044068
% Mean false negative rate: 0.058952
% Minimum error rate: 0.000000
% Maximum error rate: 0.175284
% Median error rate: 0.039342
% Mean error rate: 0.043483

%% save results 

SaveFigLotsOfWays(hfig,'GTFlyBowlPaper/GTResultsJumpWingExt20130925');
SaveFigLotsOfWays(hfigcb,'GTFlyBowlPaper/GTResultsJumpWingExt20130925ColorbarRed');
SaveFigLotsOfWays(hfigcb2,'GTFlyBowlPaper/GTResultsJumpWingExt20130925ColorbarBlue');

%% compute fraction of time estimates using importance sampling

frameweights = cell(1,nbehaviors);
for i = 1:nbehaviors,
  frameweightscurr = cell(1,nlinesgt);
  
  parfor j = 1:nlinesgt,
    
    fprintf('%d / %d, %d / %d\n',i,nbehaviors,j,nlinesgt);
    
    jabfile = gtdata(i).jabfile{j};
    intsize = gtdata(i).boutsize(j);
    frameweightscurr{j} = gtSuggestionWeights(jabfile,intsize);
    
  end
  
  frameweights{i} = frameweightscurr;
end

%% sanity check: if we re-weight the examples selected using these
% 1/weights, then we should get fractime simimlar to what is reported by
% the scores

fractimepred = nan(nbehaviors,nlinesgt);

for i = 1:nbehaviors,
  for j = 1:nlinesgt,
    
    jabfile = gtdata(i).jabfile{j};
    
    jabdata = load(jabfile,'-mat');
    scorefilename = jabdata.x.file.scorefilename;
    expdirs = jabdata.x.gtExpDirNames;
    weightpredpos = 0;
    weightpredneg = 0;
    
    labels = jabdata.x.gtLabels;
    
    for k = 1:numel(expdirs),
      scorefile = fullfile(expdirs{k},scorefilename);
      
      if ~exist(scorefile,'file'),
        sd = struct;
        fprintf('Computing scores for %s...\n',scorefile);
        [classifierinfo,sd.allScores] = JAABADetect(expdirs{k},'jabfiles',jabfile);
        sd.allScores = sd.allScores{1};
      else
        sd = load(scorefile);
      end
      
      for flyi = 1:numel(labels(k).t0s),
        
        fly = labels(k).flies(flyi);
        
        labelcurr = zeros(size(sd.allScores.postprocessed{fly}));
        for l = 1:numel(labels(k).t0s{flyi}),
          if strcmpi(labels(k).names{flyi}{l},'none'),
            labelcurr(labels(k).t0s{flyi}(l):labels(k).t1s{flyi}(l)-1) = -1;
          else
            labelcurr(labels(k).t0s{flyi}(l):labels(k).t1s{flyi}(l)-1) = 1;
          end
        end
        tpredpos = find(sd.allScores.postprocessed{fly} == 1 & labelcurr ~= 0);
        idx = frameweights{i}{j}.exp == k & ...
          frameweights{i}{j}.flies == fly & ...
          ismember(frameweights{i}{j}.t,tpredpos);
        weightpredpos = weightpredpos + sum(1./frameweights{i}{j}.wt(idx));
        
        tpredneg = find(sd.allScores.postprocessed{fly} == 0 & labelcurr ~= 0);
        idx = frameweights{i}{j}.exp == k & ...
          frameweights{i}{j}.flies == fly & ...
          ismember(frameweights{i}{j}.t,tpredneg);
        weightpredneg = weightpredneg + sum(1./frameweights{i}{j}.wt(idx));
        
        
      end
      
    end
    
    fractimepred(i,j) = weightpredpos/(weightpredpos+weightpredneg);
    
    if numel(gtdata(i).line_name{j}) <= 5,
      line_name = sprintf('GMR_%s_AE_01',gtdata(i).line_name{j});
      if ~any(strcmp({metadata.line_name},line_name)),
        line_name = sprintf('GMR_%s_AD_01',gtdata(i).line_name{j});
      end
    elseif strcmp(gtdata(i).line_name{j},'pBDPGAL4'),
      line_name = 'pBDPGAL4U';
    else
      line_name = gtdata(i).line_name{j};
    end
    idx = find(strcmp({metadata.line_name},line_name));

    switch behaviors{i},
      case 'Jump',
        fn = 'fractime_flyany_framejump';
      case 'WingExtension'
        fn = 'fractime_flyany_framewingextension';
      otherwise
        error('unknown behavior %s',behaviors{i});
    end
    fprintf('%s, %s: estimated fractime = %f, mean frac time = %f, per-experiment fractimes = %s\n',...
      behaviors{i},line_name,fractimepred(i,j),...
      linestats.means.(fn)(strcmp(linestats.line_names,line_name)),...
      mat2str(sort(allstats.(fn)(idx))));
  end
end

%% estimate frac time statistics after resampling

fractimelabel = nan(nbehaviors,nlinesgt);
fractimelabel_perexp = cell(nbehaviors,nlinesgt);
fractimelabel_std = nan(nbehaviors,nlinesgt);

for i = 1:nbehaviors,
  for j = 1:nlinesgt,
    
    jabfile = gtdata(i).jabfile{j};
    
    jabdata = load(jabfile,'-mat');
    scorefilename = jabdata.x.file.scorefilename;
    expdirs = jabdata.x.gtExpDirNames;
    weightlabelpos = zeros(1,numel(expdirs));
    weightlabelneg = zeros(1,numel(expdirs));
    
    labels = jabdata.x.gtLabels;
    
    for k = 1:numel(expdirs),
%       scorefile = fullfile(expdirs{k},scorefilename);
%       
%       if ~exist(scorefile,'file'),
%         sd = struct;
%         fprintf('Computing scores for %s...\n',scorefile);
%         [classifierinfo,sd.allScores] = JAABADetect(expdirs{k},'jabfiles',jabfile);
%         sd.allScores = sd.allScores{1};
%       else
%         sd = load(scorefile);
%       end
      
      for flyi = 1:numel(labels(k).t0s),
        
        if isempty(labels(k).t1s{flyi}),
          continue;
        end
        
        fly = labels(k).flies(flyi);
        
        labelcurr = zeros(1,labels(k).t1s{flyi}(end));
        for l = 1:numel(labels(k).t0s{flyi}),
          if strcmpi(labels(k).names{flyi}{l},'none'),
            labelcurr(labels(k).t0s{flyi}(l):labels(k).t1s{flyi}(l)-1) = -1;
          else
            labelcurr(labels(k).t0s{flyi}(l):labels(k).t1s{flyi}(l)-1) = 1;
          end
        end
        tlabelpos = find(labelcurr == 1);
        idx = frameweights{i}{j}.exp == k & ...
          frameweights{i}{j}.flies == fly & ...
          ismember(frameweights{i}{j}.t,tlabelpos);
        weightlabelpos(k) = weightlabelpos(k) + sum(1./frameweights{i}{j}.wt(idx));
        
        tlabelneg = find(labelcurr == -1);
        idx = frameweights{i}{j}.exp == k & ...
          frameweights{i}{j}.flies == fly & ...
          ismember(frameweights{i}{j}.t,tlabelneg);
        weightlabelneg(k) = weightlabelneg(k) + sum(1./frameweights{i}{j}.wt(idx));
        
        
      end
      
    end
    
    fractimelabel(i,j) = sum(weightlabelpos)/sum(weightlabelpos+weightlabelneg);
    fractimelabel_perexp{i,j} = weightlabelpos./(weightlabelpos+weightlabelneg);
    fractimelabel_std(i,j) = sqrt( nanmean( (fractimelabel_perexp{i,j}-fractimelabel(i,j)).^2 ) );
    
    if numel(gtdata(i).line_name{j}) <= 5,
      line_name = sprintf('GMR_%s_AE_01',gtdata(i).line_name{j});
      if ~any(strcmp({metadata.line_name},line_name)),
        line_name = sprintf('GMR_%s_AD_01',gtdata(i).line_name{j});
      end
    elseif strcmp(gtdata(i).line_name{j},'pBDPGAL4'),
      line_name = 'pBDPGAL4U';
    else
      line_name = gtdata(i).line_name{j};
    end
    idx = find(strcmp({metadata.line_name},line_name));

    switch behaviors{i},
      case 'Jump',
        fn = 'fractime_flyany_framejump';
      case 'WingExtension'
        fn = 'fractime_flyany_framewingextension';
      otherwise
        error('unknown behavior %s',behaviors{i});
    end
    fprintf('%s, %s: estimated fractime = %f, std = %f, per-exp = %s\nline mean frac time = %f, per-experiment fractimes = %s\n',...
      behaviors{i},line_name,fractimelabel(i,j),...
      fractimelabel_std(i,j),mat2str(fractimelabel_perexp{i,j}),...
      linestats.means.(fn)(strcmp(linestats.line_names,line_name)),...
      mat2str(sort(allstats.(fn)(idx))));
  end
end

%% estimate per-frame statistics after resampling

statfns_est = struct;
statfns_est.Jump = struct;
statfns_est.WingExtension = struct;
statfns_est.Jump(1).statfn = 'velmag_ctr';
statfns_est.Jump(1).sex = '.';
statfns_est.WingExtension(1).statfn = 'absthetadiff_nose2ell';
statfns_est.WingExtension(1).sex = 'M';
statfns_est.WingExtension(2).statfn = 'angleonclosestfly';
statfns_est.WingExtension(2).sex = 'M';
statfns_est.WingExtension(3).statfn = 'max_wing_angle';
statfns_est.WingExtension(3).sex = '.';
statfns_est.WingExtension(4).statfn = 'max_absdwing_angle';
statfns_est.WingExtension(4).sex = '.';
nstatsest = structfun(@numel,statfns_est);

statlabel = struct;
statlabel_perexp = struct;
statlabel_std = struct;

for i = 1:nbehaviors,

  behavior = behaviors{i};
  
  statlabel.(behavior) = nan(nstatsest(i),nlinesgt);
  statlabel_perexp.(behavior) = cell(nstatsest(i),nlinesgt);
  statlabel_std.(behavior) = nan(nstatsest(i),nlinesgt);

  
  for j = 1:nlinesgt,
    
    jabfile = gtdata(i).jabfile{j};
    
    jabdata = load(jabfile,'-mat');
    scorefilename = jabdata.x.file.scorefilename;
    expdirs = jabdata.x.gtExpDirNames;    
    labels = jabdata.x.gtLabels;

    wperexp = zeros(nstatsest(i),numel(expdirs));
    for stati = 1:nstatsest(i),
      statlabel_perexp.(behavior){stati,j} = zeros(1,numel(expdirs));
    end

    
    for k = 1:numel(expdirs),
%       scorefile = fullfile(expdirs{k},scorefilename);
%       
%       if ~exist(scorefile,'file'),
%         sd = struct;
%         fprintf('Computing scores for %s...\n',scorefile);
%         [classifierinfo,sd.allScores] = JAABADetect(expdirs{k},'jabfiles',jabfile);
%         sd.allScores = sd.allScores{1};
%       else
%         sd = load(scorefile);
%       end

      td = load(fullfile(expdirs{k},'registered_trx.mat'));

      for flyi = 1:numel(labels(k).t0s),
        
        if isempty(labels(k).t1s{flyi}),
          continue;
        end
        
        fly = labels(k).flies(flyi);
        
        labelcurr = zeros(1,labels(k).t1s{flyi}(end));
        for l = 1:numel(labels(k).t0s{flyi}),
          if strcmpi(labels(k).names{flyi}{l},'none'),
            labelcurr(labels(k).t0s{flyi}(l):labels(k).t1s{flyi}(l)-1) = -1;
          else
            labelcurr(labels(k).t0s{flyi}(l):labels(k).t1s{flyi}(l)-1) = 1;
          end
        end
        tlabelpos = find(labelcurr == 1);
        if isempty(tlabelpos),
          continue;
        end

        idxexpfly = find(frameweights{i}{j}.exp == k & ...
          frameweights{i}{j}.flies == fly);
        [ism,idx] = ismember(frameweights{i}{j}.t(idxexpfly),tlabelpos);
        idxcurr = idxexpfly(idx(ism));
        ilabelpos = tlabelpos - td.trx(fly).firstframe + 1;
        wt = frameweights{i}{j}.wt(idxcurr);
                
        for stati = 1:nstatsest(i),
          statfn = statfns_est.(behavior)(stati);
          
          ilabelpos = tlabelpos - td.trx(fly).firstframe + 1;
          idxsex = cellfun(@isempty,regexp(td.trx(fly).sex(ilabelpos),statfn.sex,'once'));
          tlabelpos(idxsex) = [];
          if isempty(tlabelpos),
            continue;
          end          
          ilabelpos = tlabelpos - td.trx(fly).firstframe + 1;

          % load in per-frame stat
          pfd = load(fullfile(expdirs{k},'perframe',[statfn.statfn,'.mat']));
          
          statlabel_perexp.(behavior){stati,j}(k) = statlabel_perexp.(behavior){stati,j}(k) + ...
            sum(wt.*pfd.data{fly}(ilabelpos));
          wperexp(stati,k) = wperexp(stati,k) + sum(wt);
        end
      end
      
    end

    for stati = 1:nstatsest(i),
      statlabel.(behavior)(stati,j) = sum(statlabel_perexp.(behavior){stati,j})./sum(wperexp(stati,:));
      statlabel_perexp.(behavior){stati,j} = statlabel_perexp.(behavior){stati,j}./wperexp(stati,:);
      statlabel_std.(behavior)(stati,j) = sqrt( nanmean( (statlabel_perexp.(behavior){stati,j}-statlabel.(behavior)(stati,j)).^2 ) );
    end
    
    if numel(gtdata(i).line_name{j}) <= 5,
      line_name = sprintf('GMR_%s_AE_01',gtdata(i).line_name{j});
      if ~any(strcmp({metadata.line_name},line_name)),
        line_name = sprintf('GMR_%s_AD_01',gtdata(i).line_name{j});
      end
    elseif strcmp(gtdata(i).line_name{j},'pBDPGAL4'),
      line_name = 'pBDPGAL4U';
    else
      line_name = gtdata(i).line_name{j};
    end
    idx = find(strcmp({metadata.line_name},line_name));
    
    for stati = 1:nstatsest(i),
      statfn = statfns_est.(behavior)(stati);    
      switch statfn.sex,
        case 'M',
          flys = 'male';
        case 'F',
          flys = 'female';
        otherwise
          flys = 'any';
      end
      switch behaviors{i},
        case 'Jump',
          frames = 'jump';
        case 'WingExtension'
          frames = 'wingextension';
      end
      
      fn = sprintf('%s_fly%s_frame%s',statfn.statfn,flys,frames);
      fprintf('%s, %s: estimated stat = %f, std = %f, per-exp = %s\nline mean stat = %f, per-experiment stats = %s\n',...
        fn,line_name,statlabel.(behavior)(stati,j),...
        statlabel_std.(behavior)(stati,j),mat2str(statlabel_perexp.(behavior){stati,j}),...
        linestats.means.(fn)(strcmp(linestats.line_names,line_name)),...
        mat2str(sort(allstats.(fn)(idx))));
    end
  end
end


%% collect predicted statistics

fractimepred = nan(nbehaviors,nlinesgt);
fractimepred_perexp = cell(nbehaviors,nlinesgt);
fractimepred_std = nan(nbehaviors,nlinesgt);

for i = 1:nbehaviors,
  for j = 1:nlinesgt,
    
    jabfile = gtdata(i).jabfile{j};    
    jabdata = load(jabfile,'-mat');
    expdirs = jabdata.x.gtExpDirNames;
    expnames = cell(size(expdirs));
    for k = 1:numel(expdirs),
      [~,expnames{k}] = fileparts(expdirs{k});
      expnames{k} = ['FlyBowl_',expnames{k}];
    end
    
    if numel(gtdata(i).line_name{j}) <= 5,
      line_name = sprintf('GMR_%s_AE_01',gtdata(i).line_name{j});
      if ~any(strcmp({metadata.line_name},line_name)),
        line_name = sprintf('GMR_%s_AD_01',gtdata(i).line_name{j});
      end
    elseif strcmp(gtdata(i).line_name{j},'pBDPGAL4'),
      line_name = 'pBDPGAL4U';
    else
      line_name = gtdata(i).line_name{j};
    end

    switch behaviors{i},
      case 'Jump',
        fn = 'fractime_flyany_framejump';
      case 'WingExtension'
        fn = 'fractime_flyany_framewingextension';
      otherwise
        error('unknown behavior %s',behaviors{i});
    end
    
    fractimepred(i,j) = linestats.means.(fn)(strcmp(linestats.line_names,line_name));
    [ism,idx] = ismember(expnames,{metadata.experiment_name});
    fractimepred_perexp{i,j} = allstats.(fn)(idx);
    fractimepred_std(i,j) = sqrt( nanmean( (fractimepred_perexp{i,j}-fractimepred(i,j)).^2 ) );
  end
end


statpred = struct;
statpred_perexp = struct;
statpred_std = struct;

for i = 1:nbehaviors,

  behavior = behaviors{i};
  
  statpred.(behavior) = nan(nstatsest(i),nlinesgt);
  statpred_perexp.(behavior) = cell(nstatsest(i),nlinesgt);
  statpred_std.(behavior) = nan(nstatsest(i),nlinesgt);

  
  for j = 1:nlinesgt,
    
    jabfile = gtdata(i).jabfile{j};
    
    jabdata = load(jabfile,'-mat');
    scorefilename = jabdata.x.file.scorefilename;
    expdirs = jabdata.x.gtExpDirNames;    

    expnames = cell(size(expdirs));
    for k = 1:numel(expdirs),
      [~,expnames{k}] = fileparts(expdirs{k});
      expnames{k} = ['FlyBowl_',expnames{k}];
    end
    
    if numel(gtdata(i).line_name{j}) <= 5,
      line_name = sprintf('GMR_%s_AE_01',gtdata(i).line_name{j});
      if ~any(strcmp({metadata.line_name},line_name)),
        line_name = sprintf('GMR_%s_AD_01',gtdata(i).line_name{j});
      end
    elseif strcmp(gtdata(i).line_name{j},'pBDPGAL4'),
      line_name = 'pBDPGAL4U';
    else
      line_name = gtdata(i).line_name{j};
    end
    idx = find(strcmp({metadata.line_name},line_name));
    
    for stati = 1:nstatsest(i),
      statfn = statfns_est.(behavior)(stati);    
      switch statfn.sex,
        case 'M',
          flys = 'male';
        case 'F',
          flys = 'female';
        otherwise
          flys = 'any';
      end
      switch behaviors{i},
        case 'Jump',
          frames = 'jump';
        case 'WingExtension'
          frames = 'wingextension';
      end
      
      fn = sprintf('%s_fly%s_frame%s',statfn.statfn,flys,frames);

      statpred.(behavior)(stati,j) = linestats.means.(fn)(strcmp(linestats.line_names,line_name));
      [ism,idx] = ismember(expnames,{metadata.experiment_name});
      statpred_perexp.(behavior){stati,j} = allstats.(fn)(idx);
      statpred_std.(behavior)(stati,j) = sqrt( nanmean( (statpred_perexp.(behavior){stati,j}-statpred.(behavior)(stati,j)).^2 ) );

    end
    
  end
  
end

%% plot fractime comparison

hfig = 10;
figure(hfig);
clf;
set(hfig,'Units','pixels','Position',[10, 10 1250 450]);
hax = createsubplots(1,2,[.075,.03;.2,.015]);

x0 = cumsum(skipsome*.5+1);
xl = x0'-.125;
xp = x0'+.125;

gtcolormean = [0,0,0];
gtcolor = [.5,.5,.5];
autocolormean = [.6,0,0];
autocolor = [1,.5,.5];
gtmarker = '.';
gtmarkermean = 'o';
automarker = '.';
automarkermean = 'o';

for i = 1:nbehaviors,
  
  axes(hax(i));
  hold(hax(i),'on');
  for j = 1:nlinesgt,
    hlabelexp = plot(xl(j),fractimelabel_perexp{i,j},gtmarker,'Color',gtcolor);
    hpredexp = plot(xp(j),fractimepred_perexp{i,j},automarker,'Color',autocolor);
  end
  plot([xl;xl;nan(1,nlinesgt)],[fractimelabel(i,:)+fractimelabel_std(i,:)
    fractimelabel(i,:)-fractimelabel_std(i,:)
    nan(1,nlinesgt)],'-','Color',gtcolormean);
  plot([xp;xp;nan(1,nlinesgt)],[fractimepred(i,:)+fractimepred_std(i,:)
    fractimepred(i,:)-fractimepred_std(i,:)
    nan(1,nlinesgt)],'-','Color',autocolormean);
  hlabel = plot(xl,fractimelabel(i,:),gtmarkermean,'Color',gtcolormean,'Marker',gtmarkermean,'MarkerFaceColor',gtcolormean);
  hpred = plot(xp,fractimepred(i,:),automarkermean,'Color','k','Marker',automarkermean,'MarkerFaceColor',autocolormean);

  if i == 1,
    ylabel(sprintf('Frac. time performing behavior'));
    hleg = legend([hlabel,hpred],{'Ground-truth','Prediction'});
    set(hleg,'Box','off');
  end
  
  axisalmosttight;
  ylim = get(hax(i),'YLim');
  ylim(1) = -ylim(2)*.01;
  set(hax(i),'YLim',ylim);
  
end

set(hax,'XLim',[0,x0(end)+1]);
set(hax,'XTick',x0,'XTickLabel',orderedtests(:,2)');

[orderedtesttypes,~,testtypeidx] = unique(orderedtests(:,3));
harrow = nan(numel(orderedtesttypes),nbehaviors);
htext = nan(numel(orderedtesttypes),nbehaviors);
for i = 1:numel(orderedtesttypes),
  idxcurr = find(testtypeidx == i);
  if numel(idxcurr) <= 1,
    continue;
  end
  x1 = x0(idxcurr(1))-.4;
  x2 = x0(idxcurr(end))+.4;
  for j = 1:nbehaviors,
    xfig = dsxy2figxy(hax(j),[x1,x2],[0,0]);
    pos = get(hax(j),'Position');
    yfig1 = pos(2)*.65 + [0,0];
    yfig2 = pos(2)*[.4,.65];

    harrow(i,j) = annotation('doublearrow',xfig,yfig1,...
      'Head1Style','rectangle','Head2Style','rectangle',...
      'Head1Length',1,'Head2Length',1);
    htext(i,j) = annotation('textbox',[xfig(1),yfig2(1),diff(xfig),diff(yfig2)],...
      'String',orderedtesttypes{i},'LineStyle','none','HorizontalAlignment','center');
      
  end
end

%SaveFigLotsOfWays(hfig,'GTFlyBowlPaper/GTFracTimeJumpWingExt20130925');

%% plot in a different way

hfig = 85;
figure(hfig);
clf;

i = 1;
behavior = behaviors{i};

% control behavior
controlmean = linestats.means.(statfn)(end);
controlstd = linestats.stds.(statfn)(end);

% normalize
zlabel = (fractimelabel(i,1:end-1)-controlmean)/controlstd;
zpred = (fractimepred(i,1:end-1)-controlmean)/controlstd;

minv = min(zlabel);
maxv = max(zlabel);

lim = [minv/1.05,maxv*1.05];

% equal line
plot(lim,lim,'c--');
hold on;

statfn = sprintf('fractime_flyany_frame%s',lower(behavior));

plot(max(minv,fractimelabel(i,1:end-1)),max(minv,fractimepred(i,1:end-1)),'wo','MarkerFaceColor',[.7,0,0]);
plot(repmat(fractimelabel(i,1:end-1),[2,1]),...
  [max(minv,fractimepred(i,1:end-1)-fractimepred_std(i,1:end-1))
  fractimepred(i,1:end-1)+fractimepred_std(i,1:end-1)],'-','Color',[.7,0,0]);
plot([max(minv,fractimelabel(i,1:end-1)-fractimelabel_std(i,1:end-1))
  fractimelabel(i,1:end-1)+fractimelabel_std(i,1:end-1)],...
  repmat(fractimepred(i,1:end-1),[2,1]),'-','Color',[.7,0,0]);
  
set(gca,'YScale','log','XScale','log');
xlabel('Frac. time labeled');
ylabel('Frac. time predicted');
box off;

axis([lim,lim]);

%% plot stat comparison

hfig = 12;
figure(hfig);
clf;
set(hfig,'Units','pixels','Position',[10, 10 1250 450]);
hax = createsubplots(1,sum(nstatsest),[.075,.05;.2,.015]);

order = 1:5;

x0 = cumsum(skipsome(order)*.5+1);
xl = x0'-.125;
xp = x0'+.125;

gtcolormean = [0,0,0];
gtcolor = [.5,.5,.5];
autocolormean = [.6,0,0];
autocolor = [1,.5,.5];
gtmarker = '.';
gtmarkermean = 'o';
automarker = '.';
automarkermean = 'o';

axi = 1;

for i = 1:nbehaviors,
  
  behavior = behaviors{i};
  
  isdata = ~isnan(statlabel.(behavior)) & ~isnan(statpred.(behavior));
  
  for stati = 1:nstatsest(i),
    
    statfn = statfns_est.(behavior)(stati);
    
    %isdata_order = isdata(stati,order);
  
    axes(hax(axi));
    hold(hax(axi),'on');
    
    for j = 1:numel(order),
      hlabelexp = plot(xl(j),statlabel_perexp.(behavior){stati,order(j)},gtmarker,'Color',gtcolor);
      hpredexp = plot(xp(j),statpred_perexp.(behavior){stati,order(j)},automarker,'Color',autocolor);
    end
    xcurr = [xl;xl;nan(1,numel(order))];
    ycurr = [statlabel.(behavior)(stati,order)+statlabel_std.(behavior)(stati,order)
      statlabel.(behavior)(stati,order)-statlabel_std.(behavior)(stati,order)
      nan(1,numel(order))];
    plot(xcurr,ycurr,'-','Color',gtcolormean);
    xcurr = [xp;xp;nan(1,numel(order))];
    ycurr = [statpred.(behavior)(stati,order)+statpred_std.(behavior)(stati,order)
      statpred.(behavior)(stati,order)-statpred_std.(behavior)(stati,order)
      nan(1,numel(order))];
    plot(xcurr,ycurr,'-','Color',autocolormean);
    hlabel = plot(xl,statlabel.(behavior)(stati,order),gtmarkermean,'Color',gtcolormean,'Marker',gtmarkermean,'MarkerFaceColor',gtcolormean);
    hpred = plot(xp,statpred.(behavior)(stati,order),automarkermean,'Color','k','Marker',automarkermean,'MarkerFaceColor',autocolormean);

    h = ylabel(sprintf('%s during %s',statfn.statfn,lower(behaviors{i})));
    set(h,'Interpreter','none');
%     if axi == 1,
%       hleg = legend([hlabel,hpred],{'Ground-truth','Prediction'});
%       set(hleg,'Box','off');
%     end
  
    if axi == 1,
      xlabel('Percentile of hits');
    end
    axisalmosttight;
      
    axi = axi + 1;
    
  end
end

set(hax,'XLim',[0,x0(end)+1]);
set(hax,'XTick',x0,'XTickLabel',orderedtests(order,2)');


SaveFigLotsOfWays(hfig,'GTFlyBowlPaper/GTStatsJumpWingExt20130925');
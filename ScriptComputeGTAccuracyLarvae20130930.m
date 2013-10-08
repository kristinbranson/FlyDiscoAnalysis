addpath /groups/branson/home/bransonk/behavioranalysis/code/Jdetect/Jdetect/perframe/
addpath /groups/branson/home/bransonk/behavioranalysis/code/Jdetect/Jdetect/misc;
addpath /groups/branson/home/bransonk/behavioranalysis/code/Jdetect/Jdetect/filehandling;
GTDataDir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/gtlarvaedata20130930';

%% parameters

behaviors = {'burrowing'};
jabfilestrs = {'allburrowsaug28gt.jab'};

jabfiles = cellfun(@(x) fullfile(GTDataDir,x),jabfilestrs,'UniformOutput',false);
nbehaviors = numel(behaviors);
intsize = 130;
rootexpdir = '/groups/branson/home/riveraalbam/Desktop/experiments';

%% compute fraction of time estimates using importance sampling

speciesnamesperexp = cell(1,nbehaviors);
newjabfiles = cell(1,nbehaviors);
frameweights = cell(1,nbehaviors);

for i = 1:nbehaviors,
  jabfile = jabfiles{i};
  
  % change experiment directory names
  jd = load(jabfile,'-mat');
  speciesnamesperexp{i} = cell(1,numel(jd.x.gtExpDirNames));
  for j = 1:numel(jd.x.gtExpDirNames),
    expdir = jd.x.gtExpDirNames{j};
    expdir = strrep(expdir,'\','/');
    [p,moviename] = fileparts(expdir);
    [~,speciesnamesperexp{i}{j}] = fileparts(p);
    newexpdir = fullfile(rootexpdir,speciesnamesperexp{i}{j},moviename);
    if ~exist(newexpdir,'dir'),
      error('Could not find %s (%s)',newexpdir,expdir);
    end
    jd.x.gtExpDirNames{j} = newexpdir;
  end
  
  [~,jabname] = myfileparts(jabfile);
  newjabfile = fullfile(GTDataDir,regexprep(jabname,'\.jab$','_kb.jab'));
  save(newjabfile,'-struct','jd');
  newjabfiles{i} = newjabfile;
  
  frameweights{i} = gtSuggestionWeights(jd,intsize);
end

speciesnames = unique([speciesnamesperexp{:}]);
nspecies = numel(speciesnames);

%% sanity check: if we re-weight the examples selected using these
% 1/weights, then we should get fractime simimlar to what is reported by
% the scores

fractimepred = nan(nbehaviors,nspecies);
fractimereal = nan(nbehaviors,nspecies);

for i = 1:nbehaviors,
  jabfile = newjabfiles{i};
  
  jabdata = load(jabfile,'-mat');
  scorefilename = jabdata.x.file.scorefilename;
  expdirs0 = jabdata.x.gtExpDirNames;
  labels0 = jabdata.x.gtLabels;
  
  for j = 1:nspecies,
     
    idxsp = find(strcmp(speciesnamesperexp{i},speciesnames{j}));

    expdirs = expdirs0(idxsp);
    
    weightpredpos = 0;
    weightpredneg = 0;
        
    labels = labels0(idxsp);
    
    nlabelstotal = 0;
    
    for k = 1:numel(expdirs),
      scorefile = fullfile(expdirs{k},scorefilename);
      
      expk = idxsp(k);
      
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
        idx = frameweights{i}.exp == expk & ...
          frameweights{i}.flies == fly & ...
          ismember(frameweights{i}.t,tpredpos);
        weightpredpos = weightpredpos + sum(1./frameweights{i}.wt(idx));
        
        tpredneg = find(sd.allScores.postprocessed{fly} == 0 & labelcurr ~= 0);
        idx = frameweights{i}.exp == expk & ...
          frameweights{i}.flies == fly & ...
          ismember(frameweights{i}.t,tpredneg);
        weightpredneg = weightpredneg + sum(1./frameweights{i}.wt(idx));
        
        nlabelstotal = nlabelstotal + nnz(labelcurr);
        
      end
      
      npredpos = 0;
      n = 0;
      for fly = 1:numel(sd.allScores.postprocessed),
        
        npredpos = npredpos + nnz(sd.allScores.postprocessed{fly} == 1);
        n = n + numel(sd.allScores.postprocessed{fly});
        
      end
      
    end
    
    fractimepred(i,j) = weightpredpos/(weightpredpos+weightpredneg);
    fractimereal(i,j) = npredpos/n;

    fprintf('%s: %d frames labeled, estimated fractime = %f, mean frac time = %f\n',...
      behaviors{i},nlabelstotal,fractimepred(i,j),fractimereal(i,j));
  end
end

%% estimate frac time statistics after resampling

% NOTE DID NOT MODIFY THIS PART AS THE PREVIOUS PART DID NOT WORK WELL
% ENOUGH

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
        ilabelpos = tlabelpos - td.trx(fly).firstframe + 1;
        idxsex = cellfun(@isempty,regexp(td.trx(fly).sex(ilabelpos),statfn.sex,'once'));
        tlabelpos(idxsex) = [];
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
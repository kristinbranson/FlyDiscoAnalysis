%% set up path

addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
addpath /groups/branson/bransonlab/projects/olympiad/SAGE/MATLABInterface/Trunk;

austindatadir = '/groups/branson/bransonlab/forAustin/MushroomBody_Data/for_bransonlab';
realdatadir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data';

%% parameters

minnframesperfly = 500;

plusstr = 'CTRL_CSMH_1500154_0030';
trpvk5str = 'UAS_dTrpA1_3_0062';
trpgstr = 'UAS_dTrpA1_2_0002';

effectors = {plusstr,trpvk5str,trpgstr};
PLUS = 1;
TRPVK5 = 2;
TRPG = 3;
neffectors = numel(effectors);

behaviorsplot = {
  'backup','scoresBackup.mat'
  'touch','scoresTouch.mat'
  'winggrooming','scoresWingGrooming.mat'
  'attemptedcopulation','scores_AttemptedCopulation.mat'
  'chase','scores_Chasev7.mat'
  'crabwalk','scores_Crabwalk2.mat'
  'jump','scores_Jump.mat'
  'righting','scores_Righting.mat'
  'stop','scores_Stops.mat'
  'walk','scores_Walk.mat'
  'wingextension','scores_WingExtension.mat'
  'pivotcenter','scores_pivot_center.mat'
  'pivottail','scores_pivot_tail.mat'
  'wingflick','scores_wingflick.mat'
  };

line_names = {
  'GMR_MB002B'
  'GMR_MB011B'
  'GMR_MB018B'
  'GMR_MB027B'
  'GMR_MB028B'
  'GMR_MB050B'
  'GMR_MB051B'
  'GMR_MB052B'
  'GMR_MB057B'
  'GMR_MB074C'
  'GMR_MB077B'
  'GMR_MB080C'
  'GMR_MB082C'
  'GMR_MB083C'
  'GMR_MB093C'
  'GMR_MB112C'
  'GMR_MB210B'
  'GMR_MB242A'
  'GMR_MB298B'
  'GMR_MB310C'
  'GMR_MB399B'
  'GMR_MB433B'
  'GMR_MB434B'
  'GMR_MB542B'
  'GMR_MB543B'
  'GMR_MB549C'
  'GMR_MB552B'
  };

nbehaviors = size(behaviorsplot,1);

%% locations of experiments to plot for each line

[expinfodirs,info] = mydir(fullfile(austindatadir,'*MB*'),'isdir',true);
idx = ismember({info.name},line_names);
expinfodirs = expinfodirs(idx);
info = info(idx);
missing_line_names = setdiff(line_names,{info.name});
if ~isempty(missing_line_names),
  error('Missing the following lines: %s',sprintf('%s ',missing_line_names{:}));
end
%line_names = {info.name};
nlines = numel(line_names);

expdirs = {};
metadata = [];
bdpidxperline = cell(1,nlines);

for linei = 1:nlines,
  
  expinfodir = expinfodirs{linei};
  line_name = line_names{linei};
  
  % read in lists of each type of experiment
  for effi = 1:neffectors,
  
    effector = effectors{effi};
    expfile = mydir(fullfile(expinfodir,sprintf('%s_%s*explist.txt',line_name,effector)));

    if isempty(expfile),
      error('Could not find exp list for %s/%s',line_name,effector);
    end
    expfile = expfile{1};
    expdirscurr = ReadAustinMBExpDirList(expfile);
    expdirs = [expdirs,expdirscurr]; %#ok<AGROW>
    for i = 1:numel(expdirscurr),
      
      metadatacurr = parseExpDir(expdirscurr{i},true);
      metadatacurr.effector = effector;
      metadata = structappend(metadata,metadatacurr);
      
    end

  end
  
  % also BDP
  bdpidxperline{linei} = [];
  for effi = 1:neffectors,
    effector = effectors{effi};
    expfile = mydir(fullfile(expinfodir,sprintf('pBDPGAL4U_%s*explist.txt',effector)));
    if isempty(expfile),
      error('Could not find exp list for %s/%s',line_name,effector);
    end
    expfile = expfile{1};
    expdirscurr = ReadAustinMBExpDirList(expfile);
    
    % which experiments are already on the list?
    [isold,idx] = ismember(expdirscurr,expdirs);
    
    bdpidxperline{linei} = [bdpidxperline{linei},idx(isold),numel(expdirs)+(1:nnz(~isold))];
    expdirs = [expdirs,expdirscurr(~isold)]; %#ok<AGROW>
    
    for i = find(~isold),
      
      metadatacurr = parseExpDir(expdirscurr{i},true);
      metadatacurr.effector = effector;
      metadata = structappend(metadata,metadatacurr);
      
    end
    
  end
  
end

metadata = AddSetMetadata(metadata);

nexps = numel(expdirs);

%% compute fraction of time for each experiment

expstats = struct;
expstats.Z = struct;
for behi = 1:nbehaviors,
  fn = behaviorsplot{behi,1};
  expstats.(fn) = nan(1,nexps);
  expstats.Z.(fn) = nan(1,nexps);
end

for expi = 1:nexps,
  
  fprintf('%d / %d\n',expi,nexps);

  for behi = 1:nbehaviors,
    fn = behaviorsplot{behi,1};
    scoresfile = fullfile(expdirs{expi},behaviorsplot{behi,2});
    if ~exist(scoresfile,'file'),
      warning('Scores file %s missing',scoresfile);
      break;
    end
    sd = load(scoresfile);
    nframescurr = cellfun(@numel,sd.allScores.postprocessed);
    x = [sd.allScores.postprocessed{nframescurr>=minnframesperfly}];
    x = x(~isnan(x));
    expstats.(fn)(expi) = nnz(x) / numel(x);
    expstats.Z.(fn)(expi) = numel(x);
  end
  
end

%% compute set stats

[set_names,firstidx,exp2setidx] = unique({metadata.set});
nsets = numel(set_names);

setstats = struct;
setstats.metadata = metadata(firstidx);
setstats.nexps = struct;
setstats.means = struct;
setstats.stds = struct;

for behi = 1:nbehaviors,
  fn = behaviorsplot{behi,1};
  setstats.means.(fn) = nan(1,nsets);
  setstats.stds.(fn) = nan(1,nsets);
  setstats.nexps.(fn) = nan(1,nsets);
end

for seti = 1:nsets,
  
  idx = find(exp2setidx==seti);
  if numel(idx) < 2,
    continue;
  end
  
  for behi = 1:nbehaviors,
    fn = behaviorsplot{behi,1};
    x = expstats.(fn)(idx);
    setstats.means.(fn)(seti) = nanmean(x);
    setstats.stds.(fn)(seti) = nanstd(x,1);
    setstats.nexps.(fn)(seti) = nnz(~isnan(x));
  end
  
end

bdpsetidxperline = cell(1,nlines);
for linei = 1:nlines,
  bdpsetidxperline{linei} = unique(exp2setidx(bdpidxperline{linei}));
end

%% compute genotype stats

genotypestats = struct;
genotypestats.nsets = struct;
genotypestats.means = struct;
genotypestats.stds = struct;
[genotypes,firstidx] = unique({metadata.genotype});
genotypestats.metadata = metadata(firstidx);
ngenotypes = numel(genotypes);

for behi = 1:nbehaviors,
  fn = behaviorsplot{behi,1};
  genotypestats.means.(fn) = nan(1,ngenotypes);
  genotypestats.stds.(fn) = nan(1,ngenotypes);
  genotypestats.nexps.(fn) = nan(1,ngenotypes);
end

for geni = 1:ngenotypes,
  
  genotype = genotypes{geni};
  idx = find(strcmp({setstats.metadata.genotype},genotype));
  
  for behi = 1:nbehaviors,
    fn = behaviorsplot{behi,1};
    x = setstats.means.(fn)(idx);
    genotypestats.nsets.(fn)(geni) = nnz(~isnan(x));
    genotypestats.means.(fn)(geni) = nanmean(x);
    genotypestats.stds.(fn)(geni) = nanstd(x);
  end

  idx = find(strcmp({metadata.genotype},genotype));
  
  for behi = 1:nbehaviors,
    fn = behaviorsplot{behi,1};
    x = expstats.(fn)(idx);
    genotypestats.nexps.(fn)(geni) = nnz(~isnan(x));
    genotypestats.expmeans.(fn)(geni) = nanmean(x);
    genotypestats.expstds.(fn)(geni) = nanstd(x);
  end

  
end

%% plot raw values 

maxoff = .25;
genotypeplotidx = find(strcmp({genotypestats.metadata.effector},trpgstr));
ngenotypesplot = numel(genotypeplotidx);
plusplotidx = nan(1,ngenotypesplot);
for i = 1:numel(genotypeplotidx),
  geni = genotypeplotidx(i);
  if strcmp(genotypestats.metadata(geni).line_name,'pBDPGAL4U'),
    continue;
  end
  g = sprintf('%s__%s',genotypestats.metadata(geni).line_name,plusstr);
  genj = find(strcmp({genotypestats.metadata.genotype},g));
  if isempty(genj),
    error('No gal4/plus for line %s',genotypestats.metadata(geni).line_name);
  end
  plusplotidx(i) = genj;
end

colors = [jet(ngenotypesplot-1);[0,0,0]];
shortlinenamesplot = {genotypestats.metadata(genotypeplotidx).line_name};
shortlinenamesplot = regexprep(shortlinenamesplot,'^GMR_','');


for behi = 1:nbehaviors,
  hfig = behi + 100;
  figure(hfig);
  clf;
  set(hfig,'Units','pixels','Position',[10,10,1800,800]);
  hold on;
  hax = gca;
  set(hax,'Position',[ 0.0311    0.1100    0.9439    0.8638])
  
  
  fn = behaviorsplot{behi,1};
  
  [~,genorder] = sort(genotypestats.means.(fn)(genotypeplotidx));
  
  for i = 1:ngenotypesplot,
    
    x0 = 2*i-1;
    
    geni = genotypeplotidx(genorder(i));
    genotype = genotypes{geni};
    
    % gather the per-experiment and per-set data
    
    setidx = find(strcmp({setstats.metadata.genotype},genotype));
    nsetscurr = numel(setidx);

    expy = cell(1,nsetscurr);
    sety = nan(1,nsetscurr);
    for j = 1:nsetscurr,
      seti = setidx(j);
      setname = setstats.metadata(seti).set;
      expidx = find(strcmp({metadata.set},setname));
      
      expy{j} = sort(expstats.(fn)(expidx));
      expy{j}(isnan(expy{j})) = [];
      sety(j) = setstats.means.(fn)(seti);
    end
    
    badidx = isnan(sety);
    sety(badidx) = [];
    expy(badidx) = [];
      
    nsetscurr = numel(sety);
    
    % how much should we offset it?
    if nsetscurr == 0,
      continue;
    elseif nsetscurr == 1,
      x = x0;
    else
      x = linspace(-maxoff,maxoff,nsetscurr);
      x = x - mean(x) + x0;
    end
    
    % plot experiment data
    color = colors(genorder(i),:)*.8;
    for j = 1:nsetscurr,
      hexp(i,j) = plot(repmat(x(j),size(expy{j})),expy{j},'.-','Color',(color+1)/2);
    end
    
    % plot set data
    hset(i) = plot(x,sety,'*','Color',(color+1)/2);
    
    % plot genotype stderr data
    s = genotypestats.stds.(fn)(geni) / sqrt(nsetscurr);
    hgenstd(i) = plot([x0,x0],genotypestats.means.(fn)(geni)+s*[-1,1],'-','Color',color*.8);

    % plot genotype mean data
    hgen(i) = plot(x0,genotypestats.means.(fn)(geni),'o','Color','k','MarkerFaceColor',color*.8);
    
  end
  
  % plot gal4/plus
  
  for i = 1:ngenotypesplot,
    
    x0 = 2*i-.5;
    
    geni = plusplotidx(genorder(i));
    if isnan(geni),
      continue;
    end
    genotype = genotypes{geni};
    
    % gather the per-experiment and per-set data
    
    setidx = find(strcmp({setstats.metadata.genotype},genotype));
    nsetscurr = numel(setidx);

    expy = cell(1,nsetscurr);
    sety = nan(1,nsetscurr);
    for j = 1:nsetscurr,
      seti = setidx(j);
      setname = setstats.metadata(seti).set;
      expidx = find(strcmp({metadata.set},setname));
      
      expy{j} = sort(expstats.(fn)(expidx));
      expy{j}(isnan(expy{j})) = [];
      sety(j) = setstats.means.(fn)(seti);
    end
    
    badidx = isnan(sety);
    sety(badidx) = [];
    expy(badidx) = [];
      
    nsetscurr = numel(sety);
    
    % how much should we offset it?
    if nsetscurr == 0,
      continue;
    elseif nsetscurr == 1,
      x = x0;
    else
      x = linspace(-maxoff,maxoff,nsetscurr);
      x = x - mean(x) + x0;
    end
    
    % plot experiment data
    color = [0,0,0];
    for j = 1:nsetscurr,
      hexpplus(i,j) = plot(repmat(x(j),size(expy{j})),expy{j},'.-','Color',(color+1)/2);
    end
    
    % plot set data
    hsetplus(i) = plot(x,sety,'x','Color',(color+1)/2);
    
    % plot genotype stderr data
    s = genotypestats.stds.(fn)(geni) / sqrt(nsetscurr);
    hgenstdplus(i) = plot([x0,x0],genotypestats.means.(fn)(geni)+s*[-1,1],'-','Color',color);

    % plot genotype mean data
    hgenplus(i) = plot(x0,genotypestats.means.(fn)(geni),'s','Color','k','MarkerFaceColor',color);
    
  end
  
  
  axisalmosttight;
  box off;
  ylim = get(gca,'YLim');
  ylim(1) = min(ylim(1),0);
  set(gca,'XLim',[0,2*ngenotypesplot+1.5],'YLim',ylim);
  set(gca,'XTick',1.25:2:2*ngenotypesplot+.25,'XTickLabel',shortlinenamesplot(genorder))
  htick = rotateticklabel(gca);
  ylabel(sprintf('Frac. time %s',fn));
  
  drawnow;
  %SaveFigLotsOfWays(hfig,sprintf('MBResults20131101/OutputNeurons_FracTime%s',fn));
  
end

%% collect data for plotting normalized values

maxoff = .25;
mbtrpgplotidx = find(strcmp({genotypestats.metadata.effector},trpgstr));% & ~strcmp({genotypestats.metadata.line_name},'pBDPGAL4U'));
ngenotypesplot = numel(mbtrpgplotidx);
mbtrpvk5plotidx = nan(1,ngenotypesplot);
mbplusplotidx = nan(1,ngenotypesplot);

for i = 1:numel(mbtrpgplotidx),

  mbtrpi = mbtrpgplotidx(i);
%   if strcmp(genotypestats.metadata(mbtrpi).line_name,'pBDPGAL4U'),
%     continue;
%   end

  % gal4 / dtrpa1vk5
  g = sprintf('%s__%s',genotypestats.metadata(mbtrpi).line_name,trpvk5str);
  genj = find(strcmp({genotypestats.metadata.genotype},g));
  if isempty(genj),
    error('No gal4/dtrpa1vk5 for line %s',genotypestats.metadata(mbtrpi).line_name);
  end
  mbtrpvk5plotidx(i) = genj;

  % gal4 / +
  g = sprintf('%s__%s',genotypestats.metadata(mbtrpi).line_name,plusstr);
  genj = find(strcmp({genotypestats.metadata.genotype},g));
  if isempty(genj),
    error('No gal4/plus for line %s',genotypestats.metadata(mbtrpi).line_name);
  end
  mbplusplotidx(i) = genj;
  
end

% bdp / dtrpa1g
g = sprintf('%s__%s','pBDPGAL4U',trpgstr);
genj = find(strcmp({genotypestats.metadata.genotype},g));
if isempty(genj),
  error('No BDP/dtrpa1g found');
end
bdptrpgi = genj;
bdptrpg_genotype = genotypes{bdptrpgi};

% bdp / dtrpa1vk5
g = sprintf('%s__%s','pBDPGAL4U',trpvk5str);
genj = find(strcmp({genotypestats.metadata.genotype},g));
if isempty(genj),
  error('No BDP/dtrpa1vk5 found');
end
bdptrpvk5i = genj;
bdptrpvk5_genotype = genotypes{bdptrpvk5i};

%% plot normalized values

colors = [jet(ngenotypesplot-1);0,0,0];
markers = {'o','s','d','v'};
shortlinenamesplot = {genotypestats.metadata(mbtrpgplotidx).line_name};
shortlinenamesplot = regexprep(shortlinenamesplot,'^GMR_','');

for behi = 1:nbehaviors,

  hfig = behi + 200;
  figure(hfig);
  clf;
  set(hfig,'Units','pixels','Position',[10,10,1800,800]);
  hold on;
  hax = gca;
  set(hax,'Position',[ 0.0311    0.1100    0.9439    0.8638])
  
  
  fn = behaviorsplot{behi,1};
  
  [~,genorder] = sort(genotypestats.means.(fn)(mbtrpgplotidx));
  
  for i = 1:ngenotypesplot,
        
    mbtrpgi = mbtrpgplotidx(genorder(i));
    mbtrpvk5i = mbtrpvk5plotidx(genorder(i));
    mbplusi = mbplusplotidx(genorder(i));

    mbtrpg_genotype = genotypes{mbtrpgi};
    mbtrpvk5_genotype = genotypes{mbtrpvk5i};
    mbplus_genotype = genotypes{mbplusi};

    % mb/trpg
    
    % gather the per-experiment and per-set data    
    [sety0,expy0] = CollectSetDataGivenGenotype(fn,mbtrpg_genotype,setstats,expstats,metadata);
    y0 = genotypestats.means.(fn)(mbtrpgi);
    s1 = genotypestats.stds.(fn)(mbtrpgi) / sqrt(nsetscurr);
    nsetscurr = numel(sety0);
    
    % mb/trpg - mb/plus    
    
    typei = 1;
    x0 = 4*(i-1)+typei;
    
    controly = genotypestats.means.(fn)(mbplusi);
    s2 = genotypestats.stds.(fn)(mbplusi) / sqrt(genotypestats.nsets.(fn)(mbplusi));
    s = s1 + s2;
    
    % subtract off control
    y = y0 - controly;
    sety = sety0 - controly;
    expy = cellfun(@(x) x-controly,expy0,'UniformOutput',false);
    
    % how much should we offset it?
    if nsetscurr == 0,
      continue;
    elseif nsetscurr == 1,
      x = x0;
    else
      x = linspace(-maxoff,maxoff,nsetscurr);
      x = x - mean(x) + x0;
    end
    
    % plot experiment data
    color = colors(genorder(i),:)*.8;
    for j = 1:nsetscurr,
      hexp(typei,i,j) = plot(repmat(x(j),size(expy{j})),expy{j},'.-','Color',(color+1)/2);
    end
    
    % plot set data
    hset(typei,i) = plot(x,sety,'*','Color',(color+1)/2);
    
    % plot genotype stderr data
    hgenstd(typei,i) = plot([x0,x0],y+s*[-1,1],'-','Color',color*.8);

    % plot genotype mean data
    hgen(typei,i) = plot(x0,y,markers{typei},'Color','k','MarkerFaceColor',color*.8);
    
    
    
    
    % mb/trpg - bdp/trpg
    
    typei = 2;
    
    x0 = 4*(i-1)+typei;
    
    controly = genotypestats.means.(fn)(bdptrpgi);
    s2 = genotypestats.stds.(fn)(bdptrpgi) / sqrt(genotypestats.nsets.(fn)(bdptrpgi));
    s = s1 + s2;
    
    % subtract off control
    y = y0 - controly;
    sety = sety0 - controly;
    expy = cellfun(@(x) x-controly,expy0,'UniformOutput',false);
    
    % how much should we offset it?
    if nsetscurr == 0,
      continue;
    elseif nsetscurr == 1,
      x = x0;
    else
      x = linspace(-maxoff,maxoff,nsetscurr);
      x = x - mean(x) + x0;
    end
    
    % plot experiment data
    color = colors(genorder(i),:)*.8;
    for j = 1:nsetscurr,
      hexp(typei,i,j) = plot(repmat(x(j),size(expy{j})),expy{j},'.-','Color',(color+1)/2);
    end
    
    % plot set data
    hset(typei,i) = plot(x,sety,'*','Color',(color+1)/2);
    
    % plot genotype stderr data
    hgenstd(typei,i) = plot([x0,x0],y+s*[-1,1],'-','Color',color*.8);

    % plot genotype mean data
    hgen(typei,i) = plot(x0,y,markers{typei},'Color','k','MarkerFaceColor',color*.8);

    
    
    % mb/trpvk5
    
    % gather the per-experiment and per-set data    
    [sety0,expy0] = CollectSetDataGivenGenotype(fn,mbtrpvk5_genotype,setstats,expstats,metadata);
    y0 = genotypestats.means.(fn)(mbtrpvk5i);
    s1 = genotypestats.stds.(fn)(mbtrpvk5i) / sqrt(nsetscurr);
    nsetscurr = numel(sety0);
    
    % mb/trpvk5 - mb/plus    
    
    typei = 3;
    x0 = 4*(i-1)+typei;
    
    controly = genotypestats.means.(fn)(mbplusi);
    s2 = genotypestats.stds.(fn)(mbplusi) / sqrt(genotypestats.nsets.(fn)(mbplusi));
    s = s1 + s2;
    
    % subtract off control
    y = y0 - controly;
    sety = sety0 - controly;
    expy = cellfun(@(x) x-controly,expy0,'UniformOutput',false);
    
    % how much should we offset it?
    if nsetscurr == 0,
      continue;
    elseif nsetscurr == 1,
      x = x0;
    else
      x = linspace(-maxoff,maxoff,nsetscurr);
      x = x - mean(x) + x0;
    end
    
    % plot experiment data
    color = colors(genorder(i),:)*.8;
    for j = 1:nsetscurr,
      hexp(typei,i,j) = plot(repmat(x(j),size(expy{j})),expy{j},'.-','Color',(color+1)/2);
    end
    
    % plot set data
    hset(typei,i) = plot(x,sety,'*','Color',(color+1)/2);
    
    % plot genotype stderr data
    hgenstd(typei,i) = plot([x0,x0],y+s*[-1,1],'-','Color',color*.8);

    % plot genotype mean data
    hgen(typei,i) = plot(x0,y,markers{typei},'Color','k','MarkerFaceColor',color*.8);
    
    
    
    
    % mb/trpvk5 - bdp/trpvk5
    
    typei = 4;
    
    x0 = 4*(i-1)+typei;
    
    controly = genotypestats.means.(fn)(bdptrpvk5i);
    s2 = genotypestats.stds.(fn)(bdptrpvk5i) / sqrt(genotypestats.nsets.(fn)(bdptrpvk5i));
    s = s1 + s2;
    
    % subtract off control
    y = y0 - controly;
    sety = sety0 - controly;
    expy = cellfun(@(x) x-controly,expy0,'UniformOutput',false);
    
    % how much should we offset it?
    if nsetscurr == 0,
      continue;
    elseif nsetscurr == 1,
      x = x0;
    else
      x = linspace(-maxoff,maxoff,nsetscurr);
      x = x - mean(x) + x0;
    end
    
    % plot experiment data
    color = colors(genorder(i),:)*.8;
    for j = 1:nsetscurr,
      hexp(typei,i,j) = plot(repmat(x(j),size(expy{j})),expy{j},'.-','Color',(color+1)/2);
    end
    
    % plot set data
    hset(typei,i) = plot(x,sety,'*','Color',(color+1)/2);
    
    % plot genotype stderr data
    hgenstd(typei,i) = plot([x0,x0],y+s*[-1,1],'-','Color',color*.8);

    % plot genotype mean data
    hgen(typei,i) = plot(x0,y,markers{typei},'Color','k','MarkerFaceColor',color*.8);

    
    
    
  end
  

  
  
  axisalmosttight;
  box off;
  ylim = get(gca,'YLim');
  ylim(1) = min(ylim(1),0);
  set(gca,'XLim',[0,4*ngenotypesplot+1.5],'YLim',ylim);
  set(gca,'XTick',2.25:4:4*ngenotypesplot+1.25,'XTickLabel',shortlinenamesplot(genorder))
  htick = rotateticklabel(gca);
  ylabel(sprintf('Difference in frac. time %s',fn));
  
  i = find(genorder==nlines);
  legend(hgen(:,i),{'GAL4/dTrpA1 Garrity - GAL4/+','GAL4/dTrpA1 Garrity - BDP/TrpA1 Garrity',...
    'GAL4/dTrpA1 VK5 - GAL4/+','GAL4/dTrpA1 VK5 - BDP/TrpA1 VK5'},...
    'Location','NorthWest');
  
  drawnow;
  SaveFigLotsOfWays(hfig,sprintf('MBResults20131101/OutputNeurons_FracTimeDiff%s',fn));

  
end
  
%% plot smallest difference

maxoff = .1;
colors = [jet(ngenotypesplot-1);0,0,0];
markers = {'x','o','s'};
shortlinenamesplot = {genotypestats.metadata(mbtrpgplotidx).line_name};
shortlinenamesplot = regexprep(shortlinenamesplot,'^GMR_','');
MBTRP_MINUS_MBPLUS = 1;
MBTRP_MINUS_BDPTRP = 2;
ncompares = 2;


for behi = 1:nbehaviors,

  hfig = behi + 200;
  figure(hfig);
  clf;
  set(hfig,'Units','pixels','Position',[10,10,1800,800]);
  hold on;
  hax = gca;
  set(hax,'Position',[ 0.0311    0.1100    0.9439    0.8638])
  
  
  fn = behaviorsplot{behi,1};
  
  [~,genorder] = sort(genotypestats.means.(fn)(mbtrpgplotidx));
  allminis = nan(1,ngenotypesplot);
  
  for i = 1:ngenotypesplot,
    
    ss = nan(1,ncompares);
    ys = nan(1,ncompares);
    setys = cell(1,ncompares);
    expys = cell(1,ncompares);
        
    mbtrpgi = mbtrpgplotidx(genorder(i));
    mbplusi = mbplusplotidx(genorder(i));

    mbtrpg_genotype = genotypes{mbtrpgi};
    mbplus_genotype = genotypes{mbplusi};

    % mb/trpg    
    % gather the per-experiment and per-set data    
    [sety0,expy0] = CollectSetDataGivenGenotype(fn,mbtrpg_genotype,setstats,expstats,metadata);
    y0 = genotypestats.means.(fn)(mbtrpgi);
    s1 = genotypestats.stds.(fn)(mbtrpgi) / sqrt(nsetscurr);
    nsetscurr = numel(sety0);
    if nsetscurr == 0,
      continue;
    end
    
    % mb/trpg - mb/plus    
    controly = genotypestats.means.(fn)(mbplusi);
    s2 = genotypestats.stds.(fn)(mbplusi) / sqrt(genotypestats.nsets.(fn)(mbplusi));
    s = s1 + s2;
    
    % subtract off control
    y = y0 - controly;
    sety = sety0 - controly;
    expy = cellfun(@(x) x-controly,expy0,'UniformOutput',false);
    
    ss(MBTRP_MINUS_MBPLUS) = s;
    ys(MBTRP_MINUS_MBPLUS) = y;
    setys{MBTRP_MINUS_MBPLUS} = sety;
    expys{MBTRP_MINUS_MBPLUS} = expy;

    % mb/trpg - bdp/trpg
    
    controly = genotypestats.means.(fn)(bdptrpgi);
    s2 = genotypestats.stds.(fn)(bdptrpgi) / sqrt(genotypestats.nsets.(fn)(bdptrpgi));
    s = s1 + s2;
    
    % subtract off control
    y = y0 - controly;
    sety = sety0 - controly;
    expy = cellfun(@(x) x-controly,expy0,'UniformOutput',false);
    
    ss(MBTRP_MINUS_BDPTRP) = s;
    ys(MBTRP_MINUS_BDPTRP) = y;
    setys{MBTRP_MINUS_BDPTRP} = sety;
    expys{MBTRP_MINUS_BDPTRP} = expy;

    % choose the minimum
    isambiguous = any(sign(ys)>0) && any(sign(ys)<0);
    if isambiguous,
      mini = 0;
      y = 0;
      s = nan;
      setys = {};
      expys = {};
    else
      [~,mini] = min(abs(ys./ss));
      y = ys(mini);
      s = ss(mini);
      sety = setys{mini};
      expy = expys{mini};
    end
    allminis(i) = mini;
    
    % how much should we offset it?
    x0 = i;    
    if nsetscurr == 1 || isambiguous,
      x = x0;
    else
      x = linspace(-maxoff,maxoff,nsetscurr);
      x = x - mean(x) + x0;
    end
    
    color = colors(genorder(i),:)*.8;

    if ~isambiguous,
  
      % plot experiment data
      for j = 1:nsetscurr,
        hexp(i,j) = plot(repmat(x(j),size(expy{j})),expy{j},'.-','Color',(color+1)/2);
      end
    
      % plot set data
      hset(i) = plot(x,sety,'*','Color',(color+1)/2);
    
      % plot genotype stderr data
      hgenstd(i) = plot([x0,x0],y+s*[-1,1],'-','Color',color*.8);

    end

    % plot genotype mean data
    hgen(i) = plot(x0,y,markers{mini+1},'Color','k','MarkerFaceColor',color*.8);
    
  end
  
  
  axisalmosttight;
  box off;
  ylim = get(gca,'YLim');
  ylim(1) = min(ylim(1),0);
  set(gca,'XLim',[0,ngenotypesplot+1],'YLim',ylim);
  set(gca,'XTick',1:ngenotypesplot,'XTickLabel',shortlinenamesplot(genorder))
  htick = rotateticklabel(gca);
  ylabel(sprintf('Conservative difference in frac. time %s',fn));
  
  hleg = [];
  sleg = {};
  i = find(allminis==MBTRP_MINUS_MBPLUS,1,'last');
  if ~isempty(i),
    hleg(end+1) = hgen(i);
    sleg{end+1} = 'GAL4/dTrpA1 Garrity - GAL4/+';
  end
  i = find(allminis==MBTRP_MINUS_BDPTRP,1,'last');
  if ~isempty(i),
    hleg(end+1) = hgen(i);
    sleg{end+1} = 'GAL4/dTrpA1 Garrity - BDP/TrpA1 Garrity';
  end
  i = find(allminis==0,1,'last');
  if ~isempty(i),
    hleg(end+1) = hgen(i);
    sleg{end+1} = 'Differences had opposite signs';
  end
  
  legend(hleg,sleg,'Location','NorthWest');

  
  drawnow;
  %SaveFigLotsOfWays(hfig,sprintf('MBResults20131101/OutputNeurons_FracTimeDiff%s',fn));

  
end

%% plot normalized values

colors = jet(ngenotypesplot-1);
colors = colors(randperm(ngenotypesplot-1),:);
colors = [colors;0,0,0];

markers = {'o','s','d','v'};
shortlinenamesplot = {genotypestats.metadata(mbtrpgplotidx).line_name};
shortlinenamesplot = regexprep(shortlinenamesplot,'^GMR_','');
ntypes = 2;
clear hgen;
maxoff = .2;
markersize = 7;
linewidth = 2;

for behi = 1:nbehaviors,

  hfig = behi + 200;
  figure(hfig);
  clf;
  set(hfig,'Units','pixels','Position',[10,10,630,400]);
  hold on;
  hax = gca;
  set(hax,'Position',[0.0904761904761905 0.205 0.884523809523809 0.768799999999999])
  
  
  fn = behaviorsplot{behi,1};
  
  %[~,genorder] = sort(genotypestats.means.(fn)(mbtrpgplotidx));
  genorder = 1:ngenotypesplot;
  
  plot([0,ngenotypesplot+1],[0,0],'--','Color',[.5,.5,.5]);
  
  [~,expcontrol] = CollectSetDataGivenGenotype(fn,bdptrpg_genotype,setstats,expstats,metadata);
  controly_bdptrpg = mean([expcontrol{:}]);
  s2_bdptrpg = std([expcontrol{:}],0) / sqrt(numel([expcontrol{:}]));
  
  
  for i = 1:ngenotypesplot,
    
    if ntypes == 1,
      xoffs = 0;
    else
      xoffs = linspace(-maxoff,maxoff,ntypes);
    end
        
    mbtrpgi = mbtrpgplotidx(genorder(i));
    mbtrpvk5i = mbtrpvk5plotidx(genorder(i));
    mbplusi = mbplusplotidx(genorder(i));

    mbtrpg_genotype = genotypes{mbtrpgi};
    mbtrpvk5_genotype = genotypes{mbtrpvk5i};
    mbplus_genotype = genotypes{mbplusi};

    % mb/trpg
    
    % gather the per-experiment and per-set data    
    [sety0,expy0] = CollectSetDataGivenGenotype(fn,mbtrpg_genotype,setstats,expstats,metadata);
    nsetscurr = numel(sety0);
    y0 = mean([expy0{:}]);
    s1 = std([expy0{:}],0) / sqrt(numel([expy0{:}]));

%     y0 = genotypestats.means.(fn)(mbtrpgi);
%     if genotypestats.nsets.(fn)(mbtrpgi) > 1,
%       s1 = genotypestats.stds.(fn)(mbtrpgi) / sqrt(nsetscurr);
%     else
%       s1 = genotypestats.expstds.(fn)(mbtrpgi);
%     end
    
    % mb/trpg - mb/plus    
    
    typei = 1;
    x0 = i + xoffs(typei);
    
    [~,expcontrol] = CollectSetDataGivenGenotype(fn,mbplus_genotype,setstats,expstats,metadata);
    controly = mean([expcontrol{:}]);
    s2 = std([expcontrol{:}],0) / sqrt(numel([expcontrol{:}]));
    
%     controly = genotypestats.means.(fn)(mbplusi);
%     if genotypestats.nsets.(fn)(mbplusi) > 1,
%       s2 = genotypestats.stds.(fn)(mbplusi) / sqrt(nsetscurr);
%     else
%       s2 = genotypestats.expstds.(fn)(mbplusi);
%     end

    if strcmp(mbplus_genotype,mbtrpg_genotype),
      s = 0;
    else
      s = sqrt(s1^2 + s2^2);
    end
    
    % subtract off control
    y = y0 - controly;
%     sety = sety0 - controly;
%     expy = cellfun(@(x) x-controly,expy0,'UniformOutput',false);
    
%     % how much should we offset it?
%     if nsetscurr == 0,
%       continue;
%     elseif nsetscurr == 1,
%       x = x0;
%     else
%       x = linspace(-maxoff,maxoff,nsetscurr);
%       x = x - mean(x) + x0;
%     end
%     
    % plot experiment data
    color = colors(genorder(i),:)*.8;
%     for j = 1:nsetscurr,
%       hexp(typei,i,j) = plot(repmat(x(j),size(expy{j})),expy{j},'.-','Color',(color+1)/2);
%     end
    
%     % plot set data
%     hset(typei,i) = plot(x,sety,'*','Color',(color+1)/2);
    
    % plot genotype stderr data
    hgenstd(typei,i) = plot([x0,x0],y+s*[-1,1],'-','Color',color*.8,'LineWidth',linewidth);

    % plot genotype mean data
    hgen(typei,i) = plot(x0,y,markers{typei},'Color','k','MarkerFaceColor',color*.8,'MarkerSize',markersize);
    
    
    
    
    % mb/trpg - bdp/trpg
    
    typei = 2;
    x0 = i + xoffs(typei);
    
    controly = controly_bdptrpg;
    s2 = s2_bdptrpg;
    s = sqrt(s1^2 + s2^2);
    
    % subtract off control
    y = y0 - controly;
    sety = sety0 - controly;
    expy = cellfun(@(x) x-controly,expy0,'UniformOutput',false);
    
    % how much should we offset it?
    if nsetscurr == 0,
      continue;
    elseif nsetscurr == 1,
      x = x0;
    else
      x = linspace(-maxoff,maxoff,nsetscurr);
      x = x - mean(x) + x0;
    end
    
    % plot experiment data
    color = colors(genorder(i),:)*.8;
%     for j = 1:nsetscurr,
%       hexp(typei,i,j) = plot(repmat(x(j),size(expy{j})),expy{j},'.-','Color',(color+1)/2);
%     end
    
    % plot set data
%     hset(typei,i) = plot(x,sety,'*','Color',(color+1)/2);
    
    % plot genotype stderr data
    hgenstd(typei,i) = plot([x0,x0],y+s*[-1,1],'-','Color',color*.8,'LineWidth',linewidth);

    % plot genotype mean data
    hgen(typei,i) = plot(x0,y,markers{typei},'Color','k','MarkerFaceColor',color*.8,'MarkerSize',markersize);

    
  end
  

  
  
  axisalmosttight;
  box off;
  ylim = get(gca,'YLim');
  ylim = max(abs(ylim))*[-1,1];
  %ylim(1) = min(ylim(1),0);
  set(gca,'XLim',[0,ngenotypesplot+1],'YLim',ylim);
  set(gca,'XTick',1:ngenotypesplot,'XTickLabel',shortlinenamesplot(genorder));
  htick = rotateticklabel(gca);
  for i = 1:numel(htick),
    set(htick(i),'Color',.8^2*colors(i,:));
  end
  ylabel(sprintf('Difference in frac. time %s',fn));
  
  i = find(strcmp(shortlinenamesplot(genorder),'pBDPGAL4U'));
  legend(hgen(:,i),{'GAL4/dTrpA1 Garrity - GAL4/+','GAL4/dTrpA1 Garrity - BDP/TrpA1 Garrity'},'Location','SouthEast');
  
  drawnow;
  SaveFigLotsOfWays(hfig,sprintf('MBResults20140211/OutputNeurons_FracTimeDiff%s',fn));

  
end

%% compute p-values using the rank-sum test

pvalue = nan(ngenotypesplot,nbehaviors);
p_mbplus = nan(ngenotypesplot,nbehaviors);
p_bdptrpg = nan(ngenotypesplot,nbehaviors);
d_mbplus = nan(ngenotypesplot,nbehaviors);
d_bdptrpg = nan(ngenotypesplot,nbehaviors);
nexps_mbtrpg = nan(ngenotypesplot,1);
nexps_mbplus = nan(ngenotypesplot,1);
nexps_bdptrpg = nan;

for behi = 1:nbehaviors,  
  
  fn = behaviorsplot{behi,1};
  fprintf('Behavior %s %d/%d\n',fn,behi,nbehaviors);

  % gather per-experiment data for BDP
  [~,expx_bdptrpg] = CollectSetDataGivenGenotype(fn,bdptrpg_genotype,setstats,expstats,metadata);
  mu_bdptrpg = median([expx_bdptrpg{:}]);
  if behi == 1,
    nexps_bdptrpg = numel([expx_bdptrpg{:}]);
  end
  
  for i = 1:ngenotypesplot,

    fprintf('Genotype %s %d/%d\n',mbtrpg_genotype,i,ngenotypesplot);
    
    mbtrpgi = mbtrpgplotidx(genorder(i));
    mbplusi = mbplusplotidx(genorder(i));

    mbtrpg_genotype = genotypes{mbtrpgi};
    mbplus_genotype = genotypes{mbplusi};

    % gather the per-experiment data

    % mb/trpg
    [~,expx_mbtrpg] = CollectSetDataGivenGenotype(fn,mbtrpg_genotype,setstats,expstats,metadata);
    % mb/plus
    [~,expx_mbplus] = CollectSetDataGivenGenotype(fn,mbplus_genotype,setstats,expstats,metadata);

    if behi == 1,
      nexps_mbtrpg(i) = numel([expx_mbtrpg{:}]);
      nexps_mbplus(i) = numel([expx_mbplus{:}]);
    end

    
    % gather means
    mu_mbtrpg = median([expx_mbtrpg{:}]);
    mu_mbplus = median([expx_mbplus{:}]);
    
    % see if signs of differences are the same
    d_mbplus(i,behi) = mu_mbtrpg - mu_mbplus;
    d_bdptrpg(i,behi) = mu_mbtrpg - mu_bdptrpg;
    
    p_mbplus(i,behi) = ranksum([expx_mbtrpg{:}],[expx_mbplus{:}],'method','exact');
    if numel([expx_mbtrpg{:}]) > 8 && numel([expx_bdptrpg{:}]) > 8,
      method = 'approximate';
    else
      method = 'exact';
    end
    p_bdptrpg(i,behi) = ranksum([expx_mbtrpg{:}],[expx_bdptrpg{:}],'method',method);
    
    % compute pvalue
    if sign(d_mbplus(i,behi))*sign(d_bdptrpg(i,behi)) < 0,
      pvalue(i,behi) = 1;
    else
      pvalue(i,behi) = max(p_mbplus(i,behi),p_bdptrpg(i,behi));
    end
    
  end
  
end

%% compute false discovery rate corrections

fdr = .05;
idxtest = sign(d_mbplus).*sign(d_bdptrpg) >= 0;
issig_mbplus = false(ngenotypesplot,nbehaviors);
adjp_mbplus = ones(ngenotypesplot,nbehaviors);
[issig_mbplus(idxtest),critp_mbplus,adjp_mbplus(idxtest)] = fdr_bh(p_mbplus(idxtest),fdr,'pdep');

issig_bdptrpg = false(ngenotypesplot,nbehaviors);
adjp_bdptrpg = ones(ngenotypesplot,nbehaviors);
[issig_bdptrpg(idxtest),critp_bdptrpg,adjp_bdptrpg(idxtest)] = fdr_bh(p_bdptrpg(idxtest),fdr,'pdep');

adjp_max = max(adjp_bdptrpg,adjp_mbplus);
dir = sign(d_bdptrpg);
dir(~idxtest) = 0;

figure(3442);
clf;
tmp = log10(adjp_max).*dir;
imagesc(tmp);
ntmp = ceil(abs(-2-log10(.05))/4*256);
cmtmp = linspace(0,.95,ntmp)';
cm = [ones(ntmp,1),repmat(cmtmp,[1,2])
  ones(256-2*ntmp,3)
  repmat(flipud(cmtmp),[1,2]),ones(ntmp,1)];
colormap(cm);
hcb = colorbar;
set(gca,'CLim',[-2,2]);
set(gca,'XTick',1:nbehaviors,'XTickLabel',behaviorsplot(:,1));
set(gca,'YTick',1:ngenotypesplot,'YTickLabel',shortlinenamesplot)
rotateticklabel(gca);
box off;
set(hcb,'YTick',[-2,log10(.05),-1,1,-log10(.05),2],'YTickLabel',{.01,.05,.1,.1,.05,.01});

%% output latex

shortbehaviorsplot = {'Backup','Touch','Wing groom','Att. cop.','Chase','Crab-walk','Jump','Righting','Stop','Walk','Wing ext.','Pivot-center','Pivot-tail','Wing-flick'};

fid = fopen('MBAdjPValueTable20140211.tex','w');

fprintf(fid,'\\documentclass{article}\n');
fprintf(fid,'\\usepackage[table]{xcolor}\n');
fprintf(fid,'\\usepackage[landscape]{geometry}\n');
fprintf(fid,'\\usepackage{array}\n');
fprintf(fid,'\\begin{document}\n');

fprintf(fid,'\\begin{table}\n');
fprintf(fid,'\\footnotesize\n');
fprintf(fid,'\\begin{tabular}{| ');
fprintf(fid,repmat('c|',[1,nbehaviors+1]));
fprintf(fid,'}\\hline\n');
fprintf(fid,'& \\parbox[t]{.33in}{\\centering %s} ',shortbehaviorsplot{:});
fprintf(fid,'\\\\\\hline\n');

maxc = .9;
minc = .2;
maxpplot = .05;
minpplot = .01;
for i = 1:ngenotypesplot,
  if strcmpi(shortlinenamesplot{i},'pBDPGAL4U'),
    continue;
  end
  fprintf(fid,'%s ',shortlinenamesplot{i});
  for behi = 1:nbehaviors,
    if idxtest(i,behi),
      if adjp_max(i,behi)<=fdr && adjp_max(i,behi)<=maxpplot,
        color1 = min(maxc,max(minc,(adjp_max(i,behi)-minpplot)/(maxpplot-minpplot)*(maxc-minc)+minc));
        if dir(i,behi) > 0,
          color = [1,color1,color1];
        else
          color = [color1,color1,1];
        end
        
        fprintf(fid,'& \\cellcolor[rgb]{%f,%f,%f}%.3f',color,adjp_max(i,behi));
      else
        fprintf(fid,'& %.3f ',adjp_max(i,behi));
      end
    else
      fprintf(fid,'& -- ');
    end
  end
  fprintf(fid,'\\\\\\hline\n');
end
fprintf(fid,'\\end{tabular}\n');
caption = ['Significance of behavioral differences for Fly Bowl assay. ',...
  'For each line, we compare the fraction of time the GAL4/TrpA line does each ',...
  'of 14 behaviors to two controls: the empty GAL4 crossed to TrpA ',...
  '(BDP/TrpA) and the GAL4/+. Each data point corresponds to the fraction of time ',...
  'that all 20 flies in a given video performed the given behavior. The table contains a dash if the ',...
  'behavioral differences to both controls had different signs. Otherwise, ',...
  'we use the Wilcoxon rank-sum test to compare to each control separately, ',...
  'adjust the p-value using the Benjamini-Hochberg correction for ',...
  'controlling the false discovery rate to 0.05, and report the ',...
  '(conservative) maximum adjusted p-value of the comparison to each control. ',...
  'Red indicates lines that performed the behavior significantly more than ',...
  'both controls, and blue indicates lines that performed the behavior ',...
  'significantly less than controls. The average number of videos (n) analyzed for each ',...
  'GAL4/TrpA line was ',sprintf('%.1f',mean(nexps_mbtrpg(1:end-1))),', for each GAL4/plus line ',...
  'was ',sprintf('%.1f',mean(nexps_mbplus(1:end-1))),', and for the pBDPGAL4U/TrpA line was ',...
  num2str(nexps_bdptrpg),'.'];
fprintf(fid,'\\caption{%s}',caption);
fprintf(fid,'\\end{table}\n');

fprintf(fid,'\\end{document}\n');

fclose(fid);
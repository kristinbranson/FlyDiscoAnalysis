
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
addpath /groups/branson/bransonlab/projects/olympiad/SAGE/MATLABInterface/Trunk;
addpath /groups/branson/bransonlab/projects/olympiad/anatomy/fileio;

%% parameters

expdirlist = '/groups/branson/home/bransonk/behavioranalysis/code/Jdetect/JAABA_guff/perframe/expdirs_jaabadetect20130606.txt';
control_line_names = {'pBDPGAL4U','FCF_pBDPGAL4U_1500437'};
main_control_line_name = 'pBDPGAL4U';
%requiredfiles = {'stats_perframe.mat'};
requiredfiles = {'scores_WingExtension.mat'};

% if the variance is 3x bigger, then it is not worth adding
maxfactorexp = sqrt(3);
minnframes_lowerbound = 200;

%% load in data

% list of experiments
expdirs = importdata(expdirlist);
nexps = numel(expdirs);

% get metadata
expnames = cell(size(expdirs));
for i = 1:nexps,
  [~,expnames{i}] = fileparts(expdirs{i});
end
expnames1 = cellfun(@(x) ['FlyBowl_',x],expnames,'UniformOutput',false);
metadata = SAGEListBowlExperiments('data_type',...
  {'ctrax_diagnostics_nframes_analyzed',...
  'QuickStats_BackSubStats_meanNConnComps'},...
  'checkflags',false,...
  'removemissingdata',false);
isin = ismember({metadata.file_system_path},expdirs);
metadata = metadata(isin);

[idx,metadataidx] = ismember(expdirs,{metadata.file_system_path});
ismissing = ~idx;

% try to read in data for missing exps
needsset = false(1,nexps);
for i = find(ismissing'),
  
  expdir = expdirs{i};
  metadatafile = fullfile(expdir,'Metadata.xml');
  if ~exist(metadatafile,'file'),
    continue;
  end
  try
    metadatacurr = ReadMetadataFile(metadatafile);
  catch 
    continue;
  end
  autochecksfile = fullfile(expdir,'automatic_checks_complete_results.txt');
  if exist(autochecksfile,'file'),
    try
      autochecksresults = ReadParams(autochecksfile);
    catch
      metadatacurr.automated_pf = 'U';
      metadatacurr.automated_pf_category = '';
    end
  else
    metadatacurr.automated_pf = 'U';
    metadatacurr.automated_pf_category = '';
  end
  [~,name] = fileparts(expdir);
  metadatacurr.experiment_name = ['FlyBowl_',name];
  metadatacurr.file_system_path = expdir;
  metadatacurr.line_name = metadatacurr.line;
  metadatacurr.manual_pf = 'U';
  metadatacurr.plate_bowl = sprintf('%d%s',metadatacurr.plate,metadatacurr.bowl);
  metadatacurr.rig_bowl = sprintf('%d%s',metadatacurr.rig,metadatacurr.bowl);
  needsset(i) = true;
    

  fnsremove = setdiff(fieldnames(metadatacurr),fieldnames(metadata));
  metadatacurr = rmfield(metadatacurr,fnsremove);
  metadata = structappend(metadata,metadatacurr);
  metadataidx(i) = numel(metadata);
  ismissing(i) = false;
  
end

% remove anything that could not be fixed
expdirs(ismissing) = [];
expnames(ismissing) = [];
expnames1(ismissing) = [];
metadataidx(ismissing) = [];
nexps = numel(expdirs);
metadata = metadata(metadataidx);

% add in missing sets
while true,
  if ~any(needsset),
    break;
  end
  i = find(needsset,1);
  idxcurr = find(strcmp({metadata.line_name},metadata(i).line_name));
  idxcurr = setdiff(idxcurr,i);
  if ~isempty(idxcurr),
    idxcurr1 = find((strcmp(cellfun(@(x) x(1:8),{metadata(idxcurr).exp_datetime},'UniformOutput',false),...
      metadata(i).exp_datetime(1:8))) & ([metadata(idxcurr).plate] == metadata(i).plate) & ...
      ~cellfun(@isempty,{metadata(idxcurr).set}));
    idxcurr = idxcurr(idxcurr1);
    if ~isempty(idxcurr),
      set = unique({metadata(idxcurr).set});
      if numel(set) == 1,
        metadata(i).set = set{1};
      else
        ds = datenum({metadata(idxcurr).exp_datetime},'yyyymmddTHHMMSS');
        d = datenum(metadata(i).exp_datetime,'yyyymmddTHHMMSS');
        [~,j] = min(abs(d-ds));
        metadata(i).set = set{j};
      end
    end
  end
  if isempty(idxcurr),
    set = sprintf('%s__Rig%d__%s',metadata(i).line_name,metadata(i).rig,metadata(i).exp_datetime);
    metadata(i).set = set;
  end
  needsset(i) = false;
end

% check for required files
missingfiles = false(1,nexps);
for i = 1:nexps,
  missingfiles(i) = ~all(cellfun(@(x) exist(fullfile(expdirs{i},x),'file'),requiredfiles));
end

expdirs(missingfiles) = [];
expnames(missingfiles) = [];
expnames1(missingfiles) = [];
metadata(missingfiles) = [];
nexps = numel(expdirs);

% find exps for which automated_pf and automated_pf_category disagree
isp1 = cellfun(@(x) isempty(x) || strcmpi(x,'null'),{metadata.automated_pf_category});
isp2 = strcmp({metadata.automated_pf},'P');

for i = find(isp1 ~= isp2),
  expdir = expdirs{i};
  autochecksfile = fullfile(expdir,'automatic_checks_complete_results.txt');
  if ~exist(autochecksfile,'file'),
    metadata(i).automated_pf = 'U';
    metadata(i).automated_pf_category = 'NULL';
  else
    try
      tmp = ReadParams(autochecksfile);
      metadata(i).automated_pf = tmp.automated_pf;
      if isfield(tmp,'automated_pf_category'),
        metadata(i).automated_pf_category = tmp.automated_pf_category;
      else
        metadata(i).automated_pf_category = 'NULL';
      end
    catch ME,
      warning('%d: %s: %s',i,expdir,getReport(ME));
    end
  end
end


% remove bad failures
badcategories = {...
  'bad_number_of_flies'
  'bad_number_of_flies_per_sex'
  'flag_aborted_set_to_1'
  'flag_redo_set_to_1'
  'fliesloaded_time_too_long'
  'fliesloaded_time_too_short'
  'short_video'
};

badexps = ismember({metadata.automated_pf_category},badcategories);

expdirs(badexps) = [];
expnames(badexps) = [];
expnames1(badexps) = [];
metadata(badexps) = [];
nexps = numel(expdirs);

% remove extraneous control line
badexps = strcmp({metadata.line_name},'FCF_pBDPGAL4U_1500437');
expdirs(badexps) = [];
expnames(badexps) = [];
expnames1(badexps) = [];
metadata(badexps) = [];
nexps = numel(expdirs);

% also remove exps with automated_pf = 'F' and notes_curation reflects that
% there are not enough flies
% TODO!


questionableexps = find(~strcmp({metadata.automated_pf},'P'));
dokeep = zeros(1,numel(questionableexps));
for ii = 1:numel(questionableexps),
  i = questionableexps(ii);
  fprintf('%s: automated_pf_category: %s\n',metadata(i).experiment_name,metadata(i).automated_pf_category);
  fprintf(['notes_curation:\n',metadata(i).notes_curation]);
  fprintf('\n');
  web('-browser',metadata(i).file_system_path);
  res = questdlg('Use this experiment?');
  if strcmp(res,'Cancel'),
    break;
  elseif strcmp(res,'Yes'),
    dokeep(ii) = 1;
  else
    dokeep(ii) = -1;
  end
end
badexps = questionableexps(dokeep<=0);
expdirs(badexps) = [];
expnames(badexps) = [];
expnames1(badexps) = [];
metadata(badexps) = [];
nexps = numel(expdirs);

% save to file
save CollectedPrimaryMetadata20130613.mat metadata;

%% load in the data


load CollectedPrimaryMetadata20130613.mat;
expdirs = {metadata.file_system_path};
nexps = numel(expdirs);
expnames = cell(size(expdirs));
for i = 1:nexps,
  [~,expnames{i}] = fileparts(expdirs{i});
end
expnames1 = cellfun(@(x) ['FlyBowl_',x],expnames,'UniformOutput',false);

%% load in stats

issuccess = false(1,nexps);
isquestionable = false(1,nexps);
statfns = {};
allstats = struct;
allstdstats = struct;
allnframestotal = struct;
for i = 1:nexps,

  expdir = expdirs{i};
  [~,name] = fileparts(expdir);
  fprintf('%d/%d: %s\n',i,nexps,name);
  

  statsfile = fullfile(expdir,'stats_perframe.mat');
  if ~exist(statsfile,'file'),
    fprintf('File %s does not exist\n',statsfile);
    continue;
  end
  
  statsdata = load(statsfile,'statsperexp');
  fnscurr = fieldnames(statsdata.statsperexp);
  newfns = setdiff(fnscurr,statfns);
  missingfns = setdiff(statfns,fnscurr);
  for j = 1:numel(newfns),
    fn = newfns{j};
    allstats.(fn) = nan(1,numel(expdirs));
    allstdstats.(fn) = nan(1,numel(expdirs));
    allnframestotal.(fn) = nan(1,numel(expdirs));
    fprintf('exp %d: Adding %s\n',i,fn);
  end
  if ~isempty(missingfns),
    fprintf('Experiment %s missing the following stats:\n',name);
    fprintf('%s\n',missingfns{:});
    isquestionable(i) = true;
  end
  
  statfns = union(statfns,newfns);
  for j = 1:numel(fnscurr),
    fn = fnscurr{j};
    allstats.(fn)(i) = statsdata.statsperexp.(fn).meanmean;
    allstdstats.(fn)(i) = statsdata.statsperexp.(fn).stdmean;
    allnframestotal.(fn)(i) = statsdata.statsperexp.(fn).Z;
  end
  
  issuccess(i) = true;
  
end

expdirs(~issuccess) = [];
expnames(~issuccess) = [];
expnames1(~issuccess) = [];
metadata(~issuccess) = [];
nexps = numel(expdirs);

nstats = numel(statfns);

%% remove experiments without enough frames

MINNFRAMESTOTAL = 400000;
nframestotal = allnframestotal.velmag_ctr_flyany_frameany;
badexps = nframestotal < MINNFRAMESTOTAL;

expdirs(badexps) = [];
expnames(badexps) = [];
expnames1(badexps) = [];
metadata(badexps) = [];
nexps = numel(expdirs);
allstats = structfun(@(x) x(~badexps),allstats,'UniformOutput',false);
allstdstats = structfun(@(x) x(~badexps),allstdstats,'UniformOutput',false);
allnframestotal = structfun(@(x) x(~badexps),allnframestotal,'UniformOutput',false);


%% save to file

save CollectedPrimaryPerFrameStats20130614.mat ...
  allnframestotal allstats allstdstats ...
  control_line_names expdirlist expdirs expnames expnames1 issuccess ...
  main_control_line_name metadata nexps nstats statfns;

%% compute set stats using all experiments

[set_names,firstidx,setidx] = unique({metadata.set});
nsets = numel(set_names);
setmeans = nan(nsets,nstats);
setnexps = nan(nsets,nstats);
setmetadata = metadata(firstidx);
for stati = 1:nstats,
  fn = statfns{stati};
  idxcurr = ~isnan(allstats.(fn)) & ~isinf(allstats.(fn));
  setmeans(:,stati) = accumarray(setidx(idxcurr)',allstats.(fn)(idxcurr)',[nsets,1],@mean);
  setnexps(:,stati) = hist(setidx(idxcurr),1:nsets);
end

%% compute line stats using all sets

[line_names,firstidx,set2lineidx] = unique({setmetadata.line_name});
[~,exp2lineidx] = ismember({metadata.line_name},line_names);
nlines = numel(line_names);
linemeans = nan(nlines,nstats);
linestds = nan(nlines,nstats);
linensets = nan(nlines,nstats);
linemetadata = setmetadata(firstidx);

for stati = 1:nstats,
  fn = statfns{stati};
  idxcurr = ~isnan(setmeans(:,stati)) & ~isinf(setmeans(:,stati));
  linemeans(:,stati) = accumarray(set2lineidx(idxcurr)',setmeans(idxcurr,stati)',[nlines,1],@mean);
  linestds(:,stati) = accumarray(set2lineidx(idxcurr)',setmeans(idxcurr,stati)',[nlines,1],@(x) std(x,1));
  linensets(:,stati) = hist(set2lineidx(idxcurr),1:nlines);
end

idxcontrol = find(strcmp(line_names,main_control_line_name));
controlstd = nanstd(setmeans(set2lineidx==idxcontrol,:),1,1);

%%

save CollectedPrimaryPerFrameStats20130614.mat ...
  allnframestotal allstats allstdstats control_line_names controlstd ...
  exp2lineidx expdirlist expdirs expnames expnames1 idxcontrol ...
  issuccess line_names linemeans linemetadata linensets ...
  main_control_line_name metadata nexps nlines nsets nstats ...
  set2lineidx set_names setidx setmeans setmetadata setnexps ...
  statfns;

%% how many frames do we need to compute a statistic?

%idxcurr = find(strcmp({metadata.line_name},main_control_line_name));

% statistics for which we should use every experiment, as all frames are
% used
%statidx = find(cellfun(@isempty,regexp(statfns,'^(fractime)|(boutfreq)','once')));

nbins = 50;
nframestotal = allnframestotal.velmag_ctr_flyany_frameany;
minnframes_upperbound = prctile(nframestotal,.1);
minnexps = 8;
linenexps = nan(nlines,nstats);
for stati = 1:nstats,
  linenexps(:,stati) = accumarray(set2lineidx(:),setnexps(:,stati),[nlines,1]);
end
lineidxuse = linenexps >= minnexps;

statidx = find(cellfun(@isempty,regexp(statfns,'(boutfreq)|(fractime)|(frameany)','once')));

hfig = 1;
figure(hfig);
clf;
hold on;
colors = jet(numel(statidx))*.75;

% plot all the same conditions on the same axes
statidxtmp = cellfun(@isempty,regexp(statfns,'^fractime','once'));
condition = cell(1,nstats);
condition(statidxtmp) = cellfun(@(x) x{1},regexp(statfns(statidxtmp),'_(fly.*_frame.*)$','tokens','once'),'UniformOutput',false);
condition(~statidxtmp) = cellfun(@(x) [x{1},'_frameany'],regexp(statfns(~statidxtmp),'_(fly.*)_','tokens','once'),'UniformOutput',false);
[unique_conditions,~,conditionidx] = unique(condition(statidx));
clf;
nr = floor(sqrt(numel(unique_conditions)));
nc = ceil(numel(unique_conditions)/nr);
hax = createsubplots(nr,nc,.05);
for i = 1:numel(unique_conditions),
  hold(hax(i),'on');
  title(hax(i),unique_conditions{i},'Interpreter','none')
end

% minimum number of frames required per experiment to have a reasonable
% statistic
minnframes = nan(1,numel(statfns));
std_coeffs = nan(1,numel(statfns));


for statii = randperm(numel(statidx)),
%for statii = 1:numel(statidx),
  
  stati = statidx(statii);
  fn = statfns{stati};
  
  expidxuse = lineidxuse(exp2lineidx,stati)';
  
  idxcurr1 = find(expidxuse & ~isnan(allstats.(fn)) & ~isinf(allstats.(fn)));

  % z-scored data
  mus_perline = accumarray(exp2lineidx(idxcurr1)',allstats.(fn)(idxcurr1)',[nlines,1],@mean,nan);
  sigs_perline = accumarray(exp2lineidx(idxcurr1)',allstats.(fn)(idxcurr1)',[nlines,1],@(x) std(x,1),nan);
  minsig = prctile(sigs_perline(~isnan(sigs_perline)),1);
  mus = mus_perline(exp2lineidx(idxcurr1));
  sig = max(minsig,sigs_perline(exp2lineidx(idxcurr1)));
  %sig = controlstd(stati);
  z = (allstats.(fn)(idxcurr1)'-mus)./sig;
%   medianz_perline = accumarray(exp2lineidx(idxcurr1)',abs(z),[nlines,1],@median,nan);
%   medianz_perexp = medianz_perline(exp2lineidx(idxcurr1));
  
  nframestotal = allnframestotal.(fn)(idxcurr1);
  idxignore = isnan(z') | isinf(z');

  minnframes_upperbound_curr = prctile(nframestotal,50);
  
  
%   xin = [1./sqrt(nframestotal(~idxignore))',1./nframestotal(~idxignore)'];
%   [coeffs,info] = lasso(xin,abs(z(~idxignore)),'CV',5);
%   lambdai = argmin(info.MSE);
%   coeffs = [info.Intercept(lambdai);coeffs(:,lambdai)];

  besterr = inf;
  for order = [.5],%,1,1.5],
    %xin = [ones(nnz(~idxignore),1),1./(nframestotal(~idxignore)).^order'];
    xin = 1./(nframestotal(~idxignore)).^order';
    [coeffs,err] = lsqnonneg(xin,abs(z(~idxignore)));
%     [coeffs,~,errs] = regress(abs(z(~idxignore)),xin);
%     err = sum(abs(errs));
    %disp(err-besterr)
    if err < besterr,
      besterr = err;
      bestorder = order;
      bestcoeffs = coeffs;
    end
  end
  coeffs = bestcoeffs;
  
  absz = abs(z(~isnan(z)));
  medianz = median(absz);
  minnframes_curr = (coeffs(1)./(medianz*maxfactorexp))^2;
  %minnframes_curr = (coeffs(2)/(medianz*maxfactorexp-coeffs(1)))^2;
  minnframes(stati) = max(minnframes_lowerbound,min(minnframes_curr,min(minnframes_upperbound,minnframes_upperbound_curr)));
  std_coeffs(:,stati) = coeffs;

  tmpx = linspace(min(prctile(nframestotal,.1),minnframes(stati)/2),prctile(nframestotal,95),100)';
  %tmpx1 = [ones(numel(tmpx),1),1./tmpx.^bestorder];
  tmpx1 = 1./tmpx.^bestorder;
  tmpy = tmpx1*coeffs;
  
  nframestotal1 = nframestotal(~isnan(z));
  [~,order] = sort(nframestotal1);
  [~,rank] = sort(order);
  bin = ceil(rank/numel(absz)*nbins);
  meanz = accumarray(bin(:),absz,[nbins,1],@mean);
  center = accumarray(bin(:),nframestotal1(:),[nbins,1],@mean);
  if meanz(1) < medianz*maxfactorexp,
    minnframes(stati) = minnframes_lowerbound;
  end
  fprintf('For %s, min nframes = %f (%f), median nframes = %f.\n',fn,minnframes(stati),minnframes_curr,median(nframestotal(~isnan(nframestotal))));  
  
%   clf;
%   plot(nframestotal,abs(z),'.','Color','k');
%   hold on;
%   plot(tmpx,tmpy,'-','Color','r','LineWidth',2);
%   xlim = get(gca,'XLim');
%   plot(xlim,nanmedian(abs(z(~idxignore)))+[0,0],'b-');
%   plot(center,meanz,'g.-');
%   ylim = get(gca,'YLim');
%   plot(minnframes(stati)+[0,0],ylim,'m--','LineWidth',5);
%   [~,order] = sort(nframestotal);
%   bin = ceil(order/numel(nframestotal)*nbins);
%   meanabsz = accumarray(bin(:),abs(z(:)),[nbins,1],@mean,nan);
%   xlabel('Number of frames analyzed in the experiment');
%   ylabel('Stds from the mean');
%   title(fn,'Interpreter','none');
%     input(fn);

  axi = conditionidx(statii);
  name = regexp(fn,sprintf('^(.*)_%s$',unique_conditions{conditionidx(statii)}),'once','tokens');
  if isempty(name),
    name = fn;
  else
    name = name{1};
  end
  plot(hax(axi),tmpx,tmpy,'-','Color',colors(statii,:),'LineWidth',2);
  plot(hax(axi),minnframes(stati)+[0,0],[0,max(tmpy)],'--','Color',colors(statii,:));
  text(tmpx(1),tmpy(1),name,'Interpreter','none','Color',colors(statii,:),'HorizontalAlignment','left','Parent',hax(axi));
  drawnow;
    
end

minnframes = round(minnframes);

save -append CollectedPrimaryPerFrameStats20130614.mat minnframes statidx std_coeffs;

%% remove experiments from statistics for which they don't have enough frames

[set_names,firstidx,setidx] = unique({metadata.set});
nsets = numel(set_names);

setstats = struct('means',struct,'stds',struct,'nexps',struct,'metadata',struct);
setstats.metadata = metadata(firstidx);

for stati = 1:numel(statfns),

  statfn = statfns{stati};
  
  % remove experiments that have nans or infs
  idxcurr = ~isnan(allstats.(statfn)) & ~isinf(allstats.(statfn));
  
  % remove experiments that don't have nframestotal
  nframestotal = allnframestotal.(statfn)(idxcurr);
  badexp1 = isnan(nframestotal);
  idxcurr(idxcurr) = ~badexp1;
    
  % remove experiments that don't have enough frames
  if ismember(stati,statidx),
    badexp1 = nframestotal < minnframes(stati);
    fprintf('%s: removing %d experiments that don''t have enough frames\n',statfn,nnz(badexp1));
    idxcurr(idxcurr) = ~badexp1;
  end

  % count number of experiments per set
  setstats.nexps.(statfn) = hist(setidx(idxcurr),1:numel(set_names));

  % compute statistics
  x = allstats.(statfn)(idxcurr);
  setstats.means.(statfn) = accumarray(setidx(idxcurr)',x',[nsets,1],@mean)';
  setstats.stds.(statfn) = accumarray(setidx(idxcurr)',x',[nsets,1],@(y) std(y,1))';
    
end

%% update line stats

% < 2 experiments will increase the variance of the estimate, assuming most
% sets 
minnexps = 2;

linestats = struct('means',struct,'stds',struct,'nexps',struct,'metadata',struct);
[linestats.line_names,firstidx,set2lineidx] = unique({setstats.metadata.line});
nlines = numel(linestats.line_names);
for stati = 1:numel(statfns),
  
  statfn = statfns{stati};
  
  % remove sets that don't have enough experiments or have bad values
  idxcurr = setstats.nexps.(statfn) >= minnexps & ...
    ~isnan(setstats.means.(statfn)) & ...
    ~isinf(setstats.means.(statfn));
  
  fprintf('%s: removing %d sets that don''t have enough experiments\n',statfn,nnz(~idxcurr));
  
  % count number of sets for each line
  linestats.nsets.(statfn) = hist(set2lineidx(idxcurr),1:numel(linestats.line_names));
  for i = 0:10,
    n1 = nnz(linensets(:,stati)==i);
    n2 = nnz(linestats.nsets.(statfn)==i);
    if n1 == 0 && n2 == 0,
      continue;
    end
    fprintf('Nlines with %d sets: %d -> %d\n',i,n1,n2);
  end

  % compute statistics
  x = setstats.means.(statfn)(idxcurr);
  linestats.means.(statfn) = accumarray(set2lineidx(idxcurr)',x',[nlines,1],@mean)';
  linestats.stds.(statfn) = accumarray(set2lineidx(idxcurr)',x',[nlines,1],@(y) std(y,1))';
  linestats.means.(statfn)(linestats.nsets.(statfn)==0) = nan;
  linestats.stds.(statfn)(linestats.nsets.(statfn)==0) = nan;

end

%% save

save -append CollectedPrimaryPerFrameStats20130614.mat setstats linestats nlines set2lineidx setidx nsets;

%% subtract off control means for the surrounding month

controldateradius = 15;

% which sets are control
setiscontrol = strcmp({setstats.metadata.line_name},main_control_line_name);
setidxcontrol = find(setiscontrol);

setdatenum = floor(datenum({setstats.metadata.exp_datetime},'yyyymmddTHHMMSS'));
controldatenum = setdatenum(setiscontrol);
mincontroldatenum = min(controldatenum);
maxcontroldatenum = max(controldatenum);

setstats.normmeans = struct;
setstats.controlmeans = struct;
setstats.nsetscontrolnorm = struct;
for seti = 1:nsets,
  
  fprintf('Set %s %d / %d\n',setstats.metadata(seti).set,seti,nsets);
  
  mindatenumcurr = max(mincontroldatenum,setdatenum(seti)-controldateradius);
  maxdatenumcurr = mindatenumcurr+controldateradius*2+1;
  if maxdatenumcurr > maxcontroldatenum,
    maxdatenumcurr = maxcontroldatenum;
    mindatenumcurr = maxdatenumcurr-controldateradius*2+1;
  end
  
  idxdate = (controldatenum>=mindatenumcurr & controldatenum <= maxdatenumcurr)';
  
  for stati = 1:nstats,
    
    statfn = statfns{stati};
    
    if seti == 1,
      setstats.normmeans.(statfn) = nan(1,nsets);
    end
    
    idxcurr = setidxcontrol(idxdate & setstats.nexps.(statfn)(setiscontrol) >= minnexps);
    idxcurr = setdiff(idxcurr,seti);
    setstats.controlmeans.(statfn)(seti) = mean(setstats.means.(statfn)(idxcurr));
    setstats.nsetscontrolnorm.(statfn)(seti) = nnz(idxcurr);
    setstats.normmeans.(statfn)(seti) = setstats.means.(statfn)(seti) - setstats.controlmeans.(statfn)(seti);

  end
end

%% update line stats

for stati = 1:numel(statfns),
  
  statfn = statfns{stati};
  
  % remove sets that don't have enough experiments or have bad values
  idxcurr = setstats.nexps.(statfn) >= minnexps & ...
    ~isnan(setstats.means.(statfn)) & ...
    ~isinf(setstats.means.(statfn));

  % compute statistics
  x = setstats.normmeans.(statfn)(idxcurr);
  linestats.normmeans.(statfn) = accumarray(set2lineidx(idxcurr)',x',[nlines,1],@mean,nan)';

end

%% save results

save -append CollectedPrimaryPerFrameStats20130614.mat setstats linestats nlines set2lineidx setidx nsets;


%% compare to randomly sampled control data

% which sets are control
setiscontrol = strcmp({setstats.metadata.line_name},main_control_line_name);
setidxcontrol = find(setidxcontrol);

% experiment to set index
[~,exp2setidx] = ismember({metadata.set},{setstats.metadata.set});

% experiments that are control
expiscontrol = find(ismember(exp2setidx,setiscontrol));

% number of control sets
ncontrolsets = nnz(setiscontrol);

% number of experiments in each control set
statfn = 'velmag_flyany_frameany';
setnexps_control = setstats.nexps.(statfn)(setiscontrol)';

controldatenum = floor(datenum({setstats.metadata(setiscontrol).exp_datetime},'yyyymmddTHHMMSS'));

nsamples = 1000;

fracsmaller = nan(nlines,nstats);
fracbigger = nan(nlines,nstats);

maxnexps = max(setstats.nexps.(statfn));

minnframes(isnan(minnframes)) = 0;

save ComputePValueBySamplingData20130616.mat ...
  statfns nstats setstats allstats allnframestotal minnframes minnexps ...
  setiscontrol nsets maxnexps exp2setidx nlines linestats set2lineidx ...
  nsamples;


%% run jobs on the cluster
% 
% 
% for stati = stati:nstats,
%   
%   statfn = statfns{stati};
%   
%   fprintf('Statistic %s: %d / %d\n',statfn,stati,nstats);
%   
%   isgoodset = ~isinf(setstats.normmeans.(statfn)) & ...
%     ~isnan(setstats.normmeans.(statfn)) & ...
%     setstats.nexps.(statfn) >= minnexps;
%   isgoodexp = ~isinf(allstats.(statfn)) & ...
%     ~isnan(allstats.(statfn)) & ...
%     allnframestotal.(statfn) >= minnframes(stati);
%   
%   isgoodcontrolset = isgoodset & setiscontrol;
%   if nnz(isgoodcontrolset <= 1),
%     fprintf('One or less good control sets for statistic %s, not computing p-values\n',statfn);
%     continue;
%   end
%   
%   
%   idxgoodcontrolset = find(isgoodcontrolset);
%   
%   set2expcurr = nan(nsets,maxnexps);
%   for seti = idxgoodcontrolset,
%     tmp = find(exp2setidx==seti & isgoodexp);
%     set2expcurr(seti,1:numel(tmp)) = tmp;
%   end
% 
%   parfor linei = 1:nlines,
% 
%     if mod(linei,100) == 0,
%       fprintf('Stat %s (%d / %d), line %s (%d / %d)\n',statfn,stati,nstats,...
%         linestats.line_names{linei},linei,nlines);
%     end
%     
%     %fprintf('Line %s: %d / %d\n',linestats.line_names{linei},linei,nlines);
%   
%     % number of experiments in each set
%     setidxcurr = find(set2lineidx==linei & isgoodset);
%     if isempty(setidxcurr),
%       %fprintf('No good sets found for line %s, skipping.\n',linestats.line_names{linei});
%       continue;
%     end
%     
%     nexpscurr = setstats.nexps.(statfn)(setidxcurr);
%     datenumscurr = setdatenum(set2lineidx==linei & isgoodset);
%     nexpscurr = sort(nexpscurr,'descend');
%     nsetscurr = numel(nexpscurr);
%     
%     % sample
%     setnormmu = nan(nsamples,nsetscurr);
%     for setii = 1:nsetscurr,
%       
%       setidxallowed = find(isgoodcontrolset & ...
%         setstats.nexps.(statfn) >= nexpscurr(setii));
%       for nsub = 1:nexpscurr(setii)-1,
%         if numel(setidxallowed) > 1,
%           break;
%         end
%         setidxallowed = find(isgoodcontrolset & ...
%           setstats.nexps.(statfn) >= nexpscurr(setii)-nsub);
%       end
%       nsub = nsub-1;
%       if nsub > 0,
%         fprintf('For stat %s, line %s, set %d, needed to consider sets %d exps smaller than than this set\n',statfn,linestats.line_names{linei},setii,nsub);
%         nexpscurr(setii) = nexpscurr(setii)-nsub;
%       end
%       setsampleis = randsample(setidxallowed,nsamples,true);
%             
%       % choose the experiments per set
%       controlnexpscurr = setstats.nexps.(statfn)(setsampleis);
%       
%       % by default, use the first experiments
%       expsampleis = set2expcurr(setsampleis,1:nexpscurr(setii));
%       
%       % for sets with more experiments, sample without replacement
%       tmpidx = find(controlnexpscurr > nexpscurr(setii));
%       for tmpi = tmpidx,
%         expsampleis(tmpi,:) = set2expcurr(setsampleis(tmpi),randsample(controlnexpscurr(tmpi),nexpscurr(setii)));
%       end
%       
%       tmp = allstats.(statfn)(expsampleis);
%       tmp = reshape(tmp,[nsamples,nexpscurr(setii)]);
%       mu = mean(tmp,2);
%       setnormmu(:,setii) = mu - setstats.controlmeans.(statfn)(setsampleis)';
% 
%     end
% 
%     mu = mean(setnormmu,2);
%     fracsmaller(linei,stati) = nnz(mu<linestats.normmeans.(statfn)(linei))/nsamples;
%     fracbigger(linei,stati) = nnz(mu>linestats.normmeans.(statfn)(linei))/nsamples;
%   
%   end
%   
% end

%% load results

outdir = 'ComputePValueBySamplingResults20130616';

fracsmaller = nan(nlines,nstats);
fracbigger = nan(nlines,nstats);
for stati = 1:nstats,
  
  filename = fullfile(outdir,sprintf('PvaluesForStat%03d.mat',stati));
  if ~exist(filename,'file'),
    fprintf('Stat %d failed\n',stati);
    continue;
  end
  
  tmp = load(filename);
  fracsmaller(1:end-1,stati) = tmp.fracsmaller_stat;
  fracbigger(1:end-1,stati) = tmp.fracbigger_stat;
  
end

save -append CollectedPrimaryPerFrameStats20130614.mat fracsmaller fracbigger

%% compute standard deviation of control line for each statistic
controllinestd = nan(1,nstats);
controllinemean = nan(1,nstats);
for stati = 1:nstats,
  ComputeControlLineStdBySampling;
  controllinemean(stati) = mucurr;
  controllinestd(stati) = sigcurr;
end

%% compute upper and lower bounds on these p-values at 99% confidence

alpha = .05;
[tmp1,tmp2] = binofit(fracsmaller(:)*nsamples,nsamples,alpha);
tmp1 = reshape(tmp1,[nlines,nstats]);
fracsmaller_upperbound = reshape(tmp2(:,2),[nlines,nstats]);

[tmp1,tmp2] = binofit(fracbigger(:)*nsamples,nsamples,alpha);
tmp1 = reshape(tmp1,[nlines,nstats]);
fracbigger_upperbound = reshape(tmp2(:,2),[nlines,nstats]);

%% cluster some statistics

statfnscurr = {
  'fractime_flyany_frameattemptedcopulation'
  'fractime_flyany_framebackup'
  'fractime_flyany_framebodyturn'
  'fractime_flyany_framechase'
  'fractime_flyany_framecopulation'
  'fractime_flyany_framecrabwalkall'
  'fractime_flyany_framecrabwalkextreme'
  'fractime_flyany_framejump'
  'fractime_flyany_framenotanybehavior'
  'fractime_flyany_framepivotcenter'
  'fractime_flyany_framepivottail'
  'fractime_flyany_framerighting'
  'fractime_flyany_framestop'
  'fractime_flyany_frametouch'
  'fractime_flyany_framewalk'
  'fractime_flyany_framewingextension'
  'fractime_flyany_framewingflick'
  'fractime_flyany_framewinggrooming'

  'fractime_flyany_framechase_notwingextension'
  'fractime_flyany_framestop_notwinggrooming'
  'fractime_flyany_frametouch_notchase'
  'fractime_flyany_framewingextension_notchase'
  
  'fractime_flyany_framebackup_nearfly'
  'fractime_flyany_framebackup_nearwall'
  'fractime_flyany_framebackup_notnearfly_notnearwall'
  
  'fractime_flyany_framecrabwalkextreme_nearfly'
  'fractime_flyany_framecrabwalkextreme_nearwall'
  'fractime_flyany_framecrabwalkextreme_notnearfly_notnearwall'
  
  
  'fractime_flyany_framejump_nearfly'
  'fractime_flyany_framejump_nearwall'
  'fractime_flyany_framejump_notnearfly_notnearwall'
  
  
  'fractime_flyany_framemove_nearfly'
  'fractime_flyany_framemove_nearwall'
  'fractime_flyany_framemove_notnearfly_notnearwall'
  
  
  'fractime_flyany_framepivotcenter_nearfly'
  'fractime_flyany_framepivotcenter_nearwall'
  'fractime_flyany_framepivotcenter_notnearfly_notnearwall'
  
  'fractime_flyany_framepivottail_nearfly'
  'fractime_flyany_framepivottail_nearwall'
  'fractime_flyany_framepivottail_notnearfly_notnearwall'
  
  
  'fractime_flyany_framerighting_nearfly'
  'fractime_flyany_framerighting_nearwall'
  'fractime_flyany_framerighting_notnearfly_notnearwall'
  
  'fractime_flyany_framestop_nearfly'
  'fractime_flyany_framestop_nearwall'
  'fractime_flyany_framestop_notnearfly_notnearwall'
  
  'fractime_flyany_framewalk_nearfly'
  'fractime_flyany_framewalk_nearwall'
  'fractime_flyany_framewalk_notnearfly_notnearwall'
  
  'fractime_flyfemale_framebackup'
  'fractime_flyfemale_framebodyturn'
  'fractime_flyfemale_framechase'
  'fractime_flyfemale_framecrabwalkextreme'
  'fractime_flyfemale_framejump'
  'fractime_flyfemale_framenotanybehavior'
  'fractime_flyfemale_framepivotcenter'
  'fractime_flyfemale_framepivottail'
  'fractime_flyfemale_framerighting'
  'fractime_flyfemale_framestop'
  'fractime_flyfemale_frametouch'
  'fractime_flyfemale_framewalk'
  'fractime_flyfemale_framewingflick'
  'fractime_flyfemale_framewinggrooming'
  
  'fractime_flymale_frameattemptedcopulation'
  'fractime_flymale_framebackup'
  'fractime_flymale_framebodyturn'
  'fractime_flymale_framechase'
  'fractime_flymale_framecrabwalkextreme'
  'fractime_flymale_framejump'
  'fractime_flymale_framenotanybehavior'
  'fractime_flymale_framepivotcenter'
  'fractime_flymale_framepivottail'
  'fractime_flymale_framerighting'
  'fractime_flymale_framestop'
  'fractime_flymale_frametouch'
  'fractime_flymale_framewalk'
  'fractime_flymale_framewingextension'
  'fractime_flymale_framewingflick'
  'fractime_flymale_framewinggrooming'
  
  'absangle2wall_flyany_framemove'
  'absangle2wall_flyany_framenearwall'
  'absangle2wall_flyany_framestop'
  
  'absanglefrom1to2_nose2ell_flyany_framenearfly'
  'absanglefrom1to2_nose2ell_flyany_framestop_nearfly'
  'absanglefrom1to2_nose2ell_flyany_frametouch'
  
  'absdangle2wall_flyany_framenearwall'
  
  'absdtheta_flyany_frameany'
  'absdtheta_flyany_framepivotcenter'
  'absdtheta_flyany_framepivottail'
  'absdtheta_flyany_framebodyturn'
  'absdtheta_flyany_framewalk'
  
  'absdv_cor_flyany_framecrabwalkall'
  'absdv_cor_flyany_framecrabwalkextreme'
  'absdv_cor_flyany_framenearfly'
  'absdv_cor_flyany_framenearwall'
  'absdv_cor_flyany_framenotnearfly_notnearwall'
  'absdv_cor_flyany_framewalk'
  
  'absphidiff_nose2ell_flymale_framechase'
  
  'absthetadiff_nose2ell_flyany_framenearfly'
  'absthetadiff_nose2ell_flyany_frametouch'
  
  'absthetadiff_nose2ell_flyfemale_frametouch'
  'absthetadiff_nose2ell_flymale_frametouch'
  
  'absthetadiff_nose2ell_flyfemale_framenearfly'
  'absthetadiff_nose2ell_flymale_framenearfly'
  
  'absthetadiff_nose2ell_flymale_framechase'
  'absthetadiff_nose2ell_flymale_framewingextension'
  
  'absyaw_flyany_framemove'
  
  'angleonclosestfly_flyany_framenearfly'
  'angleonclosestfly_flyany_framestop_nearfly'
  'angleonclosestfly_flyany_frametouch'
  
  'angleonclosestfly_flyfemale_framenearfly'
  'angleonclosestfly_flyfemale_frametouch'
  
  'angleonclosestfly_flymale_framechase'
  'angleonclosestfly_flymale_framenearfly'
  'angleonclosestfly_flymale_frametouch'
  'angleonclosestfly_flymale_framewingextension'
  
  'anglesub_flyany_frameany'
  
  'corfrac_maj_flyany_framepivotcenter'
  'corfrac_maj_flyany_framepivottail'
  
  'dangle2wall_flyany_framenearwall'
  'danglesub_flyany_framenearfly'
  
  'darea_flyany_frameany'
  
  'dcenter_flyany_frameany'
  'dcenter_flyany_framemove'
  'dcenter_flyany_framewingflick'
  'dcenter_flyany_framestop'
  
  'dcenter_flyfemale_frameany'
  'dcenter_flymale_frameany'
  
  'ddcenter_flyany_framenearfly'
  'ddist2wall_flyany_framenearwall'
  
  'ddnose2ell_flyany_framenearfly'
  'ddnose2ell_flymale_framechase'
    
  'dell2nose_flyany_frameany'
  'dell2nose_flymale_frameany'
  'dell2nose_flyfemale_frameany'
  'dell2nose_flyany_framewingflick'
  
  'dist2wall_flyany_frameany'
  'dist2wall_flyany_framemove'
  'dist2wall_flyany_framestop'
  
  'dist2wall_flyany_framewalk'
  
  'dist2wall_flyfemale_frameany'
  'dist2wall_flymale_frameany'
  
  'dmax_wing_angle_flyany_frameany'
  
  'dnose2ell_angle_30tomin30_flyany_frameany'
  'dnose2ell_angle_30tomin30_flyfemale_frameany'
  'dnose2ell_angle_30tomin30_flymale_frameany'
  
  'dnose2ell_flyany_frameany'
  
  'dnose2ell_flyany_framestop'
  'dnose2ell_flyany_frametouch'
  'dnose2ell_flyany_framewalk'
  
  'dnose2ell_flyfemale_frameany'
  'dnose2ell_flymale_frameany'
  
  'dnose2tail_flyany_frameany'
  'dnose2tail_flyany_framemove'
  'dnose2tail_flyfemale_frameany'
  'dnose2tail_flymale_frameany'
  
  'dtheta_flyany_frameany'
  
  'du_ctr_flyany_frameany'
  'du_ctr_flyany_framebackup'
  
  'du_ctr_flyany_framewalk'
  'du_ctr_flyfemale_frameany'
  'du_ctr_flymale_frameany'
  'du_ctr_flymale_framechase'
  
  'duration_flyany_framebackup'
  'duration_flyany_framebodyturn'
  'duration_flyany_framecrabwalkextreme'
  'duration_flyany_framejump'
  'duration_flyany_framemove'
  'duration_flyany_framepivotcenter'
  'duration_flyany_framepivottail'
  'duration_flyany_framerighting'
  'duration_flyany_framestop'
  'duration_flyany_frametouch'
  'duration_flyany_framewalk'
  'duration_flyany_framewingextension'
  'duration_flyany_framewingflick'
  'duration_flyany_framewinggrooming'
  
  'duration_flymale_frameattemptedcopulation'
  'duration_flymale_framechase'
  
  'dwing_angle_diff_flyany_frameany'
  
  'max_absdwing_angle_flyany_frameany'
  'max_absdwing_angle_flyany_framewingextension'
  'max_absdwing_angle_flyany_framewingflick'
  'max_absdwing_angle_flyany_framewinggrooming'
  
  'max_wing_angle_flyany_frameany'
  'max_wing_angle_flyany_framewingextension'
  'max_wing_angle_flyany_framewingflick'
  'max_wing_angle_flyany_framewinggrooming'
  
  'nflies_close_flyany_frameany'
  'nflies_close_flyany_framestop'
  'nflies_close_flyany_framewalk'
  'nflies_close_flyfemale_frameany'
  'nflies_close_flymale_frameany'
  'nflies_close_flymale_framechase'
  
  'velmag_ctr_flyany_frameany'
  'velmag_ctr_flyany_framejump'
  'velmag_ctr_flyany_framenearfly'
  'velmag_ctr_flyany_framenearwall'
  'velmag_ctr_flyany_framenotnearfly_notnearwall'
  'velmag_ctr_flyany_framewalk'
  'velmag_ctr_flyfemale_frameany'
  'velmag_ctr_flymale_frameany'
  'velmag_ctr_flymale_framechase'
  
  'veltoward_nose2ell_flyany_framenearfly'
  'veltoward_nose2ell_flymale_framechase'
  
  'wing_angle_diff_flyany_frameany'
  
  'wing_angle_diff_flyany_framewingextension'
  
  'wing_angle_imbalance_flyany_frameany'
  
  'wing_anglel_flyany_frameany'
  'wing_angler_flyany_frameany'
  'yaw_flyany_framemove'
  };

statidxcurr = find(ismember(statfns,statfnscurr));
statfnscurr = statfns(statidxcurr);
nstatscurr = numel(statidxcurr);

datacluster = nan(nlines,nstatscurr);
for ii = 1:nstatscurr,
  i = statidxcurr(ii);
  datacluster(:,ii) = linestats.normmeans.(statfns{i});
end

setiscontrol = strcmp({setstats.metadata.line_name},main_control_line_name);
signorm = nan(1,nstatscurr);
for ii = 1:nstatscurr,
  i = statidxcurr(ii);
  signorm(ii) = nanstd(setstats.normmeans.(statfns{i}),1);
end

zdatacluster = bsxfun(@rdivide,datacluster,signorm);

maxfraclinesmissingdata = .25;
statidxremove = find(sum(isnan(zdatacluster),1) >= nlines*maxfraclinesmissingdata);
statfnscurr0 = statfnscurr;
statidxcurr0 = statidxcurr;
datacluster0 = datacluster;
zdatacluster0 = zdatacluster;

statfnscurr(statidxremove) = [];
statidxcurr(statidxremove) = [];
datacluster(:,statidxremove) = [];
zdatacluster(:,statidxremove) = [];
signorm(statidxremove) = [];
nstatscurr = numel(statidxcurr);


shortstatnames = statfnscurr;
shortstatnames = regexprep(shortstatnames,'_flyany','');
shortstatnames = regexprep(shortstatnames,'^(.*)_fly(.*)_(.*)','$1_$3_$2');
shortstatnames = regexprep(shortstatnames,'^fractime_frame','fractime_');
shortstatnames = regexprep(shortstatnames,'^duration_frame','duration_');
shortstatnames = regexprep(shortstatnames,'_frameany','');
shortstatnames = regexprep(shortstatnames,'frame','');

shortlinenames = linestats.line_names;
shortlinenames = regexprep(shortlinenames,'GMR_','R');
shortlinenames = regexprep(shortlinenames,'_AE_01','');
shortlinenames = regexprep(shortlinenames,'_AD_01','D');

% compute pairwise distance between lines, ignoring entries for which
% either has nan
lined = zeros(nlines,nlines);
for linei = 1:nlines,
  for linej = linei+1:nlines,
    dcurr = abs(zdatacluster(linei,:)-zdatacluster(linej,:));
    lined(linei,linej) = nanmean(dcurr);
    lined(linej,linei) = lined(linei,linej);
  end
end
linedvec = squareform(lined,'tovector');

% compute pairwise distance between stats, ignoring entries for which
% either has nan
statd = zeros(nstatscurr,nstatscurr);
for stati = 1:nstatscurr,
  ignorei = isnan(zdatacluster(:,stati));
  for statj = stati+1:nstatscurr,
    ignorecurr = ignorei | isnan(zdatacluster(:,statj));
    if all(ignorecurr),
      statd(stati,statj) = nan;
    else
      r = corrcoef(zdatacluster(~ignorecurr,stati),zdatacluster(~ignorecurr,statj));
      statd(stati,statj) = 1 - abs(r(1,2));
    end
    statd(statj,stati) = statd(stati,statj);
  end
end
statdvec = squareform(statd,'tovector');



cgobj = clustergram(zdatacluster',...
  'RowLabels',shortstatnames,...
  'ColumnLabels',shortlinenames,...
  'Standardize','none',...
  'Cluster','all',...
  'RowPDist',statdvec,...
  'ColumnPDist',linedvec,...
  'Linkage','average',...
  'OptimalLeafOrder',true,...
  'ImputeFun',@ClustergramImputeFun);
  
  
%% false discovery rate correction

min_pvalue_observable = 2/nsamples;
pvalue = min(1,min(fracsmaller_upperbound(:,statidxcurr),fracbigger_upperbound(:,statidxcurr))*2);
q = .01;
goodidx = ~isnan(pvalue);
[issig,crit_p,tmp] = fdr_bh(pvalue(goodidx),q,'dep','yes');
adj_pvalue = nan(size(pvalue));
adj_pvalue(goodidx) = tmp;

%% find clustering of features and data

maxnlevels = 10;
errtol = .01;

clusterid_perlevel = cell(1,maxnlevels+1);
stats_perlevel = cell(1,maxnlevels+1);
clusterid_perlevel{1} = ones(1,nlines);
stats_perlevel{1} = {[]};
clusterisdone_perlevel = cell(1,maxnlevels+1);
clusterisdone_perlevel{1} = false;

nstatscurr0 = numel(statfnscurr0);

ncv = 10;

penaltynan = max(abs(zdatacluster0),[],1);

for level = 1:maxnlevels,
  
  nclustersprev = max(clusterid_perlevel{level});

  fprintf('Clustering stats and lines at level %d, nclusters previously = %d\n',level,nclustersprev);

  clusterids_next = zeros(1,nlines);
  
  clusterinext_off = 0;
  
  if all(clusterisdone_perlevel{level}),
    break;
  end
  
  clusterid_perlevel{level+1} = nan(1,nlines);
  
  for clusteri = 1:nclustersprev,
        
    fprintf('Analyzing cluster %d of lines...\n',clusteri);
    
    linescurr = find(clusterid_perlevel{level} == clusteri);
    nlinescurr = numel(linescurr);
    statsprev = stats_perlevel{level}{clusteri};
    inputs0 = false(1,nstatscurr0);
    inputs0(statsprev) = true;
   
    if clusterisdone_perlevel{level}(clusteri),
      
      fprintf('Clustering is done for this cluster, skipping.\n');
      
      % just copy
      clusterid_perlevel{level+1}(linescurr) = clusterinext_off+1;
      stats_perlevel{level+1}{clusterinext_off+1} = stats_perlevel{level}{clusteri};
      clusterisdone_perlevel{level+1}(clusterinext_off+1) = true;

      clusterinext_off = clusterinext_off + 1;
      continue;

    end 

    
    minerr = inf; 
    isbadline0 = any(isnan(zdatacluster0(linescurr,inputs0)),2);
    xcurr0 = [ones(nlinescurr,1),zdatacluster0(linescurr,inputs0)];
    
    for stati = 1:nstatscurr0,
      if inputs0(stati),
        continue;
      end
      
      inputs = inputs0;
      inputs(stati) = true;
      
      isbadline = isbadline0 | isnan(zdatacluster0(linescurr,stati));
      
      if nnz(~isbadline) <= 2,
        continue;
      end      

      xcurr = [xcurr0,zdatacluster0(linescurr,stati)];
      
      outputs = ~inputs;
      err = 0;
      coeffs = nan(nnz(inputs)+1,nstatscurr0);
      for statj = find(outputs),
        isbadlinecurr = isbadline|isnan(zdatacluster0(linescurr,statj));
        if nnz(~isbadlinecurr) <= 2,
          errcurr = 0;
        else
          [coeffs(:,statj),~,errcurr] = regress(zdatacluster0(linescurr(~isbadlinecurr),statj),xcurr(~isbadlinecurr,:));
        end
        err = err + sum(errcurr.^2) + nnz(isbadlinecurr)*penaltynan(statj)^2;
      end
      
      if err < minerr,
        minerr = err;
        beststati = stati;
        bestcoeffs = coeffs;
      end
    end

    fprintf('Chose to add statistic %s\n',statfnscurr0{beststati});
    
    % does including this statistic improve the regression fit?

    % compute cross-validation error    
    xprev = [ones(nlinescurr,1),zdatacluster0(linescurr,inputs0)];
    xcurr = [xprev,zdatacluster0(linescurr,beststati)];
    isbadline = isbadline0 | isnan(zdatacluster0(linescurr,beststati));
    
    tmp = ceil((1:nnz(~isbadline))/nnz(~isbadline)*ncv)';
    tmp = tmp(randperm(nnz(~isbadline)));
    cvindex = nan(nlinescurr,1);
    cvindex(~isbadline) = tmp;

    
    outputs = ~inputs0;
    outputs(beststati) = false;

    errprev = 0;
    errcurr = 0;
    
    for cvi = 1:ncv,
      idxcurr = cvindex == cvi;
            
      if nnz(~isbadline&~idxcurr) <= 2,
        continue;
      end
      
      for statj = find(outputs),
        isbadlinecurr = isbadline|isnan(zdatacluster0(linescurr,statj));
        if ~any(~isbadlinecurr&idxcurr)
          continue;
        end
        
        if nnz(~isbadlinecurr&~idxcurr) <= 2,
          continue;
        end
        
        [coeffs] = regress(zdatacluster0(linescurr(~isbadlinecurr&~idxcurr),statj),xprev(~isbadlinecurr&~idxcurr,:));
        errprev1 = zdatacluster0(linescurr(~isbadlinecurr&idxcurr),statj)-xprev(~isbadlinecurr&idxcurr,:)*coeffs;
        errprev = errprev + sum(errprev1.^2);
        [coeffs] = regress(zdatacluster0(linescurr(~isbadlinecurr&~idxcurr),statj),xcurr(~isbadlinecurr&~idxcurr,:));
        errcurr1 = zdatacluster0(linescurr(~isbadlinecurr&idxcurr),statj)-xcurr(~isbadlinecurr&idxcurr,:)*coeffs;
        errcurr = errcurr + sum(errcurr1.^2);
      end
      
    end
    
    isstatschange = errcurr < errprev*(1-errtol);

    if isstatschange,
      fprintf('Using statistic %s improves cross-validation error by %f%%, adding.\n',statfnscurr0{beststati},(errprev-errcurr)/errprev*100);      
      statscurr = [statsprev,beststati];
    else
      fprintf('Using statistic %s does not sufficiently improve cross-validation error (change = %f%%), not adding.\n',statfnscurr0{beststati},(errprev-errcurr)/errprev*100);
      statscurr = statsprev;
    end
    
    xcurr = zdatacluster0(linescurr,statscurr);
    isbadline = any(isnan(xcurr),2);
    
    % cluster using just these statistics
    fprintf('Clustering...\n');
    [clusteridx1,clustercenters,rsscurr] = mykmeans(xcurr(~isbadline,:),2,'replicates',100,'emptyaction','singleton');
    rssprev = sum(sum(bsxfun(@minus,xcurr(~isbadline,:),mean(xcurr(~isbadline,:))).^2));
    
    counts = hist(clusteridx1,1:2);
    aiccurr = sum(rsscurr) + 2*nnz(inputs)*2;
    aicprev = rssprev + 2*nnz(inputs);
    
    isclusterchange = aiccurr < aicprev*(1-errtol) && min(counts) >= 5;

    if isclusterchange,

      fprintf('Improvement in AIC by %f%%, clustering\n',(aicprev-aiccurr)/aicprev*100);
      
      clusteridx = ones(1,nlinescurr);
      if any(isbadline),
        xbad = xcurr(isbadline,:);
        idxbadline = find(isbadline);
        err1 = nansum(bsxfun(@minus,xbad,clustercenters(1,:)).^2,2);
        err2 = nansum(bsxfun(@minus,xbad,clustercenters(1,:)).^2,2);
        clusteridx(idxbadline(err1>err2)) = 2;
      end
      clusteridx(~isbadline) = clusteridx1;
      clusterid_perlevel{level+1}(linescurr) = clusterinext_off+clusteridx;
      stats_perlevel{level+1}{clusterinext_off+1} = [statsprev,beststati];
      stats_perlevel{level+1}{clusterinext_off+2} = [statsprev,beststati];
      clusterisdone_perlevel{level+1}(clusterinext_off+1:clusterinext_off+2) = false;
      clusterinext_off = clusterinext_off + 2;
      
    else
      
      fprintf('Insufficient improvement in AIC (by %f%%), not clustering\n',(aicprev-aiccurr)/aicprev*100);
      
      % just copy
      clusterid_perlevel{level+1}(linescurr) = clusterinext_off+1;
      stats_perlevel{level+1}{clusterinext_off+1} = [statsprev,beststati];
      clusterisdone_perlevel{level+1}(clusterinext_off+1) = ~isstatschange;
      clusterinext_off = clusterinext_off + 1;
      
    end

    if clusterisdone_perlevel{level+1}(clusterinext_off),
      fprintf('No change for level %d, cluster %d\n',level,clusteri);
            
    end
    
    figure;
    if numel(statscurr) == 1,
      
      edges = linspace(min(xcurr),max(xcurr),101);
      centers1 = (edges(1:end-1)+edges(2:end))/2;
      
      if isclusterchange,
        c1 = hist(xcurr(clusteridx==1),centers1);
        c2 = hist(xcurr(clusteridx==2),centers1);
        bar(centers1,[c1/sum(c1);c2/sum(c2)]');
      else
        c1 = hist(xcurr,centers1);
        bar(centers1,(c1/sum(c1))');
      end
      xlabel(statfnscurr0{statscurr(1)},'Interpreter','none');      
      title(sprintf('Level %d, cluster %d',level,clusteri));
      
    elseif numel(statscurr) == 2,
      if isclusterchange,
        plot(xcurr(clusteridx==1,1),xcurr(clusteridx==1,2),'k.');
        hold on;
        plot(xcurr(clusteridx==2,1),xcurr(clusteridx==2,2),'r.');
      else
        plot(xcurr(:,1),xcurr(:,2),'k.');
      end
      xlabel(statfnscurr0{statscurr(1)},'Interpreter','none');
      ylabel(statfnscurr0{statscurr(2)},'Interpreter','none');
      title(sprintf('Level %d, cluster %d',level,clusteri));
    else
      [~,proj] = princomp(xcurr);
      if isclusterchange,
        plot(proj(clusteridx==1,1),proj(clusteridx==1,2),'k.');
        hold on;
        plot(proj(clusteridx==2,1),proj(clusteridx==2,2),'r.');
      else
        plot(proj(:,1),proj(:,2),'k.');
      end
      xlabel('PC 1');
      title([sprintf('Level %d, cluster %d',level,clusteri),sprintf(' %s',statfnscurr0{statscurr})],'Interpreter','none');
      ylabel('PC 2');
    end
    
    drawnow;
    
  end
end

%% 

X = nan(nlines,nstatscurr0);
minnsets = 2;
for i = 1:nstatscurr0,
  
  statfn = statfnscurr0{i};
  xcurr = linestats.normmeans.(statfn);
  nsetscurr = linestats.nsets.(statfn);
  
  X(:,i) = xcurr;
  badidx = nsetscurr < minnsets;
  X(badidx,i) = nan;
  
end

% compute weighted distances between lines
Dline = nan(nlines,nlines);
for linei = 1:nlines,
  for linej = 1:nlines,

    if linei == linej,
      Dline(linei,linej) = 0;
    else
      Dline(linei,linej) = nanmean(abs(X(linei,:)-X(linej,:)));
    end
    
  end
end

tmp = Dline;
tmp(eye(size(Dline))==1) = nan;
goodidx = any(~isnan(tmp));
Dline = Dline(goodidx,goodidx);

Dlinevec = squareform(Dline,'tovector');

tree = linkage(Dlinevec,'average');
treecut = cluster(tree,'maxclust',100);

%% 

% [useallims,maxqi,filsizexy,filsizez,savedir,cm,hfig] = ...
%   myparse(varargin,'useallims',false,'maxqi',.1,...
%   'filsizexy',5,'filsizez',2,'savedir','',...
%   'colormap',kjetsmooth(256),'hfig',[]);


savedir = '/nobackup/branson/AverageAnatomyData20130618';
useallims = true;
maxqi = .1;
filsizexy = 5;
filsizez = 2;
cm = kjetsmooth(256);
save /nobackup/branson/AverageAnatomyData20130618/params.mat savedir ...
  useallims maxqi filsizexy filsizez cm anndata;

line_names = linestats.line_names(ismember(linestats.line_names,anndata.line_names));
fid = fopen('linenames_anatomy_20130618.txt','w');
fprintf(fid,'%s\n',line_names{:});
fclose(fid);

%%

datadir = '/nobackup/branson/AverageAnatomyData20130618';
cm = kjetsmooth(256);

statfn = 'fractime_flyfemale_framechase';
staticurr = find(strcmp(statfnscurr0,statfn));
stati = find(strcmp(statfns,statfn));

pvalue = min(fracbigger,fracsmaller)*2;
goodidx = ~isnan(fracbigger_upperbound(:,stati));
q = .01;
[issig,crit_p,tmp] = fdr_bh(fracbigger(goodidx,stati),q,'pdep','yes');
adj_pvalue = nan(nlines,1);
adj_pvalue(goodidx) = tmp;

clf;
plot(linestats.normmeans.(statfn),adj_pvalue,'k.');
set(gca,'YScale','log','XScale','log');
xlabel('Fraction of time chasing, females (difference from control mean)');
ylabel('Adjusted p-value');
box off;

minstds = 5;
linescurr = find(adj_pvalue<=.05 & linestats.normmeans.(statfn)'/controllinestd(stati) >= minstds);
[~,order] = sort(linestats.normmeans.(statfn)(linescurr),'descend');
linescurr = linescurr(order);
line_names = linestats.line_names(linescurr);
savedir = fullfile('/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis',sprintf('%s_anatomy_20130617',statfn));

paramsfile = fullfile('/nobackup/branson/AverageAnatomyData20130618',sprintf('%s_params.mat',statfn));
save(paramsfile,'datadir','cm','line_names','savedir');

%%

cm = kjetsmooth(256);
line_names = {
  'GMR_92B02_AE_01'    'GMR_24E06_AE_01'    'GMR_68C07_AE_01'    'GMR_73E12_AE_01'    'GMR_82D11_AE_01'    'GMR_91B01_AE_01'
  };
datadir = '/nobackup/branson/AverageAnatomyData20130618';
savedir = fullfile('/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis',sprintf('%s_anatomy_20130617','jumpcluster'));
paramsfile = fullfile('/nobackup/branson/AverageAnatomyData20130618',sprintf('%s_params.mat','jumpcluster'));
save(paramsfile,'datadir','cm','line_names','savedir');

%%

line_names = linestats.line_names(ismember(linestats.line_names,anndata.line_names));
savedir = fullfile('/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis',sprintf('alllines_anatomy_20130617'));

paramsfile = fullfile('/nobackup/branson/AverageAnatomyData20130618',sprintf('alllines_params.mat'));
save(paramsfile,'datadir','cm','line_names','savedir');

%% jumps more

statfn = 'fractime_flyany_framejump';

staticurr = find(strcmp(statfnscurr0,statfn));
stati = find(strcmp(statfns,statfn));

goodidx = ~isnan(fracbigger_upperbound(:,stati));
q = .01;
[issig,crit_p,tmp] = fdr_bh(fracbigger(goodidx,stati),q,'pdep','yes');
adj_pvalue = nan(nlines,1);
adj_pvalue(goodidx) = tmp;

clf;
plot(linestats.normmeans.(statfn)+linestats.means.(statfn)(end),adj_pvalue,'k.');
set(gca,'YScale','log','XScale','log');
xlabel('Fraction of time jumping');
ylabel('Adjusted p-value');
box off;

minstds = 5;
linescurr = find(adj_pvalue<=.05 & linestats.normmeans.(statfn)'/controllinestd(stati) >= minstds);
[~,order] = sort(linestats.normmeans.(statfn)(linescurr),'descend');
linescurr = linescurr(order);
line_names = linestats.line_names(linescurr);
savedir = fullfile('/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis',sprintf('%s_anatomy_20130617',statfn));

paramsfile = fullfile('/nobackup/branson/AverageAnatomyData20130618',sprintf('%s_params.mat',statfn));
save(paramsfile,'datadir','cm','line_names','savedir');

meanim = AverageAnatomyImageDeploy(paramsfile);

%%

statfns_short = {
  'fractime_flyany_frameattemptedcopulation'
  'fractime_flyany_framebackup'
  'fractime_flyany_framebodyturn'
  'fractime_flyany_framechase'
  'fractime_flyany_framecopulation'
  'fractime_flyany_framecrabwalkextreme'
  'fractime_flyany_framenotanybehavior'
  'fractime_flyany_framepivotcenter'
  'fractime_flyany_framepivottail'
  'fractime_flyany_framerighting'
  'fractime_flyany_framestop'
  'fractime_flyany_frametouch'
  'fractime_flyany_framewalk'
  'fractime_flyany_framewingextension'
  'fractime_flyany_framewingflick'
  'fractime_flyany_framewinggrooming'
  };

for statii = 1:numel(statfns_short),
  
  statfn = statfns_short{statii};
  
  staticurr = find(strcmp(statfnscurr0,statfn));
  stati = find(strcmp(statfns,statfn));
  
  goodidx = ~isnan(fracbigger_upperbound(:,stati));
  q = .01;
  [issig,crit_p,tmp] = fdr_bh(fracbigger(goodidx,stati),q,'pdep','yes');
  adj_pvalue = nan(nlines,1);
  adj_pvalue(goodidx) = tmp;
  
  figure(statii);
  clf;
  plot(linestats.normmeans.(statfn)+linestats.means.(statfn)(end),adj_pvalue,'k.');
  set(gca,'YScale','log','XScale','log');
  xlabel('Fraction of time jumping');
  ylabel('Adjusted p-value');
  box off;
  
  minstds = 5;
  linescurr = find(adj_pvalue<=.05 & linestats.normmeans.(statfn)'/controllinestd(stati) >= minstds);
  [~,order] = sort(linestats.normmeans.(statfn)(linescurr),'descend');
  linescurr = linescurr(order);
  line_names = linestats.line_names(linescurr);
  savedir = fullfile('/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis',sprintf('%s_anatomy_20130617',statfn));
  
  paramsfile = fullfile('/nobackup/branson/AverageAnatomyData20130618',sprintf('%s_params.mat',statfn));
  save(paramsfile,'datadir','cm','line_names','savedir');
  
  meanim = AverageAnatomyImageDeploy(paramsfile);
  
end


%% choose a random 32 lines

line_names = linestats.line_names(ismember(linestats.line_names,anndata.line_names));
idx = randsample(numel(line_names),32);
line_names = line_names(idx);
savedir = fullfile('/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis',sprintf('random32_anatomy_20130617'));

paramsfile = fullfile('/nobackup/branson/AverageAnatomyData20130618',sprintf('random32_params.mat'));
save(paramsfile,'datadir','cm','line_names','savedir');
meanim = AverageAnatomyImageDeploy(paramsfile);

%% show means for some categories

statfn = 'fractime_flymale_framechase';
staticurr = find(strcmp(statfnscurr0,statfn));
stati = find(strcmp(statfns,statfn));

minstds = 5;
linescurr = find(fracbigger_upperbound(:,stati)<=.01 & linestats.normmeans.(statfn)'/controllinestd(stati) >= minstds);
[~,order] = sort(linestats.normmeans.(statfn)(linescurr),'descend');
linescurr = linescurr(order);
linenames_curr = linestats.line_names(linescurr);

% Lines GMR_32C09_AE_01 GMR_22E02_AE_01  are not imaged

meanim = AverageAnatomyImage(linenames_curr,anndata,'useallims',true,'savedir','chaseanatomy20130617','hfig',1);

save(fullfile('chaseanatomy20130617','meanim.mat'),'meanim','imidx_perline');

% make a z-stack movie
hfig = gcf;
clf;
set(hfig,'Units','pixels','Position',[10,10,1024,512]);
hax = axes('Position',[0,0,1,1]);
z = 1;
him = imagesc(meanim(:,:,1)');
colormap(kjetsmooth(256));
set(gca,'CLim',[0,max(meanim(:))]);
hcb = colorbar('Location','East');
set(hcb,'XColor','w','YColor','w');
axis image off;
htext = text(0,0,sprintf(' z = %d',z),'Color','w','HorizontalAlignment','left',...
  'VerticalAlignment','top','FontSize',24);
set(hfig,'Visible','off');

uavifile = 'chaseanatomy20130617/meanstack.avi';
avifile = 'chaseanatomy20130617/meanstack_xvid.avi';
aviobj = VideoWriter(uavifile,'Uncompressed AVI');
aviobj.FrameRate = 5;
open(aviobj);

fr = getframe_invisible(hax);
[height,width,~] = size(fr);
fprintf('Size of frame is %d x %d\n',height,width);
gfdata = getframe_initialize(hax);
[fr,height,width] = getframe_invisible_nocheck(gfdata,[height,width],false,false);
writeVideo(aviobj,fr);

for z = 2:size(meanim,3),
  set(him,'CData',meanim(:,:,z)');
  set(htext,'String',sprintf(' z = %d',z));
  drawnow;
  fr = getframe_invisible_nocheck(gfdata,[height,width],false);
  writeVideo(aviobj,fr);
end
close(aviobj);

newheight = 4*ceil(height/4);
newwidth = 4*ceil(width/4);
cmd = sprintf('mencoder %s -o %s -ovc xvid -xvidencopts fixed_quant=4 -vf scale=%d:%d -msglevel all=2',...
  uavifile,avifile,newwidth,newheight);
status = system(cmd);
if status ~= 0,
  warning('Could not compress video');
else
  delete(uavifile);
end

%% normalize by mean image over all lines

statfn = 'fractime_flyfemale_framechase';
catdir = fullfile('/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis',sprintf('%s_anatomy_20130617',statfn));
alllinedir = fullfile('/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis',sprintf('alllines_anatomy_20130617'));

catdata = load(fullfile(catdir,'meanim.mat'));
alllinedata = load(fullfile(alllinedir,'meanim.mat'));
idx = ~alllinedata.meanim > 0;
r = randsample(numel(alllinedata.meanim),100);

napprox = 100;
minv = min(min(catdata.meanim(catdata.meanim>0)),min(alllinedata.meanim(alllinedata.meanim>0)));
maxv = max(max(catdata.meanim(:)),max(alllinedata.meanim(:)));
logminv = log10(minv);
logmaxv = log10(maxv);
samplescurr = [0,logspace(logminv,logmaxv,napprox+1)];
[binocdfin1,binocdfin2] = meshgrid(samplescurr,samplescurr);
n = numel(catdata.datafilenames);
binocdfout = 1-binocdf(binocdfin1(:)*n,n,binocdfin2(:));
idx1 = 2+round( (log10(catdata.meanim)-logminv)/(logmaxv-logminv)*napprox );
idx1(catdata.meanim==0) = 1;
idx2 = 2+round( (log10(alllinedata.meanim)-logminv)/(logmaxv-logminv)*napprox );
idx2(alllinedata.meanim==0) = 1;
pvalue_anatomy = binocdfout(sub2ind([napprox+2,napprox+2],idx1,idx2));
pvalue_anatomy(alllinedata.meanim==0) = 1;

xlabel('observed fraction');
ylabel('probability');

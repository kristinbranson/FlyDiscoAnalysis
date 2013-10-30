addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
addpath /groups/branson/bransonlab/projects/olympiad/SAGE/MATLABInterface/Trunk;
addpath /groups/branson/bransonlab/projects/olympiad/anatomy/fileio;

datadir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data';
settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
datalocparamsfilestr = 'dataloc_params.txt';

%% parameters

control_line_names = {'pBDPGAL4U','FCF_pBDPGAL4U_1500437'};
main_control_line_name = 'pBDPGAL4U';
%requiredfiles = {'stats_perframe.mat'};
requiredfiles = {'scores_WingExtension.mat'};

% if the variance is 3x bigger, then it is not worth adding
maxfactorexp = sqrt(3);
minnframes = 200;
minncontrolsets = 10;

%% read in metadata from SAGE

metadata = SAGEListBowlExperiments(...
  'checkflags',false,...
  'removemissingdata',false,...
  'screen_type','primary');

%% read in metadata from file system

dirres = mydir(datadir,'name','.*Bowl._20.*','isdir',true);

metadatadir = [];

for i = 1:numel(dirres),
  
  metadatafile = fullfile(dirres{i},'Metadata.xml');
  [~,name] = fileparts(dirres{i});
  
  fprintf('%d / %d: %s...\n',i,numel(dirres),name);
  
  if ~exist(metadatafile,'file'),
    continue;
  end
  
  try
    metadatacurr = ReadMetadataFile(metadatafile);
  catch ME,
    fprintf('Could not read metadata file %s: %s\n',metadatafile,getReport(ME));
    continue;
  end
  
  if isfield(metadatacurr,'flag_aborted') && metadatacurr.flag_aborted,
    continue;
  end
  
  if isfield(metadatacurr,'flag_redo') && isnumeric(metadatacurr.flag_redo) && metadatacurr.flag_redo,
    continue;
  end
  
  if ~isfield(metadatacurr,'screen_type'),
    continue;
  end
  
  if ~strcmpi(metadatacurr.screen_type,'primary'),
    continue;
  end
  
  metadatacurr.experiment_name = ['FlyBowl_',name];
  metadatacurr.file_system_path = dirres{i};
  
  metadatadir = structappend(metadatadir,metadatacurr);
  fprintf('  %s added\n',name);
  
end

%% compare the two to find what is missing from each

[sageindir,sage2diridx] = ismember({metadata.experiment_name},{metadatadir.experiment_name});
[dirinsage,dir2sageidx] = ismember({metadatadir.experiment_name},{metadata.experiment_name});

tmp = SAGEListBowlExperiments(...
  'checkflags',false,...
  'removemissingdata',false,...
  'experiment_name',{metadatadir(~dirinsage).experiment_name});
ism = ismember({metadatadir(~dirinsage).experiment_name},{tmp.experiment_name});
fprintf('For all experiments found using the directory search and not in SAGE, %d / %d were in SAGE with a different screen_type\n',nnz(ism),numel(ism));
disp(unique({tmp.screen_type}));

isfail = [metadata(~sageindir).flag_aborted] | [metadata(~sageindir).flag_redo];
fprintf('For all experiments found using SAGE and not the directory search, %d / %d have flag aborted or redo set\n',nnz(isfail),numel(isfail));

idx = find(~sageindir);
idx(isfail) = [];

isfile = cellfun(@(x) exist(x,'dir'),{metadata(idx).file_system_path});
fprintf('For the remaining experiments found using SAGE and not the directory search, directories exist for %d / %d \n',nnz(isfile),numel(isfile));

% on 20130911, the remaining 5 directories were those for which the file
% system path and the metadata files did not match. 

% from now on, just use metadata

%% remove weird line names, dates that are not actually our data

isquestionable = ~cellfun(@isempty,regexp({metadata.line_name},'^EXT','once'))' | ...
  datenum({metadata.exp_datetime},'yyyymmddTHHMMSS')>datenum('20120915T000000','yyyymmddTHHMMSS');

% these all are collected on 20130228, and have automated_pf = F

metadata(isquestionable) = [];

isothercontrol = strcmp({metadata.line_name},'FCF_pBDPGAL4U_1500437');
metadata(isothercontrol) = [];

%% remove failures

ismanualf = strcmp({metadata.manual_pf},'F');
fprintf('Removing %d experiments with manual pf\n',nnz(ismanualf));
% Removing 471 experiments with manual pf
metadata(ismanualf) = [];

ignorecategories = {
  'missing_perframestats_files'
};

isautomatedf = ~strcmp({metadata.automated_pf},'P');
fprintf('%d experiments have automated_pf ~= P.\n',nnz(isautomatedf));
% 565 experiments have automated_pf ~= P
idxautomatedpf = find(isautomatedf);

[categories,~,idx] = unique({metadata.automated_pf_category});
counts = hist(idx,1:numel(categories));
for i = 1:numel(categories),
  fprintf('%s: %d\n',categories{i},counts(i));
end

% bad_number_of_flies: 253 -- failures
% bad_number_of_flies_per_sex: 62 -- failures
% flag_aborted_set_to_1: 49 -- failures
% flag_redo_set_to_1: 72 -- failures
% fliesloaded_time_too_long: 20 -- failures
% fliesloaded_time_too_short: 8 -- failures
% incoming_checks_unknown: 2 -- maybe ok, rerun analysis first
% missing_capture_files: 2 -- one of these looks salvageable, one is bad
% missing_extra_diagnostics_files: 4 -- mostly ok
% missing_perframestats_files: 1 -- looks ok
% missing_results_movie_files: 7 -- these also have bad numbers of flies
% missing_video: 2 -- failures
% shiftflytemp_time_too_long: 1 -- ignore
% short_video: 80 -- failures

idx = find(strcmp({metadata.automated_pf_category},'incoming_checks_unknown'));
for i = idx,
  filename = fullfile(metadata(i).file_system_path,'automatic_checks_complete_results.txt');
  if ~exist(filename,'file'),
    fprintf('\n\n***%s: %s does not exist\n',metadata(i).experiment_name,filename);
    continue;
  end
  res = ReadParams(filename);
  fprintf('\n\n***%s: %s, %s\n',metadata(i).experiment_name,res.automated_pf,res.automated_pf_category);
  fprintf(res.notes_curation);
  fprintf('\n');
end

ismaybeok = ismember({metadata.automated_pf_category},ignorecategories);
isautomatedf(ismaybeok) = false;

metadata(isautomatedf) = [];

save CollectedPrimaryMetadata20130912.mat metadata;

%% load in statistics

version = '0.1';
analysis_protocol = '20130909';
mintimestamp = datenum('20130909','yyyymmdd');

nexps = numel(metadata);
expdirs = {metadata.file_system_path};

issuccess = false(1,nexps);
isquestionable = false(1,nexps);
isstale = false(1,nexps);

dataloc_params = ReadParams(fullfile(settingsdir,analysis_protocol,datalocparamsfilestr));
statsparamsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.statsparamsfilestr);
statsperframefeaturesfile = fullfile(settingsdir,analysis_protocol,dataloc_params.statsperframefeaturesfilestr);
stats_perframefeatures = ReadStatsPerFrameFeatures2(statsperframefeaturesfile);

statfns = cell(1,numel(stats_perframefeatures));
for i = 1:numel(stats_perframefeatures),
  
  % which per-frame feature
  field = stats_perframefeatures(i).field;
  % which frames to analyze
  frameconditionname = stats_perframefeatures(i).framecondition;
  % which flies to analyze
  flyconditionname = stats_perframefeatures(i).flycondition;
  
  fn = sprintf('%s_fly%s_frame%s',field,flyconditionname,frameconditionname);
  if numel(fn) > 63,
    fn = fn(1:63);
  end
  statfns{i} = fn;

end
statfns = unique(statfns);
nstats = numel(statfns);

allstatsmat = nan(nstats,nexps);
allstdstatsmat = nan(nstats,nexps);
allnframestotalmat = nan(nstats,nexps);

parfor i = 1:nexps,

  expdir = expdirs{i};
  [~,name] = fileparts(expdir);
  fprintf('%d/%d: %s\n',i,nexps,name);
  
  statsfile = fullfile(expdir,'stats_perframe.mat');
  if ~exist(statsfile,'file'),
    fprintf('File %s does not exist\n',statsfile);
    continue;
  end
  
  statsdata = load(statsfile,'-regexp','(statsperexp)|(version)|(real_analysis_protocol)|(timestamp)');
  
  if ~isfield(statsdata,'version') || ~strcmp(statsdata.version,version),
    isstale(i) = true;
    fprintf('%s: version missing or incorrect.\n',name);
  end
    
  if ~isfield(statsdata,'real_analysis_protocol') || ~strcmp(statsdata.real_analysis_protocol,analysis_protocol),
    isstale(i) = true;
    fprintf('%s: real_analysis_protocol missing or incorrect.\n',name);
  end

  if ~isfield(statsdata,'timestamp') || datenum(statsdata.timestamp,'yyyymmddTHHMMSS') < mintimestamp,
    isstale(i) = true;
    fprintf('%s: timestamp missing or too old\n',name);
  end
  
  
  fnscurr = fieldnames(statsdata.statsperexp);
  newfns = setdiff(fnscurr,statfns);
  missingfns = setdiff(statfns,fnscurr);
  if ~isempty(newfns),
    fprintf('exp %d: %d extra fns:\n',i,numel(newfns));
    fprintf('%s\n',newfns{:});
    isquestionable(i) = true;
  end
  if ~isempty(missingfns),
    fprintf('Experiment %s missing the following stats:\n',name);
    fprintf('%s\n',missingfns{:});
    isquestionable(i) = true;
  end
  
  ism = ismember(statfns,fnscurr);
  
  statsmatcurr = nan(1,nstats),
  stdstatsmatcurr = nan(1,nstats),
  nframestotalmatcurr = nan(1,nstats),
  for j = 1:nstats,
    if ~ism(j),
      continue;
    end
    fn = statfns{j};
    statsmatcurr(j) = statsdata.statsperexp.(fn).meanmean;
    stdstatsmatcurr(j) = statsdata.statsperexp.(fn).stdmean;
    nframestotalmatcurr(j) = statsdata.statsperexp.(fn).Z;
  end
  allstatsmat(:,i) = statsmatcurr;
  allstdstatsmat(:,i) = stdstatsmatcurr;
  allnframestotalmat(:,i) = nframestotalmatcurr;
  issuccess(i) = true;
  
end

%% investigate failures manually

idx = find(~issuccess);

issuccess(idx) = [];
isquestionable(idx) = [];
isstale(idx) = [];
allstatsmat(:,idx) = [];
allstdstatsmat(:,idx) = [];
allnframestotalmat(:,idx) = [];
metadata(idx) = [];
expdirs(idx) = [];
nexps = numel(expdirs);

% no questionable or stale files

allstats = struct;
allstdstats = struct;
allnframestotal = struct;
for j = 1:numel(statfns),
  fn = statfns{j};
  allstats.(fn) = allstatsmat(j,:);
  allstdstats.(fn) = allstdstatsmat(j,:);
  allnframestotal.(fn) = allnframestotalmat(j,:);
end

%% remove experiments without enough frames

MINNFRAMESTOTAL = 400000;
nframestotal = allnframestotal.velmag_ctr_flyany_frameany;
badexps = nframestotal < MINNFRAMESTOTAL;

expdirs(badexps) = [];
metadata(badexps) = [];
nexps = numel(expdirs);
allstats = structfun(@(x) x(~badexps),allstats,'UniformOutput',false);
allstdstats = structfun(@(x) x(~badexps),allstdstats,'UniformOutput',false);
allnframestotal = structfun(@(x) x(~badexps),allnframestotal,'UniformOutput',false);

%% save to file

save CollectedPrimaryPerFrameStats20131024.mat ...
  allnframestotal allstats allstdstats ...
  expdirs main_control_line_name metadata nexps nstats statfns;

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

%% save to file

save -append CollectedPrimaryPerFrameStats20131024.mat ...
  allnframestotal allstats allstdstats control_line_names controlstd ...
  exp2lineidx expdirs idxcontrol ...
  line_names linemeans linemetadata linensets ...
  main_control_line_name metadata nexps nlines nsets nstats ...
  set2lineidx set_names setidx setmeans setmetadata setnexps ...
  statfns;

%% code for comparing old and new versions

olddata = load('CollectedPrimaryPerFrameStats20130905.mat','allstats','allstdstats','allnframestotal','metadata');

[ism,idx] = ismember({olddata.metadata.experiment_name},{metadata.experiment_name});
idxm = find(ism);
fn = 'fractime_flyany_framewingextension';
d = olddata.allstats.(fn)(ism)-allstats.(fn)(idx(ism));
[maxdiff,ii] = max(abs(d));
i = idxm(ii);
fprintf('Biggest change for %s: experiment %s, old = %f, new = %f, difference = %f\n',fn,olddata.metadata(i).experiment_name,...
  olddata.allstats.(fn)(i),...
  allstats.(fn)(idx(i)),...
  olddata.allstats.fractime_flyany_framewinggrooming(i)-allstats.(fn)(idx(i)));


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
  badexp1 = nframestotal < minnframes;
  fprintf('%s: removing %d experiments that don''t have enough frames\n',statfn,nnz(badexp1));
  idxcurr(idxcurr) = ~badexp1;
  
  % count number of experiments per set
  setstats.nexps.(statfn) = hist(setidx(idxcurr),1:numel(set_names));

  % compute statistics
  x = allstats.(statfn)(idxcurr);
  setstats.means.(statfn) = accumarray(setidx(idxcurr)',x',[nsets,1],@mean)';
  setstats.stds.(statfn) = accumarray(setidx(idxcurr)',x',[nsets,1],@(y) std(y,1))';
  setstats.means.(statfn)(setstats.nexps.(statfn)==0) = nan;
  setstats.stds.(statfn)(setstats.nexps.(statfn)==0) = nan;
    
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

save -append CollectedPrimaryPerFrameStats20130912.mat setstats linestats nlines set2lineidx setidx nsets;

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
setstats.usedallcontrols = struct;
setstats.usedalldata = struct;

% for some stats, we should use all control sets or all the data
usealldata = false(1,nstats);
nnotenoughcontrolsets = zeros(1,nstats);
arenocontrolsets = false(1,nstats);
useallcontrolsets = false(1,nstats);
for stati = 1:nstats,
  statfn = statfns{stati};
  idxcurr = setstats.nexps.(statfn)(setiscontrol) >= minnexps;
  ncontrolsetscurr = nnz(idxcurr);
  if ncontrolsetscurr < minncontrolsets,
    usealldata(stati) = true;
    fprintf('Using all data for stat (%d control sets) %s\n',ncontrolsetscurr,statfn);
    continue;
  end

  for seti = 1:nsets,
    
    mindatenumcurr = max(mincontroldatenum,setdatenum(seti)-controldateradius);
    maxdatenumcurr = mindatenumcurr+controldateradius*2+1;
    if maxdatenumcurr > maxcontroldatenum,
      maxdatenumcurr = maxcontroldatenum;
      mindatenumcurr = maxdatenumcurr-controldateradius*2+1;
    end
  
    idxdate = (controldatenum>=mindatenumcurr & controldatenum <= maxdatenumcurr)';
    ncontrolsetscurr = nnz(idxcurr & idxdate);

    if ncontrolsetscurr == 0,
      if ncontrolsetscurr == 0,
        arenocontrolsets(stati) = true;
        break;
      end
      nnotenoughcontrolsets(stati) = nnotenoughcontrolsets(stati)+1;
      %fprintf('Need all control sets for stat %s (%d control sets for range %s to %s)\n',statfn,ncontrolsetscurr,datestr(mindatenumcurr),datestr(maxdatenumcurr));
      %break;
    end
  end
  if nnotenoughcontrolsets(stati) > nsets / 10 || arenocontrolsets(stati),
    if arenocontrolsets(stati),
      fprintf('Need all control sets for stat %s (some have 0 control sets)\n',statfn);
    else
      fprintf('Need all control sets for stat %s (%d sets do not have enough control sets)\n',statfn,nnotenoughcontrolsets(stati));
    end
    useallcontrolsets(stati) = true;
  end
  
end

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
      setstats.controlmeans.(statfn) = nan(1,nsets);
      setstats.nsetscontrolnorm.(statfn) = zeros(1,nsets);
      setstats.usedallcontrols.(statfn) = useallcontrolsets(stati);
      setstats.usedalldata.(statfn) = usealldata(stati);
    end
    
    if usealldata(stati),
      idxcurr = find(setstats.nexps.(statfn) >= minnexps);
    elseif useallcontrolsets(stati),
      idxcurr = setidxcontrol(setstats.nexps.(statfn)(setiscontrol) >= minnexps);
    else
      idxcurr = setidxcontrol(idxdate & setstats.nexps.(statfn)(setiscontrol) >= minnexps);
    end
    idxcurr = setdiff(idxcurr,seti);
      
    setstats.controlmeans.(statfn)(seti) = mean(setstats.means.(statfn)(idxcurr));
    if isnan(setstats.controlmeans.(statfn)(seti)),
      error('set %d %s, stat %d %s\n',seti,setstats.metadata(seti).set,stati,statfn);
    end
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

%% remove lines that do not have enough data

idxremove = linestats.nsets.velmag_ctr_flyany_frameany == 0;

linestatfns = ...
  {'means'
  'stds'
  'nsets'
  'normmeans'
  };

for j = 1:numel(linestatfns),
  linefn = linestatfns{j};
  for i = 1:numel(statfns),
    statfn = statfns{i};
    linestats.(linefn).(statfn)(idxremove) = [];
  end
end

linestats.line_names(idxremove) = [];

linecompfns = {
  'int_manual'
  'dist_manual'
  'int_auto'
  'dist_auto'
  };

for j = 1:numel(linecompfns),
  linefn = linecompfns{j};
  if ~isfield(linestats,linefn),
    continue;
  end
  fns = fieldnames(linestats.(linefn));
  for i = 1:numel(fns),
    fn = fns{i};
    linestats.(linefn).(fn)(idxremove) = [];
  end
end

linemetadata(idxremove) = [];
linemeans(idxremove,:) = [];
linensets(idxremove,:) = [];

line_names = linestats.line_names;
nlines = numel(line_names);
[~,exp2lineidx] = ismember({metadata.line_name},line_names);
[~,set2lineidx] = ismember({setmetadata.line_name},line_names);

%% control means

normcontrolmean = struct;
normcontrolstd = struct;
for fni = 1:nstats,
  statfn = statfns{fni};
  if ~usealldata(fni),
    idxcurr = setstats.nexps.(statfn) >= minnexps & setiscontrol & ...
      ~isnan(setstats.normmeans.(statfn)) & ~isinf(setstats.normmeans.(statfn));
  else
    idxcurr = setstats.nexps.(statfn) >= minnexps & ...
      ~isnan(setstats.normmeans.(statfn)) & ~isinf(setstats.normmeans.(statfn));
  end
  normcontrolmean.(statfn) = mean(setstats.normmeans.(statfn)(idxcurr));
  normcontrolstd.(statfn) = std(setstats.normmeans.(statfn)(idxcurr),1);
end

%% add normmean to allstats

nexps = numel(metadata);
[~,exp2setidx] = ismember({metadata.set},{setstats.metadata.set});
allstatsnorm = struct;
for fni = 1:nstats,
  statfn = statfns{fni};
  x = allstats.(statfn);
  allstatsnorm.(statfn) = allstats.(statfn) - setstats.controlmeans.(statfn)(exp2setidx);
end

%% resave

linestats.metadata = linemetadata;

save -append CollectedPrimaryPerFrameStats20131024.mat setstats linestats nlines set2lineidx nsets linemetadata linemeans linensets line_names exp2lineidx normcontrolmean normcontrolstd nexps allstatsnorm;

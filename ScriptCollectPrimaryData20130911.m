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
minnexps = 2;

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

%% compare automated_pf

automated_pf_fromfile = repmat('?',[1,numel(metadata)]);
for i = 1:numel(metadata),
  
  fprintf('%d / %d: %s\n',i,numel(metadata),metadata(i).experiment_name);
  expdir = metadata(i).file_system_path;
  apffile = fullfile(expdir,'automatic_checks_complete_results.txt');
  if ~exist(apffile,'file'),
    apffile = fullfile(expdir,'automatic_checks_incoming_results.txt');
    if ~exist(apffile,'file'),
      fprintf('No automated checks files exist for %s\n',metadata(i).experiment_name);
      continue;
    end      
  end
  tmp = ReadParams(apffile);
  if ~isfield(tmp,'automated_pf'),
    fprintf('No automated_pf field for %s\n',metadata(i).experiment_name);
    continue;
  end
  automated_pf_fromfile(i) = tmp.automated_pf;
  
end

automated_pf_fromsage = repmat('?',[1,numel(metadata)]);
automated_pf_fromsage(~cellfun(@isempty,{metadata.automated_pf})) = [metadata.automated_pf];
idxmismatch = find(automated_pf_fromfile ~= automated_pf_fromsage);

[tmp1,tmp2,tmp3] = unique([automated_pf_fromfile(idxmismatch)',automated_pf_fromsage(idxmismatch)'],'rows');
for i = 1:size(tmp1,1),
  fprintf('%d / %d experiments have automated_pf = %s on the file system and %s in SAGE:\n',nnz(tmp3==i),numel(metadata),tmp1(i,1),tmp1(i,2));
  fprintf('    %s\n',metadata(idxmismatch(tmp3==i)).experiment_name);
end

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
% 20131122: Removing 604 experiments with manual pf
% 20131127: Removing 648 experiments with manual pf
metadata(ismanualf) = [];

ignorecategories = {
  'missing_perframestats_files'
};

isautomatedf = ~strcmp({metadata.automated_pf},'P');
fprintf('%d experiments have automated_pf ~= P.\n',nnz(isautomatedf));
% 565 experiments have automated_pf ~= P
% 20131122: 553 experiments have automated_pf ~= P.
% 20131127: 553 experiments have automated_pf ~= P.
idxautomatedpf = find(isautomatedf);

[categories,~,idx] = unique({metadata(isautomatedf).automated_pf_category});
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

% 20131122:
% bad_number_of_flies: 254
% bad_number_of_flies_per_sex: 61
% flag_aborted_set_to_1: 49
% flag_redo_set_to_1: 72
% fliesloaded_time_too_long: 20
% fliesloaded_time_too_short: 8
% missing_capture_files: 1
% missing_results_movie_files: 7
% missing_video: 2
% shiftflytemp_time_too_long: 1
% short_video: 78

% 20131127
% bad_number_of_flies: 254
% bad_number_of_flies_per_sex: 61
% flag_aborted_set_to_1: 49
% flag_redo_set_to_1: 72
% fliesloaded_time_too_long: 20
% fliesloaded_time_too_short: 8
% missing_capture_files: 1
% missing_results_movie_files: 7
% missing_video: 2
% shiftflytemp_time_too_long: 1
% short_video: 78

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

save CollectedPrimaryMetadata20131127.mat metadata;

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
fprintf('Removing %d experiments without enough frames\n',nnz(badexps));
% 20131122: Removing 0 experiments without enough frames
% 20131127: Removing 0 experiments without enough frames

expdirs(badexps) = [];
metadata(badexps) = [];
nexps = numel(expdirs);
allstats = structfun(@(x) x(~badexps),allstats,'UniformOutput',false);
allstdstats = structfun(@(x) x(~badexps),allstdstats,'UniformOutput',false);
allnframestotal = structfun(@(x) x(~badexps),allnframestotal,'UniformOutput',false);

%% save to file

save CollectedPrimaryPerFrameStats20131127.mat ...
  allnframestotal allstats allstdstats ...
  expdirs main_control_line_name metadata nexps nstats statfns;

%% remove data when there are not enough samples

for stati = 1:nstats,
  
  statfn = statfns{stati};
  
  % remove experiments that have nans or infs
  % too small or bad nframestotal
  % or not enough 
  nframestotal = allnframestotal.(statfn);  
  badidx = isnan(allstats.(statfn)) | isinf(allstats.(statfn)) | ...
    isinf(nframestotal) | isnan(nframestotal) | ...
    nframestotal < minnframes;
  
  if any(badidx),
    fprintf('Removing %d experiments of data for %s\n',nnz(badidx),statfn);
    allstats.(statfn)(badidx) = nan;
    allstdstats.(statfn)(badidx) = nan;
    allnframestotal.(statfn)(badidx) = nan;
  end
  
end

%% compute set stats using all experiments

[set_names,firstidx,setidx] = unique({metadata.set});
nsets = numel(set_names);
setmeans = nan(nsets,nstats);
setstds = nan(nsets,nstats);
setnexps = nan(nsets,nstats);
setmetadata = metadata(firstidx);
for stati = 1:nstats,
  statfn = statfns{stati};
  
  idxcurr = ~isnan(allstats.(statfn)) & ~isinf(allstats.(statfn));

  setmeans(:,stati) = accumarray(setidx(idxcurr)',allstats.(statfn)(idxcurr)',[nsets,1],@mean,nan);
  setstds(:,stati) = accumarray(setidx(idxcurr)',allstats.(statfn)(idxcurr)',[nsets,1],@(x) std(x,1),nan);
  setnexps(:,stati) = hist(setidx(idxcurr),1:nsets);
end

%% compute line stats using all sets

[ line_names,firstidx,set2lineidx] = unique({setmetadata.line_name});
[~,exp2lineidx] = ismember({metadata.line_name},line_names);
nlines = numel(line_names);
linemeans = nan(nlines,nstats);
linestds = nan(nlines,nstats);
linensets = nan(nlines,nstats);
linenexps = nan(nlines,nstats);
linemetadata = setmetadata(firstidx);

for stati = 1:nstats,
  statfn = statfns{stati};
  idxcurr = ~isnan(setmeans(:,stati)) & ~isinf(setmeans(:,stati)) & ...
    setnexps(:,stati) >= minnexps;
  expidxcurr =  ~isnan(allstats.(statfn)) & ~isinf(allstats.(statfn));

  linemeans(:,stati) = accumarray(set2lineidx(idxcurr)',setmeans(idxcurr,stati)',[nlines,1],@mean,nan);

  % standard deviation over experiments!
  linestds(:,stati) = accumarray(exp2lineidx(expidxcurr)',allstats.(statfn)(expidxcurr)',[nlines,1],@(x) std(x,1),nan);

  linensets(:,stati) = hist(set2lineidx(idxcurr),1:nlines);
  linenexps(:,stati) = hist(exp2lineidx(expidxcurr),1:nlines);
end

idxcontrol = find(strcmp(line_names,main_control_line_name));
controlstd = nanstd(setmeans(set2lineidx==idxcontrol,:),1,1);
controlmean = nanmean(setmeans(set2lineidx==idxcontrol,:),1);

%% save to file

save -append CollectedPrimaryPerFrameStats20131127.mat ...
  allnframestotal allstats allstdstats control_line_names controlstd ...
  exp2lineidx expdirs idxcontrol ...
  line_names linemeans linemetadata linensets linestds linenexps ...
  main_control_line_name metadata nexps nlines nsets nstats ...
  set2lineidx set_names setidx setmeans setstds setmetadata setnexps ...
  statfns;

%% code for comparing old and new versions

olddata = load('CollectedPrimaryPerFrameStats20131122.mat','allstats','allstdstats','allnframestotal','metadata','setstats','linestats');

[ism,idx] = ismember({olddata.metadata.experiment_name},{metadata.experiment_name});
idxm = find(ism);
for stati = 1:numel(statfns),
  fn = statfns{stati};
  %fn = 'fractime_flyany_framewingextension';
  d = olddata.allstats.(fn)(ism)-allstats.(fn)(idx(ism));
  [maxdiff,ii] = max(abs(d));
  if maxdiff > 0,
    i = idxm(ii);
    fprintf('Biggest change for %s: experiment %s, old = %f, new = %f, difference = %f\n',fn,olddata.metadata(i).experiment_name,...
      olddata.allstats.(fn)(i),...
      allstats.(fn)(idx(i)),...
      olddata.allstats.(fn)(i)-allstats.(fn)(idx(i)));
  end
end

%% collect setstats

[set_names,firstidx,setidx] = unique({metadata.set});
nsets = numel(set_names);

setstats = struct('means',struct,'stds',struct,'nexps',struct,'metadata',struct);
setstats.metadata = metadata(firstidx);

for stati = 1:numel(statfns),

  statfn = statfns{stati};
  
  setstats.means.(statfn) = setmeans(:,stati)';
  setstats.stds.(statfn) = setstds(:,stati)';
  setstats.nexps.(statfn) = setnexps(:,stati)';

end

%% collect line stats

linestats = struct('means',struct,'stds',struct,'nsets',struct,'metadata',struct);
[linestats.line_names,firstidx,set2lineidx] = unique({setstats.metadata.line});
linestats.metadata = setstats.metadata(firstidx);
nlines = numel(linestats.line_names);
for stati = 1:numel(statfns),
  
  statfn = statfns{stati};
  
  linestats.means.(statfn) = linemeans(:,stati)';
  linestats.stds.(statfn) = linestds(:,stati)';
  linestats.nsets.(statfn) = linensets(:,stati)';
  linestats.nexps.(statfn) = linenexps(:,stati)';

end

%% save

save -append CollectedPrimaryPerFrameStats20131127.mat setstats linestats nlines set2lineidx setidx nsets;

%% compare old and new

[ism,idx] = ismember({olddata.linestats.metadata.line_name},{linestats.metadata.line_name});
lines_changed = {};
for stati = 1:numel(statfns),
  fn = statfns{stati};
  %fn = 'fractime_flyany_framewingextension';
  d = olddata.linestats.means.(fn)(ism)-linestats.means.(fn)(idx(ism));
  ischange = ~isnan(d) & d ~= 0;
  if any(ischange),
    fprintf('%d lines changed for %s\n',nnz(ischange),fn);
    lines_changed = union(lines_changed,{linestats.metadata(ischange).line_name});
  end
%   [maxdiff,ii] = max(abs(d));
%   if maxdiff > 0,
%     i = idxm(ii);
%     fprintf('Biggest change for %s: line %s, old = %f, new = %f, difference = %f, %f z\n',fn,olddata.linestats.metadata(i).line_name,...
%       olddata.linestats.means.(fn)(i),...
%       linestats.means.(fn)(idx(i)),...
%       (olddata.linestats.means.(fn)(i)-linestats.means.(fn)(idx(i))),...
%       (olddata.linestats.means.(fn)(i)-linestats.means.(fn)(idx(i)))/controlstd(stati));
%   end
end
fprintf('%d lines changed for some feature:\n',numel(lines_changed));
fprintf('%s\n',lines_changed{:});

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

% compute the global control mean and standard deviation
controlmean = nanmean(setmeans(set2lineidx==idxcontrol,:),1);
controlstd = nanstd(setmeans(set2lineidx==idxcontrol,:),1,1);
controlmean(usealldata) = nanmean(setmeans(:,usealldata),1);
controlstd(usealldata) = nanstd(setmeans(:,usealldata),1,1);

for stati = 1:nstats,
  
  statfn = statfns{stati};
  
  setstats.normmeans.(statfn) = nan(1,nsets);
  setstats.controlmeans.(statfn) = nan(1,nsets);
  setstats.nsetscontrolnorm.(statfn) = zeros(1,nsets);
  setstats.usedallcontrols.(statfn) = useallcontrolsets(stati);
  setstats.usedalldata.(statfn) = usealldata(stati);
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
    setstats.normmeans.(statfn)(seti) = setstats.means.(statfn)(seti) - setstats.controlmeans.(statfn)(seti) + controlmean(stati);

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

linestats.usedallcontrols = setstats.usedallcontrols;
linestats.usedalldata = setstats.usedalldata;

%% normalized exp stats

expstats = struct;
expstats.means = struct;
expstats.stds = struct;
expstats.nframes = struct;
expstats.normmeans = struct;
expstats.metadata = metadata;

for i = 1:numel(statfns),
  fn = statfns{i};
  expstats.means.(fn) = allstats.(fn);
  expstats.stds.(fn) = allstdstats.(fn);
  expstats.nframes.(fn) = allnframestotal.(fn);
end

[~,exp2setidx] = ismember({expstats.metadata.set},{setstats.metadata.set});
setiscontrol = strcmp({setstats.metadata.line_name},main_control_line_name);

for i = 1:numel(statfns),
  fn = statfns{i};
  
  expstats.normmeans.(fn) = expstats.means.(fn) - setstats.controlmeans.(fn)(exp2setidx) + controlmean(i);

end

%% normalized standard deviations over sets

for stati = 1:nstats,
  statfn = statfns{stati};
  expidxcurr =  ~isnan(expstats.normmeans.(statfn)) & ~isinf(expstats.normmeans.(statfn));

  % standard deviation over experiments!
  linestats.normstds.(statfn) = ...
    accumarray(exp2lineidx(expidxcurr)',expstats.normmeans.(statfn)(expidxcurr)',[nlines,1],@(x) std(x,1),nan)';
end


%% remove lines that do not have enough data

idxremove = linestats.nsets.velmag_ctr_flyany_frameany == 0;

linestatfns = ...
  {'means'
  'stds'
  'nsets'
  'nexps'
  'normmeans'
  'normstds'
  };

for j = 1:numel(linestatfns),
  linefn = linestatfns{j};
  for i = 1:numel(statfns),
    statfn = statfns{i};
    linestats.(linefn).(statfn)(idxremove) = [];
  end
end

linestats.line_names(idxremove) = [];
linestats.metadata(idxremove) = [];

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
linenexps(idxremove,:) = [];
linestds(idxremove,:) = [];

line_names = linestats.line_names;
nlines = numel(line_names);
[~,exp2lineidx] = ismember({metadata.line_name},line_names);
[~,set2lineidx] = ismember({setmetadata.line_name},line_names);

%% control means
% 
% normcontrolmean = struct;
% normcontrolstd = struct;
% for fni = 1:nstats,
%   statfn = statfns{fni};
%   if ~usealldata(fni),
%     idxcurr = setstats.nexps.(statfn) >= minnexps & setiscontrol & ...
%       ~isnan(setstats.normmeans.(statfn)) & ~isinf(setstats.normmeans.(statfn));
%   else
%     idxcurr = setstats.nexps.(statfn) >= minnexps & ...
%       ~isnan(setstats.normmeans.(statfn)) & ~isinf(setstats.normmeans.(statfn));
%   end
%   normcontrolmean.(statfn) = mean(setstats.normmeans.(statfn)(idxcurr));
%   normcontrolstd.(statfn) = std(setstats.normmeans.(statfn)(idxcurr),1);
% end
% 
% %% add normmean to allstats
% 
% nexps = numel(metadata);
% [~,exp2setidx] = ismember({metadata.set},{setstats.metadata.set});
% allstatsnorm = struct;
% for fni = 1:nstats,
%   statfn = statfns{fni};
%   x = allstats.(statfn);
%   allstatsnorm.(statfn) = allstats.(statfn) - setstats.controlmeans.(statfn)(exp2setidx);
% end

%% resave

save -append CollectedPrimaryPerFrameStats20131127.mat setstats linestats nlines set2lineidx nsets linemetadata linemeans linensets line_names exp2lineidx nexps expstats;

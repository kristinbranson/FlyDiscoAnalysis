function output = FlyBowlComputeFracTimePerSetAndLine(metadata,varargin)

output = struct;

%% parse parameters

[analysis_protocol,settingsdir,datalocparamsfilestr,classifierparamsfiles] = ...
  myparse(varargin,'analysis_protocol','current',...
  'settingsdir','/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings',...
  'datalocparamsfilestr','dataloc_params.txt',...
  'classifierparamsfiles',{});


%% get list of scores files

if isempty(classifierparamsfiles),
  datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
  dataloc_params = ReadParams(datalocparamsfile);
  paramsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.jaabadetectparamsfilestr);
  params = ReadParams(paramsfile);
  classifierparamsfiles = cellfun(@(x) fullfile(settingsdir,analysis_protocol,x),params.classifierparamsfiles,'UniformOutput',false);
end
classifierparams = ReadClassifierParamsFile(classifierparamsfiles,'isrelativepath',true);
scoresfilestr = cellfun(@(x) x.scorefilename,{classifierparams.file},'UniformOutput',false);
nbehaviors = numel(scoresfilestr);

behaviornames = cell(1,nbehaviors);
for i = 1:nbehaviors,
  if iscell(classifierparams(i).behaviors.names),
    behaviornames{i} = classifierparams(i).behaviors.names{1};
  else
    behaviornames{i} = classifierparams(i).behaviors.names;
  end
end
output.scoresfilestrs = scoresfilestr;
output.behaviornames = behaviornames;

%% load in scores

nexps = numel(metadata);
fractime = nan(nexps,nbehaviors);
isdata = true(nexps,nbehaviors);

for i = 1:nexps,
  if mod(i,10) == 0,
    fprintf('%d / %d\n',i,nexps);
  end
  expdir = metadata(i).file_system_path;
  parfor j = 1:nbehaviors,
    scoresfile = fullfile(expdir,scoresfilestr{j});
    if ~exist(scoresfile,'file'),
      isdata(i,j) = false;
      continue;
    end
    scoresdata = load(scoresfile);
    npos = 0;
    n = sum(scoresdata.allScores.tEnd-scoresdata.allScores.tStart+1);
    for fly = 1:numel(scoresdata.allScores.scores),
      npos = npos + sum(scoresdata.allScores.t1s{fly}-scoresdata.allScores.t0s{fly});
    end
    fractime(i,j) = npos/n;
  end
  if any(~isdata(i,:)),
    [~,name] = fileparts(expdir);
    fprintf('%s: Missing the following scores files:\n',name);
    fprintf('  %s\n',scoresfilestr{~isdata(i,:)});
  end
end

output.isdata = isdata;
output.expmetadata = metadata;
output.experiment_names = {metadata.experiment_name};
output.mean_fractime_perexp = fractime;
output.linename_perexp = {metadata.line_name};

%% compute mean and std of scores per set

[sets,~,setidx] = unique({metadata.set});
nsets = numel(sets);
mean_fractime_perset = nan(nsets,nbehaviors);
std_fractime_perset = nan(nsets,nbehaviors);
for i = 1:nsets,
  idx = setidx == i;
  mean_fractime_perset(i,:) = nanmean(fractime(idx,:),1);
  std_fractime_perset(i,:) = nanstd(fractime(idx,:),1,1);
end

setmetadata = cell(1,nsets);
linename_perset = cell(1,nsets);
for i = 1:nsets,
  setmetadata{i} = metadata(setidx==i);
  linename_perset{i} = setmetadata{i}(1).line_name;
end

nexpsperset = hist(setidx,1:nsets);

output.setmetadata = setmetadata;
output.set_names = sets;
output.mean_fractime_perset = mean_fractime_perset;
output.std_fractime_perset = std_fractime_perset;
output.linename_perset = linename_perset;
output.nexpsperset = nexpsperset;

%% per-line data

minexpsperset = 1;
enoughexps = nexpsperset >= minexpsperset;

line_names = unique(linename_perset(enoughexps));
[~,lineidx] = ismember(linename_perset,line_names);
nlines = numel(line_names);

mean_fractime_perline = nan(nlines,nbehaviors);
std_fractime_perline = nan(nsets,nbehaviors);
for i = 1:nlines
  idx = lineidx == i & enoughexps;
  mean_fractime_perline(i,:) = nanmean(mean_fractime_perset(idx,:),1);
  std_fractime_perline(i,:) = nanstd(mean_fractime_perset(idx,:),1,1);
end
nsetsperline = hist(lineidx(enoughexps),1:nlines);

output.line_names = line_names;
output.mean_fractime_perline = mean_fractime_perline;
output.std_fractime_perline = std_fractime_perline;
output.nsetsperline = nsetsperline;
output.mean_fractime_perline = mean_fractime_perline;
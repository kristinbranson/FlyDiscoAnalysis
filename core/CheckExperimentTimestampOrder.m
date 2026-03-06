function [isbroken,filesmissing,filetypes,rules] = CheckExperimentTimestampOrder(expdir,varargin)

[DEBUG,outputdir,settingsdir,analysis_protocol,datalocparamsfilestr,MAXDT,allperframefns] = ...
  myparse(varargin,'debug',false,...
  'outputdir','/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/TimestampCheckOut',...
  'settinsgdir','/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings',...
  'analysis_protocol','current',...
  'datalocparamsfilestr','dataloc_params.txt',...
  'maxdt',1/60/24,...
  'allperframefns',0);

[~,experiment_name] = fileparts(expdir);


% read in perframe features
if ~iscell(allperframefns),
  datalocparams = ReadParams(fullfile(settingsdir,analysis_protocol,datalocparamsfilestr));
  perframefnsfile = fullfile(settingsdir,analysis_protocol,datalocparams.perframefnsfilestr);
  allperframefns = importdata(perframefnsfile);
end
ignoreperframefns = {'x','y','a','b','theta','x_mm','y_mm','a_mm','b_mm','theta_mm','timestamps','dt','sex'};
perframefns = setdiff(allperframefns,ignoreperframefns);

files = struct;
files.datacap = {'movie.ufmf','SUCCESS','QuickStats.png','QuickStats.txt',...
  'temperature.txt','ufmf_diagnostics.txt','ufmf_log.txt'};
files.autoincoming = {'automatic_checks_incoming_results.txt'};
files.ctrax = {'ctrax_diagnostics.txt','ctrax_results.mat','movie.ufmf.ann'};
files.reg = {'registrationdata.mat','registrationdata.txt','registrationimage.png'};
files.sex = {'sexclassifier.mat','sexclassifier_diagnostics.txt'};
files.trx = {'registered_trx.mat'};
files.perframe = cellfun(@(x) sprintf('perframe/%s.mat',x),perframefns,'UniformOutput',false);
files.stats = {'stats_perframe.txt','stats_perframe.mat','hist_perframe.txt','hist_perframe.mat'};
%files.plots = cellfun(@(x) sprintf('analysis_plots/hist_%s.png',x),allperframefns,'UniformOutput',false);
%files.plots{end+1} = 'stats.png';
files.plots = {'analysis_plots/stats.png'};
files.diagnostics = {'temperature_diagnostics.txt',...
  'bkgd_diagnostics.txt','bkgd_diagnostics.mat','bkgd_diagnostics.png',...
  'bias_diagnostics.txt','bias_diagnostics.mat','bias_diagnostics.png',...
  'video_diagnostics.txt','video_diagnostics.mat','video_diagnostics.png'};
files.resultsmovie = {sprintf('ctrax_results_movie_%s.avi',experiment_name)};
files.autocomplete = {'automatic_checks_complete_results.txt'};

rules = {
  'datacap','autoincoming'
  'autoincoming','ctrax'
  'ctrax','reg'
  'reg','sex'
  'reg','trx'
  'sex','perframe'
  'perframe','stats'
  'stats','plots'
  'ctrax','diagnostics'
  'sex','resultsmovie'
  'resultsmovie','autocomplete'
  'plots','autocomplete'
  };

filetypes = fieldnames(files);
nfiletypes = numel(filetypes);

[~,typeidx] = ismember(rules,filetypes);

% allocate
filesmissing = {};
maxtimestamp = -inf(1,nfiletypes);
mintimestamp = inf(1,nfiletypes);

for i = 1:numel(filetypes),
  filescurr = files.(filetypes{i});
  for j = 1:numel(filescurr),
    tmp = dir(fullfile(expdir,filescurr{j}));
    if isempty(tmp),
      filesmissing{end+1} = filescurr{j}; %#ok<AGROW>
      continue;
    end
    maxtimestamp(i) = max(maxtimestamp(i),tmp.datenum);
    mintimestamp(i) = min(mintimestamp(i),tmp.datenum);
  end
end

nrules = size(rules,1);
isbroken = false(1,nrules);
for i = 1:nrules,
  file1 = typeidx(i,1);
  file2 = typeidx(i,2);
  if isinf(maxtimestamp(file1)) || isinf(mintimestamp(file2)),
    continue;
  end
  isbroken(i) = maxtimestamp(file1) > mintimestamp(file2) + MAXDT;
end

if DEBUG,
  fid = 1;
else
  if isempty(outputdir),
    return;
  end
  if ~exist(outputdir,'dir'),
    mkdir(outputdir);
  end
  outfilename = fullfile(outputdir,sprintf('%s_CheckTimestamps.txt',experiment_name));
  fid = fopen(outfilename,'w');
end

fprintf(fid,'filesmissing');
fprintf(fid,',%s',filesmissing{:});
fprintf(fid,'\n');
fprintf(fid,'brokenrules');
if ~any(isbroken),
  fprintf(',');
end
for i = find(isbroken),
  fprintf(',%s<%s',rules{i,1},rules{i,2});
end
fprintf(fid,'\n');

if fid > 1,
  fclose(fid);
end
% ScriptCheckExperiments20120726

addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
addpath /groups/branson/bransonlab/projects/olympiad/SAGE/MATLABInterface/Trunk;
addpath /groups/branson/bransonlab/projects/olympiad/FlyBowlDataCapture;
addpath /groups/branson/home/bransonk/behavioranalysis/code/Jdetect/Jdetect/experiments/starvation;
addpath /groups/branson/home/bransonk/behavioranalysis/code/Jdetect/Jdetect/perframe;
addpath /groups/branson/home/bransonk/behavioranalysis/code/Jdetect/Jdetect/experiments/age;
settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
if true,
  rootdatadirs = {'/groups/branson/bransonlab/projects/JAABA/data/FlyBowl'
    '/groups/branson/bransonlab/projects/JAABA/data/groundtruth_pBDPGAL4U_data'};
else
  rootdatadir0 = '/groups/branson/home/robiea/Projects_data/JAABA';
  tmp = dir(fullfile(rootdatadir0,'Data*'));
  rootdatadirs = cellfun(@(x) fullfile(rootdatadir0,x),{tmp.name},'UniformOutput',false);
  rootdatadirs{end+1} = '/groups/branson/home/robiea/Projects_data/JAABA/GroundTruth';
  rootdatadirs{end+1} = '/groups/branson/home/robiea/Projects_data/JAABA/GroundTruth/specific_behaviors';  
  rootdatadirs{end+1} = '/groups/branson/home/robiea/Projects_data/JAABA/test_movies';  
end
outexpdirfilename = '/groups/branson/home/bransonk/behavioranalysis/code/Jdetect/Jdetect/experiments/datamanagement/expdirs_alicejaaba20120910.txt';
analysis_protocol = '20120706';
mintimestamp = '20120501T000000';
npar = 8;


%% set data locations

expdirs = {};
clear metadata;
for i = 1:numel(rootdatadirs),
  rootdatadir = rootdatadirs{i};
  tmp = dir(fullfile(rootdatadir,'*20*'));
  tmp = tmp([tmp.isdir]);
  tmp = cellfun(@(x) fullfile(rootdatadir,x),{tmp.name},'UniformOutput',false);
  isdata = false(size(tmp));
  for j = 1:numel(tmp),
    [metadatacurr,success] = parseExpDir(tmp{j});
    if ~success, 
      fprintf('Skipping %s\n',tmp{j});
      continue; 
    end
    metadatacurr.file_system_path = tmp{j};
    metadatacurr.exp_datetime = metadatacurr.date;
    metadatacurr = rmfield(metadatacurr,'date');
    isdata(j) = true;
    if i == 1 && j == 1,
      metadata = metadatacurr;
    else
      metadata(end+1) = metadatacurr;
    end
  end
  fprintf('Found %d directories in %s\n',nnz(isdata),rootdatadir);
  expdirs = [expdirs,tmp(isdata)];
end

% fid = fopen(outexpdirfilename,'w');
% fprintf(fid,'%s\n',expdirs{:});
% fclose(fid);

for i = 1:numel(expdirs),
  [metadatacurr,success] = parseExpDir(expdirs{i});
  metadatacurr.file_system_path = expdirs{i};
  metadatacurr.exp_datetime = metadatacurr.date;
  metadatacurr = rmfield(metadatacurr,'date');
  if i == 1,
    metadata = metadatacurr;
  else
    metadata(end+1) = metadatacurr;
  end
end

%% initialize

issuccess = nan(1,numel(expdirs));
msgs = cell(1,numel(expdirs));
nproblems = 0;
ii = 0;
datenums = datenum({metadata.exp_datetime},'yyyymmddTHHMMSS');

hfig = 1;
figure(hfig);
clf;

hdata = plot(datenums,issuccess,'k.-');
set(gca,'XLim',[min(datenums)-1,max(datenums)+1],'YLim',[-.25,1.25]);
datetick('x','mmmyy','keeplimits');
hti = title(sprintf('%d problems / %d checked',nproblems,ii));

order = randperm(numel(metadata));

%% check

for ii = 1:npar:numel(metadata),
  
  %i = order(ii);
  
  idxcurr = order(ii:min(ii+npar-1,numel(metadata)));
  currissuccess = nan(1,numel(idxcurr));
  currmsgs = cell(1,numel(idxcurr));
  curr_expdirs = {metadata(idxcurr).file_system_path};
  
  parfor j = 1:numel(idxcurr),
    fprintf('***** Checking experiment %s...\n',curr_expdirs{j});
    [currissuccess(j),currmsgs{j}] = CheckExperiment20120726(curr_expdirs{j},'analysis_protocol',analysis_protocol);
  end
  
  issuccess(idxcurr) = currissuccess;
  msgs(idxcurr) = currmsgs;
  nproblems = nnz(issuccess==0);
  
  set(hdata,'YData',issuccess);
  set(hti,'String',sprintf('%d problems / %d checked',nproblems,ii+numel(idxcurr)-1));
  
  drawnow;
  
  if mod(ii,100) < npar,
    save(sprintf('CheckExperiments_robiea_%s_20120909.mat',datatype),'issuccess','msgs','order','ii','metadata');
  end
  
end

for i = find(~issuccess),
  fprintf('\n%s\n',metadata(i).file_system_path);
  fprintf('%s\n',msgs{i}{:});
end

pff_update = {'du_tail','dv_tail',...
  'magveldiff_nose2ell','magveldiff_anglesub',...
  'veltoward_nose2ell','veltoward_anglesub',...
  'closestfly_nose2ell_angle_min30to30','dnose2ell_nose2ell_angle_min30to30',...
  'closestfly_nose2ell_angle_min20to20','dnose2ell_nose2ell_angle_min20to20',...
  'closestfly_nose2ell_angle_30tomin30','dnose2ell_nose2ell_angle_30tomin30',...
  'dt'};

for i = find(~issuccess),
  for j = 1:numel(pff_update),
    tmp = any(~cellfun(@isempty,regexp(msgs{i},pff_update{j},'once')));
    if tmp,
      fprintf('%d %s: %s\n',i,metadata(i).file_system_path,pff_update{j});
    end
  end
end
    


%% check for scores files

startexp = 1;
endexp = 20001;
classifierparamsfiles = {'../JAABA_classifier_params1.txt','../JAABA_classifier_params2.txt','../JAABA_classifier_params3.txt'};

classifierfiles = {};
configfiles = {};
for k = 1:numel(classifierparamsfiles),
classifierparamsfile = classifierparamsfiles{k};

fid = fopen(classifierparamsfile,'r');
while true,
  l = fgetl(fid);
  if ~ischar(l),
    break;
  end
  if isempty(l),
    continue;
  end
  ws = regexp(l,',','split');
  classifierfiles{end+1} = ws{1};  %#ok<SAGROW>
  configfiles{end+1} = ws{2};  %#ok<SAGROW>
end
fclose(fid);
end

scorefiles = cell(1,numel(configfiles));
for i = 1:numel(configfiles),
  
  % read parameters
  fprintf('Reading classifier parameters...\n');
  configparams = ReadXMLParams(configfiles{i});
  
  
  if ~isfield(configparams.file,'scorefilename'),
    if ~isfield(configparams,'behaviors'),
      error('configparams %d does not have field behaviors',i);
    elseif ~isfield(configparams.behaviors,'names'),
      error('configparams{%d}.behaviors does not have field names',i);
    else
      scorefilename = ['scores_',configparams.behaviors.names,'.mat'];
    end
  else
    scorefilename = configparams.file.scorefilename;
  end
  
  scorefiles{i} = scorefilename;
  
end

isok = nan(1,numel(metadata));
for i = startexp:endexp,
  if i > numel(metadata), break; end
  isscore = false(1,numel(scorefiles));
  mindatenum = inf;
  for j = 1:numel(scorefiles),
    tmp = dir(fullfile(metadata(i).file_system_path,scorefiles{j}));
    isscore(j) = ~isempty(tmp);
    if ~isempty(tmp),
      mindatenum = min(mindatenum,tmp.datenum);
    end
  end
  isok(i) = all(isscore);
  [~,name] = fileparts(metadata(i).file_system_path);
  if any(~isscore),
    fprintf('%s missing',name);
    fprintf(' %s',scorefiles{~isscore});
    fprintf('\n');
  else
    fprintf('%s: %s\n',name,datestr(mindatenum));
  end
end

%% check for wing tracking results

% fid = fopen('expdirs_wingtrack_ar20120910.txt','w');
% fprintf(fid,'%s\n',metadata.file_system_path);
% fclose(fid);
for i = 1:numel(metadata),
  filename = fullfile(metadata(i).file_system_path,'wingtracking_results.mat');
  tmp = dir(filename);
  if isempty(tmp),
    fprintf('%s: MISSING\n',filename);
  else
    fprintf('%s: %s\n',filename,tmp.date);
  end
end
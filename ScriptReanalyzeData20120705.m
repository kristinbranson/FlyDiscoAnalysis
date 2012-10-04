%% experiments to reanalyze in various ways

%% set up path

addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
addpath /groups/branson/bransonlab/projects/olympiad/SAGE/MATLABInterface/Trunk;

settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
rootdir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data';

%% get list of all experiments

data = SAGEListBowlExperiments('removemissingdata',false,'checkflags',false);
for i = 1:numel(data),
  if mod(i,100) == 0,
    fprintf('%d / %d\n',i,numel(data));
  end
  filename = fullfile(data(i).file_system_path,'analysis_protocol.txt');
  didread = false;
  if exist(filename,'file'),
    try
      fid = fopen(filename,'r');
      data(i).full_analysis_protocol = fgetl(fid);
      didread = true;
      fclose(fid);
    catch ME
      error(getReport(ME));
    end
  end
  if ~didread,
    data(i).full_analysis_protocol = data(i).analysis_protocol;
  end
end

fns = {'ctrax_version','ctrax_settings','analysis_settings','analysis_version'};
for i = 1:numel(data),  
  s = regexp(data(i).full_analysis_protocol,',','split');
  ss = regexp(s,':','split');
  names = cellfun(@(x)x{1},ss,'UniformOutput',false);
  [ism,idx] = ismember(fns,names);
  for k = 1:numel(fns),
    if ism(k),
      data(i).(fns{k}) = ss{idx(k)}{2};
    else
      data(i).(fns{k}) = 'unknown';
    end
  end
end

data0 = data;

save AllExperimentMetadata20120705.mat data0;

isprimary = strcmp({data0.screen_type},'primary');

%% list of failure experiments that are not salvageable

allowed_automated_pf_categories = {
  'bad_number_of_flies'
  'bad_number_of_flies_per_sex'
  'incoming_checks_unknown'
  'missing_automated_checks_incoming_files'
  'missing_extra_diagnostics_files'
  'missing_metadata_fields'
  'missing_other_analysis_files'
  'missing_perframestats_files'
  'missing_registration_files'
  'missing_results_movie_files'
  'missing_sexclassification_files'
  'missing_tracking_files'
};
isfailure = ( ([data0.automated_pf] == 'F') & ...
  ~ismember({data0.automated_pf_category},allowed_automated_pf_categories) ) | ...
  [data0.manual_pf] == 'F';

%% get list of all archived experiments

isarchived = [data0.archived] == '1';
fid = fopen('expdirs_movefromarchive_20120710.txt','w');
fprintf(fid,'%s\n',data0(isarchived&~isfailure&isprimary).file_system_path);
fclose(fid);

%% experiments that need to be retracked

retrack_min_exp_datetime = '20110601T000000';
retrack_max_ctrax_settings = '20111221';

retrack_max_ctrax_version = '0.2';
retrack_max_maxarea = 180;

retrack_idx = lexcmp({data0.ctrax_settings},retrack_max_ctrax_settings) < 0 & ...
  lexcmp({data0.exp_datetime},retrack_min_exp_datetime) > 0 & ...
  isprimary & ~isfailure;

fid = fopen('expdirs_retrack_20120705.txt','w');
fprintf(fid,'%s\n',data0(retrack_idx).file_system_path);
fclose(fid);

% spot check annfiles to make sure these are indeed old
data = data0(retrack_idx);
ncheck = 50;
[~,order] = sort({data.exp_datetime});
idxcheck = order(unique(round(linspace(1,numel(data),ncheck))));

for ii = 1:ncheck,
  i = idxcheck(ii);
  annfile = fullfile(data(i).file_system_path,'movie.ufmf.ann');
  if ~exist(annfile,'file'),
    warning('Annotation file %s does not exist',annfile);
    continue;
  end
  [maxarea,version] = read_ann(annfile,'maxarea','version');
  if maxarea >= retrack_max_maxarea || ...
      lexcmp(version,retrack_max_ctrax_version) >= 0,
    warning('Area and version check failed for %s\n',data(i).experiment_name);
  end
end

%% make list of experiments that are archived and need to be retracked

idxmove = isarchived & retrack_idx;
fid = fopen('expdirs_movemoviefromarchive_20120710.txt','w');
fprintf(fid,'%s\n',data0(idxmove).file_system_path);
fclose(fid);

%% make list of experiments that need to have per-frame features fixed:
% fixing includes fixing the bug in storing trajectory features and
% changing how dt is computed to use just the median dt instead of the
% per-frame timestamps. 

idxfixpff = ~isfailure & isprimary & ~retrack_idx;
fid = fopen('expdirs_fixperframefeatures_20120710.txt','w');
fprintf(fid,'%s\n',data0(idxfixpff).file_system_path);
fclose(fid);


%% check for missing per-frame fns

perframefnsfile = fullfile('settings','current','perframefns.txt');
perframefns = importdata(perframefnsfile);
idxcheck = ~retrack_idx & isprimary & ~isfailure;
data = data0(idxcheck);
ismissing = cell(1,numel(data));
parfor i = 1:numel(data),
  if mod(i,100) == 0,
    fprintf('%d/%d\n',i,numel(data));
  end
  ismissing{i} = false(1,numel(perframefns));
  perframedir = fullfile(data(i).file_system_path,'perframe');
  for j = 1:numel(perframefns),
    ismissing{i}(j) = ~exist(fullfile(perframedir,[perframefns{j},'.mat']),'file');
  end
end

% experiments that need to compute all per-frame statistics (these are the
% ones that were archived or failed in sex classification)
nmissingfiles = cellfun(@(x)sum(double(x)),ismissing);
allmissing = nmissingfiles == numel(perframefns);

%% make list of experiments missing some per-frame features

tmp = find(idxcheck);
idxmissingpff = tmp(nmissingfiles > 0);
fid = fopen('expdirs_addperframefeatures_20120705.txt','w');
fprintf(fid,'%s\n',data0(idxmissingpff).file_system_path);
fclose(fid);

%% check that per-frame data is correct

% check some experiments that have everything
ncheck = 16;
idxcheck1 = find(nmissingfiles == 0);
[~,order] = sort({data(idxcheck1).exp_datetime});
idxcheck2 = idxcheck1(order(unique(round(linspace(1,numel(idxcheck1),ncheck)))));
fnscheck = {'a_mm','sex','dt','dnose2ell_angle_min20to20','dnose2ell','closestfly_angle_30tomin30','dnose2ell_angle_min30to30'};

expdirs = {data(idxcheck2).file_system_path};

isokperframedata = false(1,ncheck);
parfor ii = 1:ncheck,
  i = idxcheck2(ii);
  [~,experiment_name] = fileparts(expdirs{ii});
  outfilename = sprintf('perframefncheck_%02d_%s.txt',ii,experiment_name);
  
  isokperframedata(ii) = FlyBowlCheckPerFrameFeatures(expdirs{ii},'outfilename',outfilename,...
    'perframefns',fnscheck);
end

% check for trajectory files
fnscheck = {'a_mm','sex','dt'};
ncheck = 100;
idxcheck1 = setdiff(find([data.automated_pf]=='P' & ~allmissing),idxcheck2);
[~,order] = sort({data(idxcheck1).exp_datetime});
idxcheck2 = idxcheck1(order(unique(round(linspace(1,numel(idxcheck1),ncheck)))));

expdirs = {data(idxcheck2).file_system_path};

isokperframedata2 = false(1,ncheck);
parfor ii = 1:ncheck,
  i = idxcheck2(ii);
  [~,experiment_name] = fileparts(expdirs{ii});
  outfilename = sprintf('perframefncheck_%02d_%s.txt',ii,experiment_name);
  
  isokperframedata2(ii) = FlyBowlCheckPerFrameFeatures(expdirs{ii},'outfilename',outfilename,...
    'perframefns',fnscheck);
end

%% 

isfailure1 = strcmp({data0.automated_pf},'F') | strcmp({data0.manual_pf},'F');
idxcheck3 = find(isprimary & ~isfailure1);
[full_protocols,~,fullprotocolidx] = unique({data0(idxcheck3).full_analysis_protocol});

ncheckperprotocol = 5;
expdirs_check = {};
for i = 1:numel(full_protocols),
  idxcurr = idxcheck3(fullprotocolidx==i);
  if numel(idxcurr) <= ncheckperprotocol,
    expadd = {data0(idxcurr).file_system_path};
    expdirs_check = [expdirs_check,expadd];
  else
    [~,order] = sort({data0(idxcurr).exp_datetime});
    tmp = order(unique(round(linspace(1,numel(idxcurr),ncheckperprotocol))));
    expadd = {data0(idxcurr(tmp)).file_system_path};
    expdirs_check = [expdirs_check,expadd];
  end
  for j = 1:numel(expadd),
    [~,experiment_name] = fileparts(expadd{j});
    fprintf('%s %d: %s\n',full_protocols{i},j,experiment_name);
  end
end

fnscheck = {'a_mm','dt','anglesub','dnose2ell_angle_min20to20','dnose2ell','closestfly_angle_30tomin30','dnose2ell_angle_min30to30'};

isokperframedata3 = false(1,numel(expdirs_check));
parfor i = 1:numel(expdirs_check),
  [~,experiment_name] = fileparts(expdirs_check{i});
  outfilename = sprintf('perframefncheck_%02d_%s.txt',i,experiment_name);
  if exist(outfilename,'file'),
    continue;
  end
  try
    isokperframedata3(i) = FlyBowlCheckPerFrameFeatures(expdirs_check{i},'outfilename',outfilename,...
      'perframefns',fnscheck);
    fprintf('Check %d %s = %d\n',i,experiment_name,isokperframedata3(i));
  catch ME,
    warning('Error in check %d %s: %s\n',i,experiment_name,getReport(ME));
  end
end

% which ones passed?
[~,idxcheck4] = ismember(expdirs_check,{data0.file_system_path});
protocolcheck = {data0(idxcheck4).full_analysis_protocol};
for i = 1:numel(full_protocols),
  idxcurr = strcmp(protocolcheck,full_protocols{i});
  npass = nnz(isokperframedata3(idxcurr));
  nfail = nnz(~isokperframedata3(idxcurr));
  fprintf('*******\n%s: %d pass, %d fail\n',full_protocols{i},npass,nfail);
  idxfailure = find(idxcurr & ~isokperframedata3);
  for j = 1:numel(idxfailure),
    [~,experiment_name] = fileparts(expdirs_check{idxfailure(j)});
    outfilename = sprintf('perframefncheck_%02d_%s.txt',idxfailure(j),experiment_name);
    fprintf('\n%s:\n',outfilename);
    type(outfilename);
  end
  
  idxsuccess = find(idxcurr & isokperframedata3);
  for j = 1:numel(idxsuccess),
    [~,experiment_name] = fileparts(expdirs_check{idxsuccess(j)});
    fprintf('%s: SUCCESS\n',experiment_name);
  end

  
end

%%


isfailure1 = strcmp({data0.automated_pf},'F') | strcmp({data0.manual_pf},'F');
idxcheck3 = find(isprimary & ~isfailure1);
[full_protocols,~,fullprotocolidx] = unique({data0(idxcheck3).full_analysis_protocol});

ncheckperprotocol = 5;
expdirs_check = {};
for i = 1:numel(full_protocols),
  idxcurr = idxcheck3(fullprotocolidx==i);
  if numel(idxcurr) <= ncheckperprotocol,
    expadd = {data0(idxcurr).file_system_path};
    expdirs_check = [expdirs_check,expadd];
  else
    [~,order] = sort({data0(idxcurr).exp_datetime});
    tmp = order(unique(round(linspace(1,numel(idxcurr),ncheckperprotocol))));
    expadd = {data0(idxcurr(tmp)).file_system_path};
    expdirs_check = [expdirs_check,expadd];
  end
  for j = 1:numel(expadd),
    [~,experiment_name] = fileparts(expadd{j});
    fprintf('%s %d: %s\n',full_protocols{i},j,experiment_name);
  end
end

tmpoutputdir = '/groups/branson/bransonlab/projects/olympiad/HackHitData';
for i = 1:numel(expdirs_check),
  [~,experiment_name] = fileparts(expdirs_check{i});
  expdir = fullfile(tmpoutputdir,experiment_name);
  if ~exist(expdir,'dir'),
    fprintf('Copying %s...\n',experiment_name);
    SymbolicCopyExperimentDirectory(expdirs_check{i},tmpoutputdir);
  end
end  

isdone = false(1,numel(expdirs_check));
parfor i = 1:numel(expdirs_check),
  [~,experiment_name] = fileparts(expdirs_check{i});
  expdir = fullfile(tmpoutputdir,experiment_name);
  fprintf('Fixing %s...\n',expdir);
  try
    FlyBowlFixPerFrameFeatures(expdir,'analysis_protocol','20120706');
    isdone(i) = true;
  catch ME
    warning(getReport(ME));
  end
end

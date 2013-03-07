% ScriptCheckExperiments20120726

addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
addpath /groups/branson/bransonlab/projects/olympiad/SAGE/MATLABInterface/Trunk;
addpath /groups/branson/bransonlab/projects/olympiad/FlyBowlDataCapture;
addpath /groups/branson/home/bransonk/behavioranalysis/code/Jdetect/Jdetect/experiments/starvation;
addpath /groups/branson/home/bransonk/behavioranalysis/code/Jdetect/Jdetect/perframe;
addpath /groups/branson/home/bransonk/behavioranalysis/code/Jdetect/Jdetect/experiments/age;
settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
rootdatadir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlData';

%% parameters

analysis_protocol = 'current_non_olympiad_mushroombody';

%% get list of data

outexpdirfilename = '/groups/branson/home/bransonk/behavioranalysis/code/Jdetect/Jdetect/experiments/mbscreen/expdirs_mushroombody_20121008.txt';
if ~exist(outexpdirfilename,'file'),
  metadata = SAGEListBowlExperiments('screen_type','non_olympiad_mushroombody','checkflags',false,'rootdir',rootdatadir);
  fid = fopen(outexpdirfilename,'w');
  fprintf(fid,'%s\n',metadata.file_system_path);
  fclose(fid);
  expdirs = {metadata.file_system_path};
else
  expdirs = importdata(outexpdirfilename);
end

%% check that protocols are the same

for i = 1:numel(metadata),
  if mod(i,100) == 0,
    fprintf('%d / %d\n',i,numel(metadata));
  end
  filename = fullfile(metadata(i).file_system_path,'analysis_protocol.txt');
  didread = false;
  if exist(filename,'file'),
    try
      fid = fopen(filename,'r');
      metadata(i).full_analysis_protocol = fgetl(fid);
      didread = true;
      fclose(fid);
    catch ME
      error(getReport(ME));
    end
  end
  if ~didread,
    metadata(i).full_analysis_protocol = metadata(i).analysis_protocol;
  end
end

fns = {'ctrax_version','ctrax_settings','analysis_settings','analysis_version'};
for i = 1:numel(metadata),  
  s = regexp(metadata(i).full_analysis_protocol,',','split');
  ss = regexp(s,':','split');
  names = cellfun(@(x)x{1},ss,'UniformOutput',false);
  [ism,idx] = ismember(fns,names);
  for k = 1:numel(fns),
    if ism(k),
      metadata(i).(fns{k}) = ss{idx(k)}{2};
    else
      metadata(i).(fns{k}) = 'unknown';
    end
  end
end


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
isfailure = ( ([metadata.automated_pf] == 'F') & ...
  ~ismember({metadata.automated_pf_category},allowed_automated_pf_categories) ) | ...
  [metadata.manual_pf] == 'F';

%% find experiments that need to be retracked

retrack_min_exp_datetime = '20110601T000000';
retrack_max_ctrax_settings = '20111221';

retrack_max_ctrax_version = '0.2';
retrack_max_maxarea = 180;

retrack_idx = lexcmp({metadata.ctrax_settings},retrack_max_ctrax_settings) < 0 & ...
  lexcmp({metadata.exp_datetime},retrack_min_exp_datetime) > 0 & ...
  ~isfailure;

fid = fopen('expdirs_retrack_mushroombody_20121008.txt','w');
fprintf(fid,'%s\n',metadata(retrack_idx).file_system_path);
fclose(fid);

% spot check annfiles to make sure these are indeed old
metadata_check = metadata(retrack_idx);
ncheck = 50;
[~,order] = sort({metadata_check.exp_datetime});
idxcheck = order(unique(round(linspace(1,numel(metadata_check),ncheck))));

for ii = 1:ncheck,
  i = idxcheck(ii);
  annfile = fullfile(metadata_check(i).file_system_path,'movie.ufmf.ann');
  if ~exist(annfile,'file'),
    warning('Annotation file %s does not exist',annfile);
    continue;
  end
  [maxarea,version] = read_ann(annfile,'maxarea','version');
  if maxarea >= retrack_max_maxarea || ...
      lexcmp(version,retrack_max_ctrax_version) >= 0,
    warning('Area and version check failed for %s\n',metadata_check(i).experiment_name);
  end
end

%% find experiments that need to be fixed

idxfixpff = ~isfailure & ~retrack_idx;
fid = fopen('expdirs_fixperframefeatures_mushroombody_20121008.txt','w');
fprintf(fid,'%s\n',metadata(idxfixpff).file_system_path);
fclose(fid);

%% for now, also fix experiments that need to retracked

%% check that experiments are fixed

% check all experiments
mindatenum = 0;
% mindatenum = datenum('20120201','yyyymmdd');
for i = 1:numel(metadata),
  metadata(i).exp_datenum = datenum(metadata(i).exp_datetime,'yyyymmddTHHMMSS');
end
expidx = [metadata.exp_datenum] >= mindatenum;

CheckExperiments(expdirs(expidx),'checkprotocol',true,'checkfix',true,...
  'analysis_protocol',analysis_protocol);


%% check that scores are generated

classifierparamsfiles = {
  '/groups/branson/home/bransonk/behavioranalysis/code/Jdetect/Jdetect/experiments/JAABA_classifier_params1.txt'
  '/groups/branson/home/bransonk/behavioranalysis/code/Jdetect/Jdetect/experiments/JAABA_classifier_params2.txt'
  '/groups/branson/home/bransonk/behavioranalysis/code/Jdetect/Jdetect/experiments/JAABA_classifier_params3.txt'
  };

CheckExperiments(expdirs(expidx),'checkprotocol',false,'checkfix',false,...
  'classifierparamsfiles',classifierparamsfiles,...
  'analysis_protocol',analysis_protocol);

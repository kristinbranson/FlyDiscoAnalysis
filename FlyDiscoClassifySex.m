function [trx,summary_diagnostics,X] = FlyDiscoClassifySex(expdir,varargin)

version = '0.2';

%% parse parameters
[analysis_protocol,...
  settingsdir,...
  datalocparamsfilestr,...
  outdir,...
  dosave,...
  override_gender,...
  override_sex,...
  ~, ...
  ~, ...
  ~] = ...
    myparse(varargin,...
            'analysis_protocol','current_bubble',...
            'settingsdir',default_settings_folder_path(),...
            'datalocparamsfilestr','dataloc_params.txt',...
            'outdir',[],...
            'dosave',true,...
            'override_gender','', ...
            'override_sex','', ...
            'forcecompute', false, ...
            'debug', false, ...
            'do_run', []);

fprintf('Classifying sex for %s\n',expdir);

%% read in the data locations
datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
dataloc_params = ReadParams(datalocparamsfile);

% %%
% logger = PipelineLogger(expdir,mfilename(),...
%         dataloc_params,'classifysex_logfilestr',...
%         settingsdir,analysis_protocol,'versionstr',version);     

%% load the data
fprintf('Loading data...\n');

trxfile = fullfile(expdir,dataloc_params.trxfilestr);
load(trxfile,'trx');
% first and end frame
firstframe = min([trx.firstframe]);
endframe = max([trx.endframe]);

%% data locations
sexclassifierparamsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.sexclassifierparamsfilestr);
sexclassifierintxtfile = fullfile(settingsdir,analysis_protocol,dataloc_params.sexclassifiertxtfilestr);
sexclassifier_params = ReadParams(sexclassifierparamsfile);
sexclassifierin = ReadParams(sexclassifierintxtfile);

%% output files
sexclassifieroutmatfile = fullfile(expdir,outdir,...
  dataloc_params.sexclassifiermatfilestr);
sexclassifierdiagnosticsfile = fullfile(expdir,outdir,...
  dataloc_params.sexclassifierdiagnosticsfilestr);
trxfileout = fullfile(expdir,outdir,dataloc_params.trxfilestr);
tftrxappend = strcmp(trxfile,trxfileout);  % logical

%% compute areas
fprintf('Computing area...\n');

nflies = numel(trx); 

if nflies <= 0,
  error('No flies tracked in this video');
end

X = cell(1,nflies);
for fly = 1:nflies,
  % compute area
  area = (2*trx(fly).a_mm).*(2*trx(fly).b_mm)*pi;
  badidx = isinf(area) | isnan(area);
  if any(isnan(area)),
    warning('FlyBubbleClasisfySex:area',...
      'NaNs found in area for fly %d of experiment %s',fly,expdir);
    area(badidx) = inf;
  end
  areacurr = SmoothAreaOutliers(area,...
    sexclassifierin.areasmooth_filterorder,...
    sexclassifierin.areasmooth_maxfreq,...
    sexclassifierin.areasmooth_maxerrx);
  X{fly} = areacurr(:);
end

%% perform sex classification, if needed

% read metadata
metadatafile = fullfile(expdir,dataloc_params.metadatafilestr);
metadata = ReadMetadataFile(metadatafile);

% metadata.sex should be one of { 'b', 'm', 'f', 'x' }.  (Can also be
% uppercase.)  This field describes how the flies were prepared for the
% experiment.
%
% 'b' means there are both male and female flies, the user would like each 
%     fly's sex to be identified.  In this case, we use machine vision to
%     determine the sex of each fly.
% 'm' means all flies are known to be male.  In this case we simply label 
%     them all that way in registered_trx.mat, without using machine vision 
%     to determine the sex of each.
% 'f' means all flies are known to be female.  In this case we simply label 
%     them all that way in registered_trx.mat, without using machine vision 
%     to determine the sex of each.
% 'x' means flies were not presorted for sex, and the user does not want to 
%     perform sex classification.  In this case, we label each
%     fly's sex as 'x' in registered_trx.mat, meaning the sex is unknown.
%
% if the metadata lacks a 'sex' field, the deprecated 'gender' field is used as a
% fallback, for backward compatibility.  If neither field is present, a warning is
% raised, and 'b' is assumed,

% Sanitize sex field
sex = sanitize_metadata_sex(metadata, override_sex, override_gender) ;  % guaranteed to be one of {'b', 'm', 'f', 'x'}

% Run the machine vision algo to determine sex of each fly, or use the
% predetermined values if called for.
if strcmp(sex,'b') ,
  [trx, diagnostics, mu_area, var_area, state2sex, ll] = ...
    add_sex_field_to_trx_computed_using_machine_vision(trx, X, version, analysis_protocol, sexclassifier_params, sexclassifierin, ...
                                                       dosave, sexclassifieroutmatfile) ; 
else
  [trx, diagnostics, mu_area, var_area, state2sex, ll] = add_predetermined_sex_field_to_trx(trx, sex, X) ;
end

%% count number of flies, females, males
counts = struct;
counts.nfemales = zeros(1,endframe);
counts.nmales = zeros(1,endframe);
counts.nflies = zeros(1,endframe);
for fly = 1:numel(trx),
  isfemale = strcmp(trx(fly).sex,'F');
  ismale = strcmp(trx(fly).sex,'M');
  counts.nfemales(trx(fly).firstframe:trx(fly).endframe) = ...
    counts.nfemales(trx(fly).firstframe:trx(fly).endframe) + double(isfemale);
  counts.nmales(trx(fly).firstframe:trx(fly).endframe) = ...
    counts.nmales(trx(fly).firstframe:trx(fly).endframe) + double(ismale);
  counts.nflies(trx(fly).firstframe:trx(fly).endframe) = ...
    counts.nflies(trx(fly).firstframe:trx(fly).endframe) + 1;
end

% ignore part of video with no flies
fns = fieldnames(counts);
for i = 1:numel(fns),
  fn = fns{i};
  counts.(fn)(1:firstframe-1) = nan;
  counts.(fn)(endframe+1:end) = nan;
end

%% summary diagnostics

summary_diagnostics = struct;
ifemale = find(strcmp(state2sex,'F'));
imale = find(strcmp(state2sex,'M'));
summary_diagnostics.classifier_mu_area_female = mu_area(ifemale);
summary_diagnostics.classifier_mu_area_male = mu_area(imale);
summary_diagnostics.classifier_mu_area_female = mu_area(ifemale);
summary_diagnostics.classifier_mu_area_male = mu_area(imale);
summary_diagnostics.classifier_var_area_female = var_area(ifemale);
summary_diagnostics.classifier_var_area_male = var_area(imale);
if isempty(ll),
  summary_diagnostics.classifier_loglik = nan;
else
  summary_diagnostics.classifier_loglik = ll(end);
end
summary_diagnostics.classifier_niters = numel(ll);

fns = fieldnames(diagnostics);
for i = 1:numel(fns),
  fn = fns{i};
  summary_diagnostics.(sprintf('mean_%s',fn)) = nanmean([diagnostics.(fn)]); %#ok<NANMEAN> 
  summary_diagnostics.(sprintf('median_%s',fn)) = nanmedian([diagnostics.(fn)]); %#ok<NANMEDIAN> 
  summary_diagnostics.(sprintf('std_%s',fn)) = nanstd([diagnostics.(fn)],1); %#ok<NANSTD> 
  summary_diagnostics.(sprintf('min_%s',fn)) = min([diagnostics.(fn)]);
  summary_diagnostics.(sprintf('max_%s',fn)) = max([diagnostics.(fn)]);
end

fns = fieldnames(counts);
for i = 1:numel(fns),
  fn = fns{i};
  summary_diagnostics.(sprintf('mean_%s',fn)) = nanmean([counts.(fn)]); %#ok<NANMEAN> 
  summary_diagnostics.(sprintf('median_%s',fn)) = nanmedian([counts.(fn)]); %#ok<NANMEDIAN> 
  summary_diagnostics.(sprintf('std_%s',fn)) = nanstd([counts.(fn)],1); %#ok<NANSTD> 
  summary_diagnostics.(sprintf('min_%s',fn)) = min([counts.(fn)]);
  summary_diagnostics.(sprintf('max_%s',fn)) = max([counts.(fn)]);
end

fns1 = fieldnames(summary_diagnostics);
fprintf('Summary diagnostics:\n');
for i = 1:numel(fns1),
  fprintf('  %s: %f\n',fns1{i},summary_diagnostics.(fns1{i}));
end

%% write diagnostics

% sexclassifierinfo = logger.runInfo;
% sexclassifierinfo.version = version;

if dosave,
  if exist(sexclassifierdiagnosticsfile,'file'),
    try %#ok<TRYNC>
      delete(sexclassifierdiagnosticsfile);
    end
  end      
  fid = fopen(sexclassifierdiagnosticsfile,'w');
  if fid < 0,
    warning('FlyDiscoClassifySex:diags',...
      'Could not open file %s for writing, printing diagnostics to stdout instead',sexclassifierdiagnosticsfile);
    fid = 1;
  else
    cleaner = onCleanup(@()(fclose(fid))) ;
  end
else
  fid = 1;
end
fns = fieldnames(summary_diagnostics);
for i = 1:numel(fns),
  fprintf(fid,'%s,%f\n',fns{i},summary_diagnostics.(fns{i}));
end

%% resave

if dosave,
  try
    if tftrxappend
      save('-append',trxfileout,'trx');
    else
      tmp = load(trxfile);
      tmp.trx = trx;
      %tmp.sexclassifierinfo = sexclassifierinfo;
      save(trxfileout,'-struct','tmp');
    end
  catch
    try
      tmp = load(trxfile);
      trxfilebak = [trxfileout,'.bak'];
      unix(sprintf('mv %s %s',trxfileout,trxfilebak));
      tmp.trx = trx;
%       tmp.sexclassifierinfo = sexclassifierinfo;
      save(trxfileout,'-struct','tmp');
    catch ME,
      warning('FlyDiscoClassifySex:save',...
        'Could not save to file %s: %s',trxfileout,getReport(ME));
      fprintf('Could not save to file %s: %s',trxfileout,getReport(ME));
    end
  end
end

%logger.close();

end

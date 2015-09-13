function [trx,summary_diagnostics,X] = FlyBubbleClassifySex(expdir,varargin)

version = '0.2';

%% parse parameters
[analysis_protocol,...
  settingsdir,...
  datalocparamsfilestr,...
  outdir,...
  dosave] = myparse(varargin,...
  'analysis_protocol','current_bubble',...
  'settingsdir','/groups/branson/home/robiea/Code_versioned/FlyBubbleAnalysis/settings',...
  'datalocparamsfilestr','dataloc_params.txt',...
  'outdir',[],...
  'dosave',true);

fprintf('Classifying sex for %s\n',expdir);

%% read in the data locations
datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
dataloc_params = ReadParams(datalocparamsfile);

%%
logger = PipelineLogger(expdir,mfilename(),...
        dataloc_params,'classifysex_logfilestr',...
        settingsdir,analysis_protocol,'versionstr',version);     

%% load the data
logger.log('Loading data...\n');

trxfile = fullfile(expdir,dataloc_params.trxfilestr);
load(trxfile,'trx');
% first and end frame
firstframe = min([trx.firstframe]); %#ok<NODEF>
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
tftrxappend = strcmp(trxfile,trxfileout);

%% compute areas
logger.log('Computing area...\n');

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

%% fix sex if gender is not 'b'

% read gender
metadatafile = fullfile(expdir,dataloc_params.metadatafilestr);
metadata = ReadMetadataFile(metadatafile);

if isfield(metadata,'gender') && ~strcmpi(metadata.gender,'b'),
  logger.log('gender is not "b", not doing sex classification, just setting sex to %s for all flies',upper(metadata.gender));

  % set sex to metadata.gender for all flies
  % also set diagnostics that we can
  clear diagnostics
  mean_area_all = nanmean(cat(1,X{:}));
  var_area_all = nanstd(cat(1,X{:}),1);
  for fly = 1:numel(trx),
    trx(fly).sex = repmat({upper(metadata.gender)},[1,trx(fly).nframes]); %#ok<AGROW>
    diagnostics_curr = struct;
    diagnostics_curr.normhmmscore = nan;
    diagnostics_curr.nswaps = 0;
    diagnostics_curr.meanabsdev = nanmean(abs(X{fly}-mean_area_all));
    diagnostics(fly) = diagnostics_curr; %#ok<AGROW>
  end
  
  mu_area = nan(1,2);
  var_area = nan(1,2);
  if strcmpi(metadata.gender,'f'),
    mu_area(2) = mean_area_all;
    var_area(2) = var_area_all;
  elseif strcmpi(metadata.gender,'m'),
    mu_area(1) = mean_area_all;
    var_area(1) = var_area_all;
  end
  state2sex = {'M','F'};
  ll = [];
  
else
  
  %% learn a 2-state HMM for area in an unsupervised manner
  
  logger.log('gender = "b", learning 2-state HMM...\n');
  
  % initialize parameters
  nstates = 2;
  Psame = sexclassifierin.psame;
  ptrans = ones(nstates)-Psame;
  ptrans(eye(nstates)==1) = Psame;
  if isfield(sexclassifier_params,'frac_female'),
    prior = [1-sexclassifier_params.frac_female,sexclassifier_params.frac_female];
  else
    prior = ones(1,nstates)/nstates;
  end
  state2sex = cell(1,nstates);
  
  % em for hmm
  [mu_area,var_area,ll]=hmm_multiseq_1d(X,nstates,Psame,...
    sexclassifier_params.niters_em,sexclassifier_params.tol_em,...
    [],prior);
  if mu_area(1) > mu_area(2),
    mu_area = mu_area(end:-1:1);
    var_area = var_area(end:-1:1);
  end
  state2sex{argmax(mu_area)} = 'F';
  state2sex{argmin(mu_area)} = 'M';
  
  %% save classifier
  
  filterorder = sexclassifierin.areasmooth_filterorder; %#ok<NASGU>
  maxfreq = sexclassifierin.areasmooth_maxfreq; %#ok<NASGU>
  maxerrx = sexclassifierin.areasmooth_maxerrx; %#ok<NASGU>
  
  if dosave,
    try
      if exist(sexclassifieroutmatfile,'file'),
        delete(sexclassifieroutmatfile);
      end
    end
    try
      save(sexclassifieroutmatfile,'mu_area','var_area','ptrans','prior','ll',...
        'nstates','state2sex','maxerrx','maxfreq','filterorder','version','analysis_protocol');
    catch ME,
      warning('FlyBubbleClassifySex:save',...
        'Could not save to file %s: %s',sexclassifieroutmatfile,getReport(ME));
      logger.log('Could not save to file %s: %s\n',sexclassifieroutmatfile,getReport(ME));
    end
  end
  
  %% classify sex
  
  clear diagnostics;
  for fly = 1:numel(trx),
    
    % Viterbi to classify per-frame
    [trx(fly).sex,diagnostics(fly)] = ClassifySex(X{fly}',...
      mu_area,var_area,ptrans,state2sex); %#ok<AGROW>
    
  end
  
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
  summary_diagnostics.(sprintf('mean_%s',fn)) = nanmean([diagnostics.(fn)]);
  summary_diagnostics.(sprintf('median_%s',fn)) = nanmedian([diagnostics.(fn)]);
  summary_diagnostics.(sprintf('std_%s',fn)) = nanstd([diagnostics.(fn)],1);
  summary_diagnostics.(sprintf('min_%s',fn)) = min([diagnostics.(fn)]);
  summary_diagnostics.(sprintf('max_%s',fn)) = max([diagnostics.(fn)]);
end

fns = fieldnames(counts);
for i = 1:numel(fns),
  fn = fns{i};
  summary_diagnostics.(sprintf('mean_%s',fn)) = nanmean([counts.(fn)]);
  summary_diagnostics.(sprintf('median_%s',fn)) = nanmedian([counts.(fn)]);
  summary_diagnostics.(sprintf('std_%s',fn)) = nanstd([counts.(fn)],1);
  summary_diagnostics.(sprintf('min_%s',fn)) = min([counts.(fn)]);
  summary_diagnostics.(sprintf('max_%s',fn)) = max([counts.(fn)]);
end

fns1 = fieldnames(summary_diagnostics);
logger.log('Summary diagnostics:\n');
for i = 1:numel(fns1),
  logger.log('  %s: %f\n',fns1{i},summary_diagnostics.(fns1{i}));
end

%% write diagnostics

sexclassifierinfo = logger.runInfo;
sexclassifierinfo.version = version;

if dosave,
  if exist(sexclassifierdiagnosticsfile,'file'),
    try %#ok<TRYNC>
      delete(sexclassifierdiagnosticsfile);
    end
  end      
  fid = fopen(sexclassifierdiagnosticsfile,'w');
  if fid < 0,
    warning('FlyBubbleClassifySex:diags',...
      'Could not open file %s for writing, printing diagnostics to stdout instead',sexclassifierdiagnosticsfile);
    fid = 1;
  end
else
  fid = 1;
end
fns = fieldnames(summary_diagnostics);
for i = 1:numel(fns),
  fprintf(fid,'%s,%f\n',fns{i},summary_diagnostics.(fns{i}));
end
fns = fieldnames(sexclassifierinfo);
for i = 1:numel(fns),
  val = sexclassifierinfo.(fns{i});
  if ischar(val)
    fprintf(fid,'%s,%s\n',fns{i},val);
  end
end

if dosave && fid > 1,
  fclose(fid);
end

%% resave

if dosave,
  try
    if tftrxappend
      save('-append',trxfileout,'trx','sexclassifierinfo');
    else
      save(trxfileout,'trx','sexclassifierinfo');
    end
  catch %#ok<CTCH>
%     try
%       tmp = load(trxfile);
%       %delete(trxfile); % AL20150913 now have trxfileout
%       tmp.trx = trx;
%       tmp.sexclassifierinfo = sexclassifierinfo;
%       save(trxfileout,'-struct','tmp');
%     catch ME,
      warning('FlyBubbleClassifySex:save',...
        'Could not save to file %s: %s',trxfileout,getReport(ME));
      logger.log('Could not save to file %s: %s',trxfileout,getReport(ME));
%     end
  end
end

logger.close();

function [trx,summary_diagnostics,X] = FlyBowlClassifySex2(expdir,varargin)

%% parse parameters
[analysis_protocol,settingsdir,datalocparamsfilestr,dosave] = ...
  myparse(varargin,...
  'analysis_protocol','current',...
  'settingsdir','/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings',...
  'datalocparamsfilestr','dataloc_params.txt',...
  'dosave',true);

%% read in the data locations
datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
dataloc_params = ReadParams(datalocparamsfile);

%% load the data

trxfile = fullfile(expdir,dataloc_params.trxfilestr);
load(trxfile,'trx');
% first and end frame
firstframe = min([trx.firstframe]);
endframe = max([trx.endframe]);

%% data locations

sexclassifierparamsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.sexclassifierparamsfilestr);
sexclassifierintxtfile = fullfile(settingsdir,analysis_protocol,dataloc_params.sexclassifiertxtfilestr);
sexclassifieroutmatfile = fullfile(expdir,dataloc_params.sexclassifiermatfilestr);

sexclassifier_params = ReadParams(sexclassifierparamsfile);
sexclassifierin = ReadParams(sexclassifierintxtfile);

%% learn a 2-state HMM for area in an unsupervised manner

% initialize parameters
nstates = 2;
Psame = sexclassifierin.psame;
ptrans = ones(nstates)-Psame;
ptrans(eye(nstates)==1) = Psame;
prior = ones(1,nstates)/nstates; %#ok<NASGU>
state2sex = cell(1,nstates);

% break into smaller sequences where we're sure of area estimate
nflies = numel(trx); 
X = cell(1,nflies);
for fly = 1:nflies,
  % compute area
  area = (2*trx(fly).a_mm).*(2*trx(fly).b_mm)*pi;
  badidx = isinf(area) | isnan(area);
  if any(isnan(area)),
    warning('NaNs found in area for fly %d of experiment %s',fly,expdir);
    area(badidx) = inf;
  end
  areacurr = SmoothAreaOutliers(area,...
    sexclassifierin.areasmooth_filterorder,...
    sexclassifierin.areasmooth_maxfreq,...
    sexclassifierin.areasmooth_maxerrx);
  X{fly} = areacurr(:);
end

% em for hmm
[mu_area,var_area,ll]=hmm_multiseq_1d(X,nstates,Psame,...
  sexclassifier_params.niters_em,sexclassifier_params.tol_em);
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
  save(sexclassifieroutmatfile,'mu_area','var_area','ptrans','prior','ll',...
    'nstates','state2sex','maxerrx','maxfreq','filterorder');
end

%% classify sex

clear diagnostics;
for fly = 1:numel(trx), 
  
  % Viterbi to classify per-frame
  [trx(fly).sex,diagnostics(fly)] = ClassifySex(X{fly}',...
    mu_area,var_area,ptrans,state2sex); %#ok<AGROW>
  
end

%% count number of flies, females, males
counts = struct;
counts.nfemales = zeros(1,endframe);
counts.nmales = zeros(1,endframe);
for fly = 1:numel(trx),
  isfemale = strcmp(trx(fly).sex,'F');
  counts.nfemales(trx(fly).firstframe:trx(fly).endframe) = ...
    counts.nfemales(trx(fly).firstframe:trx(fly).endframe) + double(isfemale);
  counts.nmales(trx(fly).firstframe:trx(fly).endframe) = ...
    counts.nmales(trx(fly).firstframe:trx(fly).endframe) + double(~isfemale);
end
counts.nflies = counts.nfemales + counts.nmales;

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
summary_diagnostics.classifier_loglik = ll(end);
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

%% write diagnostics

sexclassifierdiagnosticsfile = fullfile(expdir,dataloc_params.sexclassifierdiagnosticsfilestr);
if dosave,
  fid = fopen(sexclassifierdiagnosticsfile,'w');
else
  fid = 1;
end
fns = fieldnames(summary_diagnostics);
for i = 1:numel(fns),
  fprintf(fid,'%s,%f\n',fns{i},summary_diagnostics.(fns{i}));
end

if dosave,
  fclose(fid);
end

%% resave

if dosave,
  save('-append',trxfile,'trx');
end
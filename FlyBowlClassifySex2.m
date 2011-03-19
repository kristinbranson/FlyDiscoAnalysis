function trx = FlyBowlClassifySex2(expdir,varargin)

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
nflies = numel(trx); %#ok<NODEF>
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
counts.nfemales = zeros(1,max([trx.endframe]));
counts.nmales = zeros(1,max([trx.endframe]));
for fly = 1:numel(trx),
  isfemale = strcmp(trx(fly).sex,'F');
  counts.nfemales(trx(fly).firstframe:trx(fly).endframe) = ...
    counts.nfemales(trx(fly).firstframe:trx(fly).endframe) + double(isfemale);
  counts.nmales(trx(fly).firstframe:trx(fly).endframe) = ...
    counts.nmales(trx(fly).firstframe:trx(fly).endframe) + double(~isfemale);
end
counts.nflies = counts.nfemales + counts.nmales;

%% write diagnostics
sexclassifierdiagnosticsfile = fullfile(expdir,dataloc_params.sexclassifierdiagnosticsfilestr);
if dosave,
  fid = fopen(sexclassifierdiagnosticsfile,'w');
else
  fid = 1;
end

ifemale = find(strcmp(state2sex,'F'));
imale = find(strcmp(state2sex,'M'));
fprintf(fid,'classifier_mu_area_female,%f\n',mu_area(ifemale));
fprintf(fid,'classifier_mu_area_male,%f\n',mu_area(imale));
fprintf(fid,'classifier_var_area_female,%f\n',var_area(ifemale));
fprintf(fid,'classifier_var_area_male,%f\n',var_area(imale));
fprintf(fid,'classifier_loglik,%f\n',ll(end));
fprintf(fid,'classifier_niters,%f\n',numel(ll));

fns = fieldnames(diagnostics);
for i = 1:numel(fns),
  fn = fns{i};
  fprintf(fid,'mean_%s,%f\n',fn,nanmean([diagnostics.(fn)]));
  fprintf(fid,'median_%s,%f\n',fn,nanmedian([diagnostics.(fn)]));
  fprintf(fid,'std_%s,%f\n',fn,nanstd([diagnostics.(fn)],1));
  fprintf(fid,'min_%s,%f\n',fn,min([diagnostics.(fn)]));
  fprintf(fid,'max_%s,%f\n',fn,max([diagnostics.(fn)]));
end

fns = fieldnames(counts);
for i = 1:numel(fns),
  fn = fns{i};
  fprintf(fid,'mean_%s,%f\n',fn,nanmean([counts.(fn)]));
  fprintf(fid,'median_%s,%f\n',fn,nanmedian([counts.(fn)]));
  fprintf(fid,'std_%s,%f\n',fn,nanstd([counts.(fn)],1));
  fprintf(fid,'min_%s,%f\n',fn,min([counts.(fn)]));
  fprintf(fid,'max_%s,%f\n',fn,max([counts.(fn)]));
end

if dosave,
  fclose(fid);
end

%% resave

if dosave,
  save('-append',trxfile,'trx');
end
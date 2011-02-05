function trx = FlyBowlClassifySex(expdir,varargin)

[analysis_protocol,settingsdir,datalocparamsfilestr] = ...
  myparse(varargin,...
  'analysis_protocol','current',...
  'settingsdir','/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings',...
  'datalocparamsfilestr','dataloc_params.txt');

%% read in the data locations
datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
dataloc_params = ReadParams(datalocparamsfile);
      
%% read in the sex classification parameters

sexclassifiermatfile = fullfile(settingsdir,analysis_protocol,dataloc_params.sexclassifiermatfilestr);
sexclassifier = load(sexclassifiermatfile);

%% load trx

trxfile = fullfile(expdir,dataloc_params.trxfilestr);
load(trxfile,'trx');

%% classify

clear diagnostics;
for fly = 1:numel(trx), %#ok<NODEF>

  % compute area
  area = (2*trx(fly).a_mm).*(2*trx(fly).b_mm)*pi;
  
  % smooth area
  areasmooth = SmoothAreaOutliers(area,...
    sexclassifier.filterorder,...
    sexclassifier.maxfreq,...
    sexclassifier.maxerrx);
  
  % Viterbi to classify per-frame
  [trx(fly).sex,diagnostics(fly)] = ClassifySex(areasmooth,...
    sexclassifier.mu_area,...
    sexclassifier.var_area,...
    sexclassifier.ptrans,...
    sexclassifier.state2sex); %#ok<AGROW>
  
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
fid = fopen(sexclassifierdiagnosticsfile,'w');
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


fclose(fid);

%% resave

save('-append',trxfile,'trx');
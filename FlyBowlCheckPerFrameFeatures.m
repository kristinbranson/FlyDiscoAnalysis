function success = FlyBowlCheckPerFrameFeatures(expdir,varargin)

success = true;

traj_fns = {'x','y','theta','a','b','timestamps',...
  'x_mm','y_mm','a_mm','b_mm','theta_mm','dt','sex'};
charfns = {'sex'};

perframefns = [traj_fns,{'dnose2ell_angle_min30to30','dnose2ell_angle_min20to20',...
  'dnose2ell_angle_30tomin30','dnose2ell','dell2nose','closestfly_nose2ell'}];

[analysis_protocol,settingsdir,datalocparamsfilestr,perframefns,perframefnsset,outfilename,fnskeep] = ...
  myparse(varargin,...
  'analysis_protocol','current',...
  'settingsdir','/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings',...
  'datalocparamsfilestr','dataloc_params.txt',...
  'perframefns', perframefns,...
  'perframenfnsset',[],...
  'outfilename','',...
  'fnskeep',{}...
  );

if ~isempty(perframefnsset),
  switch perframefnsset,
    case '20120810',
      perframefns = {'a_mm','dt','velmag_ctr','dnose2ell_angle_min20to20','magveldiff','dnose2ell','closestfly_angle_30tomin30','dnose2ell_angle_min30to30','du_tail'};
  end
end

if isempty(outfilename),
  fid = 1;
else
  fid = fopen(outfilename,'w');
end


%% load in trx


fprintf('Initializing trx...\n');

trx = Trx('analysis_protocol',analysis_protocol,'settingsdir',settingsdir,...
  'datalocparamsfilestr',datalocparamsfilestr);
rootdir = 'tmp_perframecheck';
if ~exist(rootdir,'dir'),
  mkdir(rootdir);
end
[~,experiment_name] = fileparts(expdir);
SymbolicCopyExperimentDirectory(expdir,rootdir);
tmpexpdir = fullfile(rootdir,experiment_name);
tmpperframedir = fullfile(tmpexpdir,trx.dataloc_params.perframedir);
tmp = dir(fullfile(tmpperframedir,'*.mat'));
for i = 1:numel(tmp),
  fn = regexprep(tmp(i).name,'\.mat$','');
  if ~ismember(fn,fnskeep),
    delete(fullfile(tmpperframedir,tmp(i).name));
  end
end
%delete(fullfile(tmpperframedir,'*.mat'));
trx.AddExpDir(tmpexpdir,'openmovie',false);
perframedir = fullfile(expdir,trx.dataloc_params.perframedir);

if ~iscell(perframefns) && ischar(perframefns) && strcmpi(perframefns,'all'),
  perframefns = dir(fullfile(perframedir,'*.mat'));
  perframefns = regexprep({perframefns.name},'\.mat$','');
  isscoreorlabel = ~cellfun(@isempty,regexp(perframefns,'([sS]core)|([lL]abel)','once'));
  perframefns(isscoreorlabel) = [];
end

for fni = 1:numel(perframefns),
  fn = perframefns{fni};
  filename = fullfile(perframedir,[fn,'.mat']);
  if ~exist(filename,'file'),
    fprintf(fid,'%s: file %s does not exist\n',fn,filename);
    continue;
  end
  tmp = load(filename,'data'); 
  dataold = tmp.data;
  fprintf('Computing %s...\n',fn);
  datanew = trx.(fn);
  if ndims(dataold) ~= ndims(datanew),
    success = false;
    fprintf(fid,'%s: cell ndims do not match.\n',fn);
    continue;
  end
  if ~all(size(dataold) == size(datanew)),
    success = false;
    fprintf(fid,'%s: cell size does not match.\n',fn);
    continue;
  end
  for i = 1:numel(dataold),
    if numel(dataold{i}) ~= numel(datanew{i}),
      success = false;
      fprintf(fid,'%s(%d) number of elements does not match.\n',fn,i);
      continue;
    end
    if ismember(fn,charfns),
      if ~all(strcmp(dataold{i},datanew{i})),
        success = false;
        fprintf(fid,'%s(%d): string mismatch\n',fn,i);
      end
    else
        
      maxerr = max(abs(dataold{i}(:)-datanew{i}(:)));
      maxv = max(max(dataold{i}(:)),max(datanew{i}(:)));
      minv = min(min(dataold{i}(:)),min(datanew{i}(:)));
      den = maxv - minv;
      if den <= 0,
        den = 1;
      end
      if maxerr/den > .000001,
        success = false;
        fprintf(fid,'%s(%d): maximum mismatch = %f\n',fn,i,maxerr);
      end
    end
  end  
end

if ~isempty(outfilename),
  fclose(fid);
end
  
rmdir(tmpexpdir,'s');

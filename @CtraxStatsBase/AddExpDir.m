% obj.AddExpDir(expdir)
% Adds the trajectories and all other data associated with experiment
% expdir to the data represented by obj. The movie associated with
% expdir is also opened for reading. If the trajectories for expdir
% have not yet been registered, then registration is performed and the
% results are saved to the trxfile. If expdir is already represented,
% it is removed then added again.  
function AddExpDir(obj,expdir)
if ismember(expdir,obj.expdir_bases),
  obj.RemoveExpDir(expdir);
end

obj.nexpdirs = obj.nexpdirs + 1;
n = obj.nexpdirs;

obj.expdir_bases{n} = expdir;
obj.expdirs{n} = fullfile(obj.rootdatadir,expdir);
if obj.NOWRITEACCESS,
  obj.write_expdirs{n} = fullfile(obj.resultsdir,expdir);
else
  obj.write_expdirs{n} = obj.expdirs{n};
end

% open movie
obj.moviefiles{n} = fullfile(obj.expdirs{n},obj.moviefilestr);
if ~exist(obj.moviefiles{n},'file'),
  error('Movie %s does not exist',obj.moviefiles{n});
end
[obj.readframes{n},obj.nframes(n),obj.moviefids(n),obj.headerinfos{n}] = ...
  get_readframe_fcn(obj.moviefiles{n});
im = obj.readframes{n}(1);
[obj.nrs(n),obj.ncs(n),obj.ncolors(n)] = size(im);

% read metadata
obj.metadatafiles{n} = fullfile(obj.expdirs{n},obj.metadatafilestr);
obj.metadata{n} = loadXMLDataTree(obj.metadatafiles{n});

% read temperature
obj.temperaturefiles{n} = fullfile(obj.expdirs{n},obj.temperaturefilestr);
if ~exist(obj.temperaturefiles{n},'file'),
  warning('Temperature file %s does not exist',obj.temperaturefiles{n});
  obj.temperaturefiles{n} = '';
  obj.temperaturestreams{n} = [];
else
  obj.temperaturestreams{n} = dlmread(obj.temperaturefiles{n});
end

% read annotation
obj.annfiles{n} = fullfile(obj.write_expdirs{n},obj.annfilestr);
if ~exist(obj.annfiles{n},'file');
  error('Annotation file %s does not exist',obj.annfiles{n});
end
obj.anns{n} = read_ann(obj.annfiles{n});
if isempty(fieldnames(obj.anns{n})),
  warning('Could not parse annotation file %s',obj.annfiles{n});
end

% resize images read from annotation
for i = 1:length(obj.annfile_images),
  fn = obj.annfile_images{i};
  if isfield(obj.anns{n},fn),
    obj.anns{n}.(fn) = reshape(obj.anns{n}.(fn),[obj.nrs(n),obj.ncs(n),numel(obj.anns{n}.(fn))/(obj.nrs(n)*obj.ncs(n))]);
  end
end

% read trajectories
obj.ctraxfiles{n} = fullfile(obj.write_expdirs{n},obj.ctraxfilestr);
obj.trxfiles{n} = fullfile(obj.write_expdirs{n},obj.trxfilestr);
obj.registrationfiles{n} = fullfile(obj.write_expdirs{n},obj.registrationfilestr);
if exist(obj.trxfiles{n},'file'),
  fprintf('Registration already done, loading from file %s\n',obj.registrationfiles{n});
  trx = load_tracks(obj.trxfiles{n});
else
  % register
  fprintf('Registering ...\n');
  
  % which bowl is this?
  if iscell(obj.bowl2MarkerAngle),
    bowl = obj.getMetaDataField('bowl','n',n);
    fprintf('Bowl = %s\n',bowl);
    i = find(strcmpi(bowl,obj.bowl2MarkerAngle(:,1)),1);
    if isempty(i),
      error('Bowl %s not in bowl2MarkerAngle',bowl);
    end
    bowlMarkerPairTheta_true = obj.bowl2MarkerAngle{i,2};
  end
  
  params = [fieldnames(obj.detectregistrationparams),struct2cell(obj.detectregistrationparams)]';
  params = params(:)';
  params(end+1:end+2) = {'bowlMarkerPairTheta_true',bowlMarkerPairTheta_true};
  if exist(obj.registrationfiles{n},'file'),
    inregistrationfile = obj.registrationfiles{n};
  else
    inregistrationfile = '';
  end
  trx = RegisterTrx(obj.ctraxfiles{n},...
    'annname',obj.annfiles{n},...
    'inregistrationfile',inregistrationfile,...
    'outregistrationfile',obj.registrationfiles{n},...
    'detectregistrationparams',params,...
    'outtrxname',obj.trxfiles{n}); 
  save(obj.trxfiles{n},'trx');
  fprintf('Registered trx stored to file %s\n',obj.registrationfiles{n});
end
registrationData = load(obj.registrationfiles{n});
obj.registrationData{n} = extractRegistrationFunctionAndTransform(registrationData);

% compute any necessary derived measurements
obj.landmarksfiles{n} = fullfile(obj.write_expdirs{n},obj.landmarksfilestr);
obj.closestflyfiles{n} = fullfile(obj.write_expdirs{n},obj.closestflyfilestr);
obj.speedfiles{n} = fullfile(obj.write_expdirs{n},obj.speedfilestr);

if n > 1 && obj.didComputeLandmarkMeasurements,
  trx = obj.ComputeLandmarkMeasurements(trx);
end
if n > 1 && obj.didComputeClosestFlyMeasurements,
  trx = obj.ComputeClosestFlyMeasurements(trx);
end
if n > 1 && obj.didComputeSpeedMeasurements,
  trx = obj.ComputeSpeedMeasurements(trx);
end

obj.trx = structappend(obj.trx,trx);
clear trx;

% number of flies
nfliesold = obj.nflies;
obj.nflies = length(obj.trx);
obj.nfliespermovie(n) = obj.nflies - nfliesold;

% indexing flies by movie
obj.movie2flies{n} = nfliesold+1:obj.nflies;
obj.fly2movie(nfliesold+1:obj.nflies) = n;

% movie size in mm
obj.width_mms(n) = obj.ncs(n) / obj.trx(obj.movie2flies{n}(1)).pxpermm;
obj.height_mms(n) = obj.nrs(n) / obj.trx(obj.movie2flies{n}(1)).pxpermm;
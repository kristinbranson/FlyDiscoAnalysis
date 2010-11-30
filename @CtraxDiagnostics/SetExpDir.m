function SetExpDir(obj,expdir)

if obj.isopen,
  obj = CloseExpDir(obj,expdir);
end

obj.expdir = expdir;

% open movie
obj.moviefile = fullfile(obj.expdir,obj.moviefilestr);
if ~exist(obj.moviefile,'file'),
  error('Movie %s does not exist',obj.moviefile);
end
[obj.readframe,obj.nframes,obj.moviefid,obj.headerinfo] = get_readframe_fcn(obj.moviefile);
im = obj.readframe(1);
obj.nr = size(im,1);
obj.nc = size(im,2);
obj.ncolors = size(im,3);

% read annotation
if obj.NOWRITEACCESS,
  tmp = regexp(obj.expdir,'/?([^/]+)/?$','tokens','once');
  obj.expdir_base = tmp{1};
  obj.annfile = fullfile(obj.resultsdir,obj.expdir_base,obj.annfilestr);
else
  obj.annfile = fullfile(obj.expdir,obj.annfilestr);
end
if ~exist(obj.annfile,'file');
  error('Annotation file %s does not exist',obj.annfile);
end
obj.ann = read_ann(obj.annfile);
if isempty(fieldnames(obj.ann)),
  warning('Could not parse annotation file %s',obj.annfile);
end

% resize images read from annotation
for i = 1:length(obj.ann_images),
  fn = obj.ann_images{i};
  if isfield(obj.ann,fn),
    obj.ann.(fn) = reshape(obj.ann.(fn),[obj.nr,obj.nc,numel(obj.ann.(fn))/(obj.nr*obj.nc)]);
  end
end

% read trajectories
obj.ctraxfile = fullfile(obj.expdir,obj.ctraxfilestr);
obj.trxfile = fullfile(obj.expdir,obj.trxfilestr);
if exist(trxfile,'file'),
  obj.trx = load_tracks(trxfile);
else
  obj.trx = load_tracks(ctraxfile);
  % TODO: implement calibration
  %obj.trx = calibrate_trx(obj.trx,'savefile',obj.trxfile);
end
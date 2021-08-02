% obj.AddExpDir(expdir)
% Adds data associated with experiment expdir to the data represented by obj.
function AddExpDir(obj,expdir,varargin)

[dooverwrite,openmovie,traj,tryloadwingtrx] = myparse(varargin,...
  'dooverwrite',true,...
  'openmovie',true,...
  'traj',[],...
  'tryloadwingtrx',true... % if true, load traj from dataloc_params.wingtrxfilestr
  );   

% remove trailing /
if expdir(end) == '/' || (ispc && expdir(end) == '\'),
  expdir = expdir(1:end-1);
end

if ismember(expdir,obj.expdirs),
  obj.RemoveExpDir(expdir);
end

obj.nexpdirs = obj.nexpdirs + 1;
n = obj.nexpdirs;

obj.expdirs{n} = expdir;

if openmovie && ~isempty(obj.dataloc_params.moviefilestr),
  moviename = fullfile(obj.expdirs{n},obj.dataloc_params.moviefilestr);
  [~,nframes,fid,vidinfo] = get_readframe_fcn(moviename);
  if ~isempty(fid) && ~isnan(fid) && fid > 1,
    fclose(fid);
  end

  % store video info
  obj.nrs(n) = vidinfo.nr;
  obj.ncs(n) = vidinfo.nc;
  obj.ncolors(n) = vidinfo.ncolors;
  obj.movie_nframes(n) = nframes;
end

% read trajectories
if tryloadwingtrx && isfield(obj.dataloc_params,'wingtrxfilestr') && ...
    exist(fullfile(obj.expdirs{n},obj.dataloc_params.wingtrxfilestr),'file'),
  trxfilestr = obj.dataloc_params.wingtrxfilestr;
else
  trxfilestr = obj.dataloc_params.trxfilestr;
end
obj.trxfiles{n} = fullfile(obj.expdirs{n},trxfilestr);
if isempty(traj),
  if ~exist(obj.trxfiles{n},'file'),
    error('Trajectory file %s does not exist',obj.trxfiles{n});
  end
  [traj,~,~,obj.movie_timestamps{n}] = load_tracks(obj.trxfiles{n});
end

% set movie properties when there is no movie
if ~openmovie || isempty(obj.dataloc_params.moviefilestr),
  obj.nrs(n) = max([traj.y]);
  obj.ncs(n) = max([traj.x]);
  obj.ncolors(n) = 0;
  obj.movie_nframes(n) = max([traj.endframe]);
end

% number of flies
obj.nfliespermovie(n) = length(traj);
nfliesold = obj.nflies;
obj.nflies = obj.nflies + obj.nfliespermovie(n);

% indexing flies by movie
obj.exp2flies{n} = nfliesold+1:obj.nflies;
obj.fly2exp(nfliesold+1:obj.nflies) = n;

% initialize data cache
obj.ndatacachedperexp(n) = 0;
obj.datacached{n} = {};
obj.fnscached{n} = {};
obj.nfnscached(n) = 0;

% store trajectories
obj.StoreTrajectories(n,traj,dooverwrite);

% movie size in mm
if isfield(obj,'pxpermm'),
  obj.width_mms(n) = obj.ncs(n) / obj.pxpermm(n);
  obj.height_mms(n) = obj.nrs(n) / obj.pxpermm(n);
else
  obj.width_mms(n) = nan;
  obj.height_mms(n) = nan;
end


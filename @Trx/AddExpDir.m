% obj.AddExpDir(expdir)
% Adds data associated with experiment expdir to the data represented by obj.
function AddExpDir(obj,expdir,varargin)

[dooverwrite] = myparse(varargin,'dooverwrite',true);

if ismember(expdir,obj.expdirs),
  obj.RemoveExpDir(expdir);
end

obj.nexpdirs = obj.nexpdirs + 1;
n = obj.nexpdirs;

obj.expdirs{n} = expdir;

moviename = fullfile(obj.expdirs{n},obj.dataloc_params.moviefilestr);
[~,nframes,fid,vidinfo] = get_readframe_fcn(moviename);
fclose(fid);

% store video info
obj.nrs(n) = vidinfo.nr;
obj.ncs(n) = vidinfo.nc;
obj.ncolors(n) = vidinfo.ncolors;
obj.movie_nframes(n) = nframes;

% read trajectories
obj.trxfiles{n} = fullfile(obj.expdirs{n},obj.dataloc_params.trxfilestr);
if ~exist(obj.trxfiles{n},'file'),
  error('Trajectory file %s does not exist',obj.trxfiles{n});
end
traj = load_tracks(obj.trxfiles{n});

% number of flies
obj.nfliespermovie(n) = length(traj);
nfliesold = obj.nflies;
obj.nflies = obj.nflies + obj.nfliespermovie(n);

% indexing flies by movie
obj.exp2flies{n} = nfliesold+1:obj.nflies;
obj.fly2exp(nfliesold+1:obj.nflies) = n;

% initialize data cache
obj.ndatacachedperexp(n) = 0;
obj.datacached{n} = struct;

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


% obj.AddExpDir(expdir)
% Adds data associated with experiment expdir to the data represented by obj.
function AddExpDir(obj,expdir,vidinfo)

if ismember(expdir,obj.expdir_bases),
  obj.RemoveExpDir(expdir);
end

obj.nexpdirs = obj.nexpdirs + 1;
n = obj.nexpdirs;

obj.expdir_bases{n} = expdir;
obj.expdir_reads{n} = fullfile(obj.rootreaddir,expdir);
obj.expdir_writes{n} = fullfile(obj.rootwritedir,expdir);

% store video info
obj.nrs(n) = vidinfo.nr;
obj.ncs(n) = vidinfo.nc;
obj.ncolors(n) = vidinfo.ncolors;
obj.movie_nframes(n) = vidinfo.nframes;

% read trajectories
obj.trxfiles{n} = fullfile(obj.expdir_writes{n},obj.trxfilestr);
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
obj.StoreTrajectories(n,traj);

% movie size in mm
obj.width_mms(n) = obj.ncs(n) / obj.pxpermm(n);
obj.height_mms(n) = obj.nrs(n) / obj.pxpermm(n);


% obj.RemoveExpDir(expdir)
% Removes the trajectories and all other data associated with experiment
% expdir to the data represented by obj. The movie file associated with
% expdir is also closed.

function RemoveExpDir(obj,expdir)

n = find(strcmp(expdir,obj.expdir_bases),1);
if isempty(n),
  error('Experiment directory %s is not currently loaded',expdir);
end

if numel(obj.moviefids) >= n,
  fclose(obj.moviefids(n));
end

if numel(obj.expdirs) >= n,
  obj.expdirs(n) = [];
end
if numel(obj.write_expdirs) >= n,
  obj.write_expdirs(n) = [];
end
if numel(obj.expdir_bases) >= n,
  obj.expdir_bases(n) = [];
end

if numel(obj.moviefiles) >= n,
  obj.moviefiles(n) = [];
end
if numel(obj.metadatafiles) >= n,
  obj.metadatafiles(n) = [];
end
if numel(obj.annfiles) >= n,
  obj.annfiles(n) = [];
end
if numel(obj.ctraxfiles) >= n,
  obj.ctraxfiles(n) = [];
end
if numel(obj.trxfiles) >= n,
  obj.trxfiles(n) = [];
end
if numel(obj.registrationfiles) >= n,
  obj.registrationfiles(n) = [];
end
if numel(obj.landmarksfiles) >= n,
  obj.landmarksfiles(n) = [];
end
if numel(obj.closestflyfiles) >= n,
  obj.closestflyfiles(n) = [];
end
if numel(obj.speedfiles) >= n,
  obj.speedfiles(n) = [];
end

% function handle for reading frames from current movie
if numel(obj.readframes) >= n,
  obj.readframes(n) = [];
end

% number of frames in current movie
if numel(obj.nframes) >= n,
  obj.nframes(n) = [];
end

% file pointer for current movie
if numel(obj.moviefids) >= n,
  obj.moviefids(n) = [];
end

% header info read from current movie
if numel(obj.headerinfos) >= n,
  obj.headerinfos(n) = [];
end

% metadata
if numel(obj.metadata) >= n,
  obj.metadata(n) = [];
end

% registration data
if numel(obj.registrationData) >= n,
  obj.registrationData(n) = [];
end

% movie frame size
if numel(obj.nrs) >= n,
  obj.nrs(n) = [];
end
if numel(obj.ncs) >= n,
  obj.ncs(n) = [];
end
if numel(obj.ncolors) >= n,
  obj.ncolors(n) = [];
end
if numel(obj.width_mms) >= n,
  obj.width_mms(n) = [];
end
if numel(obj.height_mms) >= n,
  obj.height_mms(n) = [];
end

% annotation info
if numel(obj.anns) >= n,
  obj.anns(n) = [];
end

if numel(obj.movie2flies) >= n,
  flies = obj.movie2flies{n};
  nfliescurr = length(flies);
else
  flies = [];
  nfliescurr = [];
end

% current trajectories
if numel(obj.movie2flies) >= n,
  obj.trx(flies) = [];

  % index maps
  obj.fly2movie(flies) = [];

  % indexes for everything after current movie need to be decremented
  for i = n+1:obj.nexpdirs,
    obj.movie2flies{i} = obj.movie2flies{i} - nfliescurr;
  end
  obj.movie2flies(n) = [];
end
if numel(obj.nfliespermovie) >= n,
  obj.nfliespermovie(n) = [];
end

obj.nexpdirs = obj.nexpdirs - 1;
obj.nflies = obj.nflies - nfliescurr;
function RemoveExpDir(obj,expdir)

% which experiment
n = find(strcmp(expdir,obj.expdir_bases),1);
if isempty(n),
  error('Experiment %s not loaded',expdir);
end

% clear file names
obj.expdir_bases(n) = [];
obj.expdir_reads(n) = [];
obj.expdir_writes(n) = [];

% clear video info
obj.nrs(n) = [];
obj.ncs(n) = [];
obj.ncolors(n) = [];
obj.nframes(n) = [];

% clear trajectory files
obj.trxfiles(n) = [];

% update number of flies
obj.nflies = obj.nflies - obj.nfliespermovie(n);
obj.nfliespermovie(n) = [];

% clear indexing info
obj.fly2exp(obj.exp2flies{n}) = [];
obj.exp2flies(n) = [];

% clear conversion from pixels to mm
obj.pxpermm(n) = [];

% clear movie size in mm
obj.width_mms(n) = [];
obj.height_mms(n) = [];

% clear cache
obj.datacached(n) = [];
obj.ndatacached = obj.ndatacached - obj.ndatacachedperexp(n);
obj.ndatacachedperexp(n) = [];

% update number of experiments
obj.nexpdirs = obj.nexpdirs - 1;

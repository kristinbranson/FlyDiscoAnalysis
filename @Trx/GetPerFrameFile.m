function filename = GetPerFrameFile(obj,fn,n)

dirname = fullfile(obj.expdir_writes{n},obj.dataloc_params.perframedir);
if ~exist(dirname,'file'),
  [success,msg,~] = mkdir(obj.expdir_writes{n},obj.dataloc_params.perframedir);
  if ~success,
    error('Error creating per-frame directory: %s',msg);
  end
end
filename = fullfile(dirname,[fn,'.mat']);

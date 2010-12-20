 % obj.DeleteRegistrationFiles([expdirs])
 % Deletes the files created during the registration step -- the
 % registerd trajectories (<expdir>/<trxfilestr>) and the registration
 % parameters (<expdir>/<trxfilestr>). If expdirs is input, then the
 % registration files for the specified expdirs will be deleted. If not
 % specified, then registration files for all loaded experiments will be
 % deleted.
 
function DeleteRegistrationFiles(obj,expdirs)

if ~exist('expdirs','var'),
  ns = 1:obj.nexpdirs;
else
  ns = obj.expdir2n(expdirs);
end

fprintf('Deleting the following registration files:\n');
for n = ns,
  if exist(obj.registrationfiles{n},'file'),
    fprintf('  %s\n',obj.registrationfiles{n});
    delete(obj.registrationfiles{n});
  end
  if exist(obj.trxfiles{n},'file') && ~strcmp(obj.trxfiles{n},obj.ctraxfiles{n}),
    fprintf('  %s\n',obj.trxfiles{n});
    delete(obj.trxfiles{n});    
  end
end
idx = find([data.automated_pf] == 'F');
[failed_screen_types,~,screenidx] = unique({data(idx).screen_type});


screeni = screeni+1;
fprintf('*** Screen type = %s ***\n',failed_screen_types{screeni});
tmp = data(idx(screenidx==screeni));

for i = 1:numel(tmp),
  [~,experiment_name] = fileparts(tmp(i).file_system_path);
  fprintf('%s\n',experiment_name);
  fprintf('automated_pf_category = %s\n',tmp(i).automated_pf_category);
  filename = fullfile(tmp(i).file_system_path,'automatic_checks_incoming_results.txt');
  if ~exist(filename,'file'),
    fprintf('automatic_checks_incoming_results.txt does not exist\n');
  else
    tmp2 = ReadParams(filename);
    if tmp2.automated_pf == 'U',
      fprintf('incoming checks passed\n');
    else
      fprintf(['incoming check notes = ',tmp2.notes_curation,'\n']);
    end
  end
  filename = fullfile(tmp(i).file_system_path,'automatic_checks_complete_results.txt');
  if ~exist(filename,'file'),
    fprintf('automatic_checks_complete_results.txt does not exist\n');
  else
    tmp2 = ReadParams(filename);
    if tmp2.automated_pf == 'U',
      fprintf('complete checks passed\n');
    else
      fprintf(['complete check notes = ',tmp2.notes_curation,'\n']);
    end
  end
  filename = fullfile(tmp(i).file_system_path,sprintf('ctrax_results_movie_%s.avi',experiment_name));
  fprintf('results movie exists = %d\n',exist(filename,'file'));
end
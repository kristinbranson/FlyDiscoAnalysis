function reset_test_experiment_folder(varargin)
  if length(varargin) >= 1 ,
    do_single_experiment = true ;
    experiment_name = varargin{1} ;
  else
    do_single_experiment = false ;
  end
  this_script_path = mfilename('fullpath') ;
  this_folder_path = fileparts(this_script_path) ;
  if do_single_experiment ,
    template_test_experiment_folder_path = fullfile(this_folder_path, 'analysis-test-template', experiment_name) ;
    test_experiment_folder_path = fullfile(this_folder_path, 'analysis-test-folder', experiment_name) ;
  else
    template_test_experiment_folder_path = fullfile(this_folder_path, 'analysis-test-template') ;
    test_experiment_folder_path = fullfile(this_folder_path, 'analysis-test-folder') ;
  end
  if exist(test_experiment_folder_path, 'file') ,
    try
      rmdir(test_experiment_folder_path, 's') ;  % often fails on issues with NFS lock files
    catch
    end
  end
  if exist(test_experiment_folder_path, 'file') ,
    %rmdir(test_experiment_folder_path, 's') ;  % issues with NFS lock files
    % Just delete all the non-hiddden files, this should be good enough
    command_line = ['find ' test_experiment_folder_path ' -not -path ''*/\.*'' -type f,l -exec rm ''{}'' \;'] ;
    system_with_error_handling(command_line) ;
  end
  if ~exist(test_experiment_folder_path, 'file') ,
    mkdir(test_experiment_folder_path) ;
  end
  copyfile(template_test_experiment_folder_path, test_experiment_folder_path) ;
end

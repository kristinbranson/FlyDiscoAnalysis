function result = did_acc_pass_for_experiment_folder(experiment_folder_name)

acc_output_file_name = fullfile(experiment_folder_name, 'automatic_checks_complete_info.mat') ;
if exist(acc_output_file_name, 'file') ,
  s = load('-mat', acc_output_file_name) ;
  if isfield(s, 'success') ,
    result = logical(s.success) ;
  else
    result = false ;
  end
else
  result = false ;
end

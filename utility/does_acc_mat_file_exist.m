function result = does_acc_mat_file_exist(experiment_folder_name)

acc_output_file_name = fullfile(experiment_folder_name, 'automatic_checks_complete_info.mat') ;
result = logical(exist(acc_output_file_name, 'file')) ;

function result = did_analysis_of_experiment_folder_succeed(experiment_folder_name)
% This function is deprecated.  Use did_acc_pass_for_experiment_folder(),
% which does the same thing but has a more-precise name.

warning('The function %s is deprecated.  Use did_acc_pass_for_experiment_folder.m, which does the same thing but has a more-precise name.', mfilename()) ;
result = did_acc_pass_for_experiment_folder(experiment_folder_name) ;

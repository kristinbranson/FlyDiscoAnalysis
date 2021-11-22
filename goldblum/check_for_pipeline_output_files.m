function all_tests_passed_from_experiment_index = ...
        check_for_pipeline_output_files(relative_path_to_folder_from_experiment_index, goldblum_destination_folder_path)
    % This function prints stuff to stdout when there are test failures
    successful_test_file_names = {'registered_trx.mat' 'wingtracking_results.mat' 'automatic_checks_complete_results.txt' } ;
    failed_test_file_names = { 'automatic_checks_complete_results.txt' } ;
    %excluded_test_file_names = { 'registered_trx.mat' 'wingtracking_results.mat' } ;    
    experiment_count = length(relative_path_to_folder_from_experiment_index) ;
    all_tests_passed_from_experiment_index = true(size(relative_path_to_folder_from_experiment_index)) ;
    for experiment_index = 1 : experiment_count ,
        relative_path_to_experiment_folder = relative_path_to_folder_from_experiment_index{experiment_index} ;
        experiment_folder_path = fullfile(goldblum_destination_folder_path, relative_path_to_experiment_folder) ;
        % First check for the file that indicates we wouldn't expect a experiment to
        % process fully
        is_experiment_expected_to_succeed = ~exist(fullfile(experiment_folder_path, 'flydisco-pipeline-should-fail-for-this-experiment'), 'file') ;
        if is_experiment_expected_to_succeed ,
            for i = 1 : length(successful_test_file_names) ,
                test_file_name = successful_test_file_names{i} ;
                test_file_path = ...
                    fullfile(experiment_folder_path, test_file_name) ;
                if ~exist(test_file_path, 'file') ,
                    fprintf('Test failure: No output file at %s\n', test_file_path) ;
                    all_tests_passed_from_experiment_index(experiment_index) = false ;
                end
            end
        else
            % Check for the files that *should* be there
            for i = 1 : length(failed_test_file_names) ,
                test_file_name = failed_test_file_names{i} ;
                test_file_path = ...
                    fullfile(experiment_folder_path, test_file_name) ;
                if ~exist(test_file_path, 'file') ,
                    fprintf('Test failure: No output file at %s\n', test_file_path) ;
                    all_tests_passed_from_experiment_index(experiment_index) = false ;
                end
            end
%             % Check for the files that should *not* be there
%             for i = 1 : length(excluded_test_file_names) ,
%                 test_file_name = excluded_test_file_names{i} ;
%                 test_file_path = ...
%                     fullfile(experiment_folder_path, test_file_name) ;
%                 if exist(test_file_path, 'file') ,
%                     fprintf('Test failure: Output file is at %s, but it shouldn''t be!\n', test_file_path) ;
%                     all_tests_passed_from_experiment_index(experiment_index) = false ;
%                 end
%             end            
        end
    end
end

function analyze_experiment_folders(folder_path_from_experiment_index, settings_folder_path, ...
                                    lab_head_last_name, ...
                                    do_use_bqueue, do_actually_submit_jobs, analysis_parameters)     

    % Figure out which experiments still need to be analyzed                                  
    experiment_count = length(folder_path_from_experiment_index) ;
    is_to_be_analyzed_from_experiment_index = true(experiment_count, 1) ;
    for i = 1 : length(folder_path_from_experiment_index) ,
        experiment_folder_path = folder_path_from_experiment_index{i} ;
        analysis_successful_file_path = fullfile(experiment_folder_path, 'ANALYSIS-COMPLETED-SUCCESSFULLY') ;
        analysis_failed_file_path = fullfile(experiment_folder_path, 'ANALYSIS-FAILED') ;        
        analysis_in_progress_file_path = fullfile(experiment_folder_path, 'ANALYSIS-IN-PROGRESS') ;
        aborted_file_path = fullfile(experiment_folder_path, 'ABORTED') ;
        is_to_be_analyzed_from_experiment_index(i) = ...
            ~logical(exist(analysis_successful_file_path, 'file')) && ~logical(exist(analysis_failed_file_path, 'file')) && ...
            ~logical(exist(analysis_in_progress_file_path, 'file')) && ~logical(exist(aborted_file_path, 'file')) ;
    end

    % Report how many experiments are left to be analyzed
    folder_path_from_to_be_analyzed_experiment_index = folder_path_from_experiment_index(is_to_be_analyzed_from_experiment_index) ;
    to_be_analyzed_experiment_count = length(folder_path_from_to_be_analyzed_experiment_index) ;
    fprintf('There are %d experiments that will be analyzed.\n', to_be_analyzed_experiment_count) ;
    if to_be_analyzed_experiment_count > 0 ,
        fprintf('Submitting these for analysis...\n') ;
    end
                                  
    if do_use_bqueue ,
        maxiumum_slot_count = 400 ;
        slots_per_job = 4 ;
        bqueue = bqueue_type(do_actually_submit_jobs, maxiumum_slot_count) ;

        % Queue the jobs
        for i = 1 : to_be_analyzed_experiment_count ,
            experiment_folder_path = folder_path_from_to_be_analyzed_experiment_index{i} ;
            [~, experiment_folder_name] = fileparts2(experiment_folder_path) ;
            stdouterr_file_path = find_free_log_file_name(experiment_folder_path) ;       
            bsub_options = sprintf('-P %s -J %s-flydisco-%s', lab_head_last_name, lab_head_last_name, experiment_folder_name) ;
            bqueue.enqueue(slots_per_job, ...
                           stdouterr_file_path, ...
                           bsub_options, ...
                           @fly_disco_analysis_pipeline_wrapper, ...
                               experiment_folder_path, ...
                               settings_folder_path, ...
                               analysis_parameters) ;
        end

        % Actually run the jobs
        maximum_wait_time = inf ;
        do_show_progress_bar = true ;
        tic_id = tic() ;
        job_statuses = bqueue.run(maximum_wait_time, do_show_progress_bar) ;
        toc(tic_id)
        
        % Report on any failed runs
        if all(job_statuses==1) ,
            % All is well
            fprintf('All jobs completed successfully.\n') ;
        else
            % Print the folders that had errors
            had_error = (job_statuses==-1) ;
            folder_path_from_errored_experiment_index = folder_path_from_to_be_analyzed_experiment_index(had_error) ;
            if ~isempty(folder_path_from_errored_experiment_index) ,                
                fprintf('The following experiment folders had errors:\n') ;
                for i = 1 : length(folder_path_from_errored_experiment_index) ,
                    experiment_folder_path = folder_path_from_errored_experiment_index{i} ;
                    fprintf('    %s\n', experiment_folder_path) ;
                end
                fprintf('\n') ;
            end
            
            % Print the folders that did not finish
            did_not_finish = (job_statuses==0) ;
            folder_path_from_unfinished_experiment_index = folder_path_from_to_be_analyzed_experiment_index(did_not_finish) ;
            if ~isempty(folder_path_from_unfinished_experiment_index) ,                
                fprintf('The following experiment folders did not finish processing in the alloted time:\n') ;
                for i = 1 : length(folder_path_from_unfinished_experiment_index) ,
                    experiment_folder_path = folder_path_from_unfinished_experiment_index{i} ;
                    fprintf('    %s\n', experiment_folder_path) ;
                end
                fprintf('\n') ;
            end
        end
    else
        for i = 1 : to_be_analyzed_experiment_count ,
            experiment_folder_path = folder_path_from_to_be_analyzed_experiment_index{i} ;
            fly_disco_analysis_pipeline_wrapper(experiment_folder_path, settings_folder_path, analysis_parameters) ;
        end
    end
end

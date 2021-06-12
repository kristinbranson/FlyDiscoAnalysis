function analyze_experiment_folders(folder_path_from_experiment_index, settings_folder_path, lab_head_last_name, ...
                                    do_force_analysis, do_use_bqueue, do_actually_submit_jobs, analysis_parameters)

    % Proces arguments                                
    if ~exist('do_force_analysis', 'var') || isempty(do_force_analysis) ,
        do_force_analysis = false ;
    end
    if ~exist('do_use_bqueue', 'var') || isempty(do_use_bqueue) ,
        do_use_bqueue = true ;
    end
    if ~exist('do_actually_submit_jobs', 'var') || isempty(do_actually_submit_jobs) ,
        do_actually_submit_jobs = true ;
    end
    if ~exist('analysis_parameters', 'var') || isempty(analysis_parameters) ,
        analysis_parameters = cell(1,0) ;
    end

    % If do_force_analysis is true, clear any files indicating ongoing run
    experiment_count = length(folder_path_from_experiment_index) ;
    if do_force_analysis ,
        for i = 1 : experiment_count ,
            experiment_folder_path = folder_path_from_experiment_index{i} ;
            analysis_in_progress_file_path = fullfile(experiment_folder_path, 'ANALYSIS-IN-PROGRESS') ;
            if exist(analysis_in_progress_file_path, 'file') ,
                delete(analysis_in_progress_file_path) ;
            end
        end
    end
    
    % We don't analyze experiments that are already being analyzed, ones where
    % the experiment was aborted during data-taking
    is_to_be_analyzed_from_experiment_index = true(experiment_count, 1) ;
    for i = 1 : experiment_count ,
        experiment_folder_path = folder_path_from_experiment_index{i} ;
        analysis_in_progress_file_path = fullfile(experiment_folder_path, 'ANALYSIS-IN-PROGRESS') ;
        aborted_file_path = fullfile(experiment_folder_path, 'ABORTED') ;
        is_to_be_skipped = ...
          exist(analysis_in_progress_file_path, 'file') || exist(aborted_file_path, 'file') ;
        is_to_be_analyzed_from_experiment_index(i) = ~is_to_be_skipped ;
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
        slots_per_job = 2 ;
        bqueue = bqueue_type(do_actually_submit_jobs, maxiumum_slot_count) ;

        % Queue the jobs
        for i = 1 : to_be_analyzed_experiment_count ,
            experiment_folder_path = folder_path_from_to_be_analyzed_experiment_index{i} ;
            [~, experiment_folder_name] = fileparts2(experiment_folder_path) ;
            % We use the options list to pass the stdout/stderr file, b/c the
            % usual mechanism doesn't support appending.
            stdouterr_file_path = fullfile(experiment_folder_path, 'flydisco-analysis-log.txt') ;
            bsub_options = sprintf('-P %s -J %s-flydisco-%s -o %s -e %s', ...
                                   lab_head_last_name, ...
                                   lab_head_last_name, ...
                                   experiment_folder_name, ...
                                   stdouterr_file_path, ...
                                   stdouterr_file_path) ;
            bqueue.enqueue(slots_per_job, ...
                           [], ...
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

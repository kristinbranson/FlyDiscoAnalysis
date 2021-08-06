function goldblum_analyze_experiment_folders(folder_path_from_experiment_index, settings_folder_path, lab_head_last_name, ...
                                             do_force_analysis, do_use_bqueue, do_actually_submit_jobs, analysis_parameters)

    % Process arguments                                
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

    % Specify bsub parameters
    maxiumum_slot_count = 400 ;
    slots_per_job = 4 ;
    
%     % If do_force_analysis is true, clear any files indicating ongoing run
%     experiment_count = length(folder_path_from_experiment_index) ;
%     if do_force_analysis ,
%         for i = 1 : experiment_count ,
%             experiment_folder_path = folder_path_from_experiment_index{i} ;
%             analysis_in_progress_file_path = fullfile(experiment_folder_path, 'PIPELINE-IN-PROGRESS') ;
%             try
%                 ensure_file_does_not_exist(analysis_in_progress_file_path) ;
%             catch me
%                 fprintf('Tried to delete the file %s (if it exists), but something went wrong.  Proceeding nevertheless.\n', analysis_in_progress_file_path) ;
%                 fprintf('Here''s some information about what went wrong:\n') ;
%                 fprintf('%s\n', me.getReport()) ;                
%             end
%         end
%     end
    
%     % We don't analyze experiments that are already being analyzed
%     is_to_be_analyzed_from_experiment_index = true(experiment_count, 1) ;
%     for i = 1 : experiment_count ,
%         experiment_folder_path = folder_path_from_experiment_index{i} ;
%         analysis_in_progress_file_path = fullfile(experiment_folder_path, 'PIPELINE-IN-PROGRESS') ;
%         is_to_be_skipped = ...
%           exist(analysis_in_progress_file_path, 'file') ;
%         is_to_be_analyzed_from_experiment_index(i) = ~is_to_be_skipped ;
%     end

    % Report how many experiments are to be analyzed
    experiment_count = length(folder_path_from_experiment_index) ;
    fprintf('There are %d experiments that will be analyzed.\n', experiment_count) ;
    if experiment_count > 0 ,
        fprintf('Submitting these for analysis...\n') ;
    end

    % Run goldblum_FlyDiscoPipeline_wrapper() on all experiments
    if do_use_bqueue ,
        bqueue = bqueue_type(do_actually_submit_jobs, maxiumum_slot_count) ;

        % Queue the jobs
        for i = 1 : experiment_count ,
            experiment_folder_path = folder_path_from_experiment_index{i} ;
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
                           @goldblum_FlyDiscoPipeline_wrapper, ...
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
        successful_job_count = sum(job_statuses==1) ;
        errored_job_count = sum(job_statuses==-1) ;
        did_not_finish_job_count = sum(job_statuses==0) ;
        if experiment_count == successful_job_count ,
            % All is well
            fprintf('All %d jobs completed successfully.\n', successful_job_count) ;
        else
            % Print the folders that completed successfully
            did_complete_successfully = (job_statuses==+1) ;
            folder_path_from_successful_experiment_index = folder_path_from_experiment_index(did_complete_successfully) ;
            if ~isempty(folder_path_from_successful_experiment_index) ,                
                fprintf('These %d jobs completed successfully:\n', successful_job_count) ;
                for i = 1 : length(folder_path_from_successful_experiment_index) ,
                    experiment_folder_path = folder_path_from_successful_experiment_index{i} ;
                    fprintf('    %s\n', experiment_folder_path) ;
                end
                fprintf('\n') ;
            end            
            
            % Print the folders that had errors
            had_error = (job_statuses==-1) ;
            folder_path_from_errored_experiment_index = folder_path_from_experiment_index(had_error) ;
            if ~isempty(folder_path_from_errored_experiment_index) ,                
                fprintf('These %d experiment folders had errors:\n', errored_job_count) ;
                for i = 1 : length(folder_path_from_errored_experiment_index) ,
                    experiment_folder_path = folder_path_from_errored_experiment_index{i} ;
                    fprintf('    %s\n', experiment_folder_path) ;
                end
                fprintf('\n') ;
            end
            
            % Print the folders that did not finish
            did_not_finish = (job_statuses==0) ;
            folder_path_from_unfinished_experiment_index = folder_path_from_experiment_index(did_not_finish) ;
            if ~isempty(folder_path_from_unfinished_experiment_index) ,                
                fprintf('These %d experiment folders did not finish processing in the alloted time:\n', did_not_finish_job_count) ;
                for i = 1 : length(folder_path_from_unfinished_experiment_index) ,
                    experiment_folder_path = folder_path_from_unfinished_experiment_index{i} ;
                    fprintf('    %s\n', experiment_folder_path) ;
                end
                fprintf('\n') ;
            end
        end
    else
        % If not using bqueue, just run them normally (usually just for debugging)
        job_statuses = nan(1, experiment_count) ;
        for i = 1 : experiment_count ,
            experiment_folder_path = folder_path_from_experiment_index{i} ;
            goldblum_FlyDiscoPipeline_wrapper(experiment_folder_path, settings_folder_path, analysis_parameters) ;
            job_statuses(i) = +1 ;  % Indicates completed sucessfully
        end
    end

    
    
    %
    % "Caboose" phase
    %
        
    % Run the caboose jobs
    if do_use_bqueue ,
        bqueue = bqueue_type(do_actually_submit_jobs, maxiumum_slot_count) ;

        % Queue the jobs
        for i = 1 : experiment_count ,
            experiment_folder_path = folder_path_from_experiment_index{i} ;
            [~, experiment_folder_name] = fileparts2(experiment_folder_path) ;
            % We use the options list to pass the stdout/stderr file, b/c the
            % usual mechanism doesn't support appending.
            stdouterr_file_path = fullfile(experiment_folder_path, 'flydisco-analysis-log.txt') ;
            bsub_options = sprintf('-P %s -J %s-flydisco-caboose-%s -o %s -e %s', ...
                                   lab_head_last_name, ...
                                   lab_head_last_name, ...
                                   experiment_folder_name, ...
                                   stdouterr_file_path, ...
                                   stdouterr_file_path) ;
            bqueue.enqueue(slots_per_job, ...
                           [], ...
                           bsub_options, ...
                           @goldblum_FlyDiscoCaboose_wrapper, ...
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
        successful_job_count = sum(job_statuses==1) ;
        errored_job_count = sum(job_statuses==-1) ;
        did_not_finish_job_count = sum(job_statuses==0) ;
        if experiment_count == successful_job_count ,
            % All is well
            fprintf('All %d caboose jobs completed successfully.\n', successful_job_count) ;
        else
            % Print the folders that completed successfully
            did_complete_successfully = (job_statuses==+1) ;
            folder_path_from_successful_experiment_index = folder_path_from_experiment_index(did_complete_successfully) ;
            if ~isempty(folder_path_from_successful_experiment_index) ,                
                fprintf('These %d caboose jobs completed successfully:\n', successful_job_count) ;
                for i = 1 : length(folder_path_from_successful_experiment_index) ,
                    experiment_folder_path = folder_path_from_successful_experiment_index{i} ;
                    fprintf('    %s\n', experiment_folder_path) ;
                end
                fprintf('\n') ;
            end            
            
            % Print the folders that had errors
            had_error = (job_statuses==-1) ;
            folder_path_from_errored_experiment_index = folder_path_from_experiment_index(had_error) ;
            if ~isempty(folder_path_from_errored_experiment_index) ,                
                fprintf('These %d experiment folders had errors during the caboose phase:\n', errored_job_count) ;
                for i = 1 : length(folder_path_from_errored_experiment_index) ,
                    experiment_folder_path = folder_path_from_errored_experiment_index{i} ;
                    fprintf('    %s\n', experiment_folder_path) ;
                end
                fprintf('\n') ;
            end
            
            % Print the folders that did not finish
            did_not_finish = (job_statuses==0) ;
            folder_path_from_unfinished_experiment_index = folder_path_from_experiment_index(did_not_finish) ;
            if ~isempty(folder_path_from_unfinished_experiment_index) ,                
                fprintf('These %d experiment folders did not finish processing in the alloted time during the caboose phase:\n', did_not_finish_job_count) ;
                for i = 1 : length(folder_path_from_unfinished_experiment_index) ,
                    experiment_folder_path = folder_path_from_unfinished_experiment_index{i} ;
                    fprintf('    %s\n', experiment_folder_path) ;
                end
                fprintf('\n') ;
            end
        end
    else
        % If not using bqueue, just run them normally (usually just for debugging)
        job_statuses = nan(1, experiment_count) ;
        for i = 1 : experiment_count ,
            experiment_folder_path = folder_path_from_experiment_index{i} ;
            goldblum_FlyDiscoCaboose_wrapper(experiment_folder_path, settings_folder_path, analysis_parameters) ;
            job_statuses(i) = +1 ;  % Indicates completed sucessfully
        end
    end    
end

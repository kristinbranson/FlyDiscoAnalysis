function test_transfero_FlyDiscoPipeline_wrapper_on_experiment_list(folder_path_from_experiment_index, ...
                                                                    cluster_billing_account_name, ...
                                                                    user_name_for_configuration_purposes, ...
                                                                    do_use_bqueue, ...
                                                                    do_actually_submit_jobs, ...
                                                                    ssh_host_name, ...
                                                                    varargin)
         
    %goldblum_analyze_experiment_folders  Runs the FlyDisco pipeline on a set of
    %                                     experiment folders.
    %
    %   goldblum_analyze_experiment_folders(experiment_folder_list) runs the
    %   FlyDisco pipeline on the experiment folders given in the cell array of
    %   strings experiment_folder_list.  This is the function used by Goldblum to
    %   run experiments through the FlyDisco pipeline.  Note in particular that this
    %   function first runs the pipeline on all the experiments with the
    %   automaticcheckscomplete stage turned off, then runs them all with *only* the
    %   automaticcheckscomplete stage turned on.  By default, each run of the
    %   pipeline is submitted as a single LSF cluster job.  The function
    %   FlyDiscoPipeline() is used to run the pipeline.  See the documentation for
    %   that function for more details.
    %
    %   goldblum_analyze_experiment_folders(experiment_folder_list, settings_folder_path)
    %   uses analysis-protocol folders drawn from settings_folder_path instead of
    %   the default setting folder path.  Again, see the documentation of
    %   FlyDiscoPipeline() for more details.
    %
    %   goldblum_analyze_experiment_folders(..., cluster_billing_account_name)
    %   bills any jobs submitted to the LSF cluster to the account specified by 
    %   the string cluster_billing_account_name.  Examples might be 'branson',
    %   'rubin', and 'scicompsoft'.  
    %   
    %   goldblum_analyze_experiment_folders(..., do_use_bqueue), if do_use_bqueue is
    %   true, uses the bqueue_type() framework for submitting jobs to the LSF
    %   cluster.  If do_use_bqueue is false, jobs are run locally, and without first
    %   submitting them to a bqueue.  This option is mostly useful for debugging.
    %   If missing or empty, do_use_bqueue defaults to true.
    %
    %   goldblum_analyze_experiment_folders(..., do_actually_submit_jobs), if 
    %   do_use_bqueue is true, determines whether jobs are actually submitted to the
    %   LSF queue or are simply run locally.  Again, this is mostly useful for
    %   debugging purposes.  If do_use_bqueue is false, this argument is ignored.  
    %   If missing or empty, do_actually_submit_jobs defaults to true.
    %
    %   goldblum_analyze_experiment_folders(..., analysis_parameters_as_name_value_list)
    %   allows the caller to specify additional (key, value) pairs that are passed to
    %   FlyDiscoPipeline().  For instance, a caller might set this argument to 
    %   {'docomputeperframestats', 'off'} to disable the computeperframestats stage
    %   of FlyDiscoPipeline().  See the documentation of FlyDiscoPipeline() for more
    %   details, including a complete list of the supported keys.  If missing or
    %   empty, analysis_parameters_as_name_value_list defaults to cell(1,0).

    % Process arguments                                
    if ~exist('do_use_bqueue', 'var') || isempty(do_use_bqueue) ,
        do_use_bqueue = true ;
    end
    if ~exist('do_actually_submit_jobs', 'var') || isempty(do_actually_submit_jobs) ,
        do_actually_submit_jobs = true ;
    end
    if ~exist('ssh_host_name', 'var') || isempty(ssh_host_name) ,
        ssh_host_name = '' ;
    end
    
    optional_arguments_as_name_value_list = varargin ;

    % Specify bsub parameters
    maxiumum_slot_count = 400 ;
    slots_per_job = 4 ;
    do_use_xvfb = true ;  % Matlab on linux leaks memory when you call getframe() without an X11 server   

    % Report how many experiments are to be analyzed
    experiment_count = length(folder_path_from_experiment_index) ;
    fprintf('There are %d experiments that will be analyzed.\n', experiment_count) ;
    if experiment_count > 0 ,
        fprintf('Submitting these for analysis...\n') ;
    end

    % Run transfero_FlyDiscoPipeline_wrapper() on all experiments
    if do_use_bqueue ,
        bqueue = bqueue_type(do_actually_submit_jobs, maxiumum_slot_count, do_use_xvfb, ssh_host_name) ;

        % Queue the jobs
        for i = 1 : experiment_count ,
            experiment_folder_path = folder_path_from_experiment_index{i} ;
            [~, experiment_folder_name] = fileparts2(experiment_folder_path) ;
            % We use the options list to pass the stdout/stderr file, b/c the
            % usual mechanism doesn't support appending.
            stdouterr_file_path = fullfile(experiment_folder_path, 'flydisco-analysis-log.txt') ;
            bsub_options = sprintf('-P %s -J %s-flydisco-%s -o %s -e %s', ...
                                   cluster_billing_account_name, ...
                                   cluster_billing_account_name, ...
                                   experiment_folder_name, ...
                                   stdouterr_file_path, ...
                                   stdouterr_file_path) ;
            bqueue.enqueue(slots_per_job, ...
                           [], ...
                           bsub_options, ...
                           @transfero_FlyDiscoPipeline_wrapper, ...
                               experiment_folder_path, ...
                               user_name_for_configuration_purposes, ...
                               optional_arguments_as_name_value_list{:}) ;
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
            transfero_FlyDiscoPipeline_wrapper(experiment_folder_path, user_name_for_configuration_purposes, optional_arguments_as_name_value_list{:}) ;
            job_statuses(i) = +1 ;  % Indicates completed sucessfully
        end
    end
end

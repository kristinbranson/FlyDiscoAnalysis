function run_transfero_FlyDiscoPipeline_wrapper_on_experiment_list(folder_path_from_experiment_index, ...
                                                                   cluster_billing_account_name, ...
                                                                   user_name_for_configuration_purposes, ...
                                                                   do_use_bqueue, ...
                                                                   do_actually_submit_jobs, ...
                                                                   do_try, ...
                                                                   ssh_host_name, ...
                                                                   varargin)
         
    %run_transfero_FlyDiscoPipeline_wrapper_on_experiment_list
    %   Runs the FlyDisco pipeline on a set of experiment folders, optionally on
    %   an LSF compute cluster.
    %
    %   run_transfero_FlyDiscoPipeline_wrapper_on_experiment_list(experiment_folder_list)
    %   runs the FlyDisco pipeline on the experiment folders given in the cell array
    %   of strings experiment_folder_list.  This function uses the function
    %   transfero_FlyDiscoPipeline_wrapper(), which is typically used by Transfero
    %   to run experiments through the FlyDisco pipeline.  Note in particular that
    %   this function first runs the pipeline on all the experiments with the
    %   automaticcheckscomplete stage turned off, then runs them all with *only* the
    %   automaticcheckscomplete stage turned on.  By default, each run of the
    %   pipeline is submitted as a single LSF cluster job.  The function
    %   FlyDiscoPipeline() is used to run the pipeline.  See the documentation for
    %   that function for more details.
    %
    %   run_transfero_FlyDiscoPipeline_wrapper_on_experiment_list(...,
    %   cluster_billing_account_name) bills any jobs submitted to the LSF cluster to
    %   the account specified by the string cluster_billing_account_name.  Examples
    %   might be 'branson', 'rubin', and 'scicompsoft'.  If empty, the cluster
    %   billing account name will be looked for in the configuration file
    %   specified by user_name_for_configuration_purposes.
    %   
    %   run_transfero_FlyDiscoPipeline_wrapper_on_experiment_list(...,
    %   user_name_for_configuration_purposes) specifies the user name used for
    %   configuration purposes.  This is normally the name of '*lab' or '*robot'
    %   account used for automated runs.  It is used to look up the cluster
    %   billing account name and the settings folder.  If
    %   cluster_billing_account_name is non-empty, it will override the setting
    %   specified in the lab configuration file.
    %   
    %   run_transfero_FlyDiscoPipeline_wrapper_on_experiment_list(...,
    %   do_use_bqueue), if do_use_bqueue is true, uses the bqueue_type() framework
    %   for submitting jobs to the LSF cluster.  If do_use_bqueue is false, jobs are
    %   run locally, and without first submitting them to a bqueue.  This option is
    %   mostly useful for debugging. If missing or empty, do_use_bqueue defaults to
    %   true.
    %
    %   run_transfero_FlyDiscoPipeline_wrapper_on_experiment_list(...,
    %   do_actually_submit_jobs), if do_use_bqueue is true, determines whether jobs
    %   are actually submitted to the LSF queue or are simply run locally.  Again,
    %   this is mostly useful for debugging purposes.  If do_use_bqueue is false,
    %   this argument is ignored. If missing or empty, do_actually_submit_jobs
    %   defaults to true.
    %
    %   run_transfero_FlyDiscoPipeline_wrapper_on_experiment_list(..., 
    %   do_try), if do_use_bqueue is false, determines whether jobs are run inside
    %   a try-catch block or not.  Setting do_try==false is mostly useful for
    %   debugging purposes.  If do_use_bqueue is true, this argument is ignored.
    %   If missing or empty, do_try defaults to true.
    %
    %   run_transfero_FlyDiscoPipeline_wrapper_on_experiment_list(...,
    %   ssh_host_name), if ssh_host_name is nonempty, issues all LSF commands (bsub,
    %   bjobs, etc) through the named ssh_host_name.  This makes it possible to
    %   submit and monitor cluster jobs even if the local machine is not an LSF
    %   submit host.
    %
    %   run_transfero_FlyDiscoPipeline_wrapper_on_experiment_list(..., 'param1',
    %   val1) allows additional parameter name/value pairs to be specified.  These
    %   are passed directly through to FlyDiscoPipeline().  For instance, a caller
    %   might set the parameter 'docomputeperframestats' to 'off' to disable the
    %   computeperframestats stage of FlyDiscoPipeline().  See the documentation of
    %   FlyDiscoPipeline() for more details, including a complete list of the
    %   supported optional parameters.  Some other commonly-used parameters are
    %   'settingsdir', to specify the settings directory to use; and
    %   'analysis_protocol', to specify the subdirectory within the settings
    %   directory to use (overriding the screen_type set in the experiment
    %   folders).

    % Process arguments                                
    if ~exist('do_use_bqueue', 'var') || isempty(do_use_bqueue) ,
        do_use_bqueue = true ;
    end
    if ~exist('do_actually_submit_jobs', 'var') || isempty(do_actually_submit_jobs) ,
        do_actually_submit_jobs = true ;
    end
    if ~exist('do_try', 'var') || isempty(do_try) ,
        do_try = true ;
    end
    if ~exist('ssh_host_name', 'var') || isempty(ssh_host_name) ,
        ssh_host_name = '' ;
    end

    % Extract any bsub arguments passed in
    [slots_per_job, maximum_slot_count, do_use_xvfb, optional_pipeline_arguments_as_name_value_list] = ...
      myparse_nocheck(varargin, ...
                      'slots_per_job', 4, ...
                      'maximum_slot_count', 400, ...
                      'do_use_xvfb', true) ;
      % About do_use_xvfb: Matlab on linux leaks memory when you call getframe()
      % without an X11 server, so we typically wrap the matlab call in the 'xvfb'
      % command to give it an in-memory-only X11 server.

    % We get passed ssh_host_name as a normal arg,
    % But we want to pass that through to transfero_FlyDiscoPipeline_wrapper()
    % so that it's used for the submission of the APT part of each pipeline run.
    if ~isempty(ssh_host_name) ,
        optional_pipeline_arguments_as_name_value_list(end+1:end+2) = {'sshhost', ssh_host_name} ;
    end
    optional_pipeline_arguments_as_name_value_list(end+1:end+2) = {'do_try', do_try} ;

    % if cluster_billing_account_name is empty, look it up in the configuration
    % file
    if isempty(cluster_billing_account_name) ,
      configuration_function_name = sprintf('%s_configuration', user_name_for_configuration_purposes) ;
      configuration = feval(configuration_function_name) ;
      % Unpack the per-lab configuration file
      cluster_billing_account_name = configuration.cluster_billing_account_name ;  
    end

    % Report how many experiments are to be analyzed
    experiment_count = length(folder_path_from_experiment_index) ;
    fprintf('There are %d experiments that will be analyzed.\n', experiment_count) ;
    if experiment_count > 0 ,
        fprintf('Submitting these for analysis...\n') ;
    end

    % Run transfero_FlyDiscoPipeline_wrapper() on all experiments
    if do_use_bqueue ,
        bqueue = bqueue_type(do_actually_submit_jobs, do_try, maximum_slot_count, do_use_xvfb, ssh_host_name) ;

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
                               optional_pipeline_arguments_as_name_value_list{:}) ;
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
    else
        % If not using bqueue, just run them normally (usually just for debugging)
        job_statuses = nan(1, experiment_count) ;
        for i = 1 : experiment_count ,
            experiment_folder_path = folder_path_from_experiment_index{i} ;
            transfero_FlyDiscoPipeline_wrapper(experiment_folder_path, user_name_for_configuration_purposes, optional_pipeline_arguments_as_name_value_list{:}) ;
            job_statuses(i) = +1 ;  % Indicates completed sucessfully
        end
    end
end

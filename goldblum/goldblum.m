function goldblum(do_transfer_data_from_rigs, do_run_analysis, do_use_bqueue, do_actually_submit_jobs, analysis_parameters, configuration)
    %GOLDBLUM Transfer FlyDisco experiment folders from rig computers and analyze them.
    %   goldblum() transfers FlyDisco experiment folders from the specified rig
    %   computers and analyzes them on an LSF cluster.  What rig computers to
    %   searched, and a variety of other settings, are determined from the username of
    %   the user running goldblum().    
    
    % Deal with arguments
    if ~exist('do_transfer_data_from_rigs', 'var') || isempty(do_transfer_data_from_rigs) ,
        do_transfer_data_from_rigs = true ;
    end
    if ~exist('do_run_analysis', 'var') || isempty(do_run_analysis) ,
        do_run_analysis = true ;
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
    if ~exist('configuration', 'var') || isempty(configuration) ,
        % Load the per-lab configuration file
        user_name = get_user_name() ;
        configuration_function_name = sprintf('%s_configuration', user_name) ;
        configuration = feval(configuration_function_name) ;
    end
    
    % Unpack the per-lab configuration file
    cluster_billing_account_name = configuration.cluster_billing_account_name ;
    host_name_from_rig_index = configuration.host_name_from_rig_index ;
    rig_user_name_from_rig_index = configuration.rig_user_name_from_rig_index ;
    data_folder_path_from_rig_index = configuration.data_folder_path_from_rig_index ;
    destination_folder = configuration.destination_folder ;    
    settings_folder_path = configuration.settings_folder_path ;
    %does_use_per_user_folders = configuration.does_use_per_user_folders ;
    to_process_folder_name = 'to-process' ;

    % Add a "banner" to the start of the log
    start_time_as_char = char(datetime('now','TimeZone','local','Format','y-MM-dd HH:mm Z')) ;
    fprintf('\n') ;
    fprintf('********************************************************************************\n') ;
    fprintf('\n') ;
    fprintf('Goldblum run starting at %s\n', start_time_as_char) ;
    fprintf('\n') ;
    fprintf('********************************************************************************\n') ;
    fprintf('\n') ;    

    % Report the Matlab version
    matlab_ver_string = version() ;
    fprintf('Matlab version:\n%s\n\n', matlab_ver_string) ;
    
    % Get info about the state of the repo, output to log
    fprintf('FlyDiscoAnalysis repository state:\n')
    this_script_path = mfilename('fullpath') ;
    source_folder_path = realpath(fileparts(fileparts(this_script_path))) ;
    git_report = get_git_report(source_folder_path) ;
    fprintf('%s', git_report) ;

    % If the settings folder is not part of the FDA repo, print a report about it
    internal_settings_folder_path = fullfile(source_folder_path, 'settings-internal') ;
    canonical_settings_folder_path = realpath(settings_folder_path) ;
    fprintf('Settings folder repository state:\n')
    if strcmp(canonical_settings_folder_path, internal_settings_folder_path) ,
        fprintf('(Using internal settings folder)\n\n\n')
    else
        settings_git_report = get_git_report(settings_folder_path) ;
        fprintf('%s', settings_git_report) ;        
    end
    
%     % Convert e.g. flybowl-ww1.hhmi.org to flybowl-ww1    
%     short_host_name_from_rig_index = cellfun(@short_host_name_from_host_name, host_name_from_rig_index, 'UniformOutput', false) ;
%     
%     % Destination folder is different for each rig, to avoid name collisions
%     destination_folder_from_rig_index = ...
%         cellfun(@(short_host_name)(fullfile(destination_folder, short_host_name)), ...
%                 short_host_name_from_rig_index, ...
%                 'UniformOutput', false) ;
    
    % % Get the full path to the Python script that copies data from the rig machines
    % this_script_path = mfilename('fullpath') ;
    % this_folder_path = fileparts(this_script_path) ;
    % sync_script_path = fullfile(this_folder_path, 'remote_sync_and_delete_contents.py') ;
    
    % For each rig, copy the data over to the Janelia filesystem, and delete the
    % original data
    if do_transfer_data_from_rigs ,
        rig_count = length(host_name_from_rig_index) ;
        for rig_index = 1 : rig_count ,
            rig_host_name = host_name_from_rig_index{rig_index} ;
            rig_user_name = rig_user_name_from_rig_index{rig_index} ;
            lab_data_folder_path = data_folder_path_from_rig_index{rig_index} ;

            try
                relative_path_from_synched_experiment_folder_index = ...
                    remote_sync_verify_and_delete_experiment_folders(rig_user_name, ...
                                                                     rig_host_name, ...
                                                                     lab_data_folder_path, ...
                                                                     destination_folder, ...
                                                                     to_process_folder_name) ;                
                add_links_to_to_process_folder(destination_folder, to_process_folder_name, relative_path_from_synched_experiment_folder_index) ;
            catch me 
                fprintf('There was a problem doing the sync from %s:%s as %s to %s:\n', ...
                        rig_host_name, lab_data_folder_path, rig_user_name, destination_folder) ;
                disp(me.getReport()) ;    
            end                        
        end
    else
        fprintf('Skipping transfer of data from rigs.\n') ;
    end
    
    % Run the analysis script on links in the to-process folder
    if do_run_analysis ,
        % Get the links from the to_process_folder_name folder
        to_process_folder_path = fullfile(destination_folder, to_process_folder_name) ;
        folder_name_from_experiment_index = simple_dir(to_process_folder_path) ;
        link_path_from_experiment_index = ...
            cellfun(@(folder_name)(fullfile(to_process_folder_path, folder_name)), ...
                                   folder_name_from_experiment_index, ...
                                   'UniformOutput', false) ;
        canonical_path_from_experiment_index = cellfun(@realpath, link_path_from_experiment_index, 'UniformOutput', false) ;
        goldblum_analyze_experiment_folders(canonical_path_from_experiment_index, settings_folder_path, cluster_billing_account_name, ...
                                            do_use_bqueue, do_actually_submit_jobs, analysis_parameters)
        
        % Whether those succeeded or failed, remove the links from the
        % to-process folder
        experiment_folder_count = length(link_path_from_experiment_index) ;
        for i = 1 : experiment_folder_count ,
            experiment_folder_link_path = link_path_from_experiment_index{i} ; 
            % experiment_folder_link_path is almost certainly a symlink, but check
            % anyway
            if is_symbolic_link(experiment_folder_link_path) ,
                delete(experiment_folder_link_path) ;
            end
        end
    else
        fprintf('Skipping analysis.\n') ;
    end       
    
    % Want the start and end of a single goldblum run to be clear in the log
    fprintf('\n') ;
    fprintf('********************************************************************************\n') ;
    fprintf('\n') ;
    fprintf('Goldblum run started at %s is ending\n', start_time_as_char) ;
    fprintf('\n') ;
    fprintf('********************************************************************************\n') ;
    fprintf('\n') ;    
end

function goldblum(do_transfer_data_from_rigs, do_run_analysis, do_use_bqueue, do_actually_submit_jobs, analysis_parameters, configuration)
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
        index_of_ell_in_lab = regexp(user_name, 'lab$', 'once') ;
        if isempty(index_of_ell_in_lab) ,
            error('No per-lab configuration specified, and username "%s" does not seem to be a shared lab user', user_name) ;
        end
        pi_last_name = user_name(1:index_of_ell_in_lab-1) ;
        configuration_function_name = sprintf('%s_configuration', pi_last_name) ;
        configuration = feval(configuration_function_name) ;
    end
    
    % Unpack the per-lab configuration file
    lab_head_last_name = configuration.lab_head_last_name ;
    host_name_from_rig_index = configuration.host_name_from_rig_index ;
    rig_user_name_from_rig_index = configuration.rig_user_name_from_rig_index ;
    data_folder_path_from_rig_index = configuration.data_folder_path_from_rig_index ;
    destination_folder = configuration.destination_folder ;    
    settings_folder_path = configuration.settings_folder_path ;
    does_use_per_user_folders = configuration.does_use_per_user_folders ;
    
%     % Convert e.g. flybowl-ww1.hhmi.org to flybowl-ww1    
%     short_host_name_from_rig_index = cellfun(@short_host_name_from_host_name, host_name_from_rig_index, 'UniformOutput', false) ;
%     
%     % Destination folder is different for each rig, to avoid name collisions
%     destination_folder_from_rig_index = ...
%         cellfun(@(short_host_name)(fullfile(destination_folder, short_host_name)), ...
%                 short_host_name_from_rig_index, ...
%                 'UniformOutput', false) ;
    
    % The data for. e.g. the Branson lab, will be in <rig_data_folder_path>/branson
    lab_data_folder_path_from_rig_index = ...
        cellfun(@(data_folder_path)(fullfile(data_folder_path, lab_head_last_name)), ...
                data_folder_path_from_rig_index, ...
                'UniformOutput', false) ;
    
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
            lab_data_folder_path = lab_data_folder_path_from_rig_index{rig_index} ;

            try
                remote_sync_and_verify_and_delete_contents(rig_user_name, rig_host_name, lab_data_folder_path, destination_folder, does_use_per_user_folders) ;
            catch me 
                fprintf('There was a problem doing the sync from %s:%s as %s to %s:\n', ...
                        rig_host_name, lab_data_folder_path, rig_user_name, destination_folder) ;
                disp(me.getReport()) ;    
            end
        end
    else
        fprintf('Skipping transfer of data from rigs.\n') ;
    end
    
    % Run the analysis script on the destination folder
    if do_run_analysis ,
        find_experiments_and_analyze(destination_folder, settings_folder_path, lab_head_last_name, ...
                                     do_use_bqueue, do_actually_submit_jobs, analysis_parameters) ;
    else
        fprintf('Skipping analysis.\n') ;
    end        
end

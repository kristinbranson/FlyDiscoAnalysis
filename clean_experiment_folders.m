function clean_experiment_folders(root_folder_name, do_delete_tracking_files, do_wet_run)
    % Delete the pipeline output files from any experiment folders under root_folder_name.
    % If do_delete_tracking_files is true, delete the tracking files also.
    % If do_delete_tracking_files is missing or empty, default is to *not* delete
    % tracking files.
    % If do_wet_run is true, actually do the deleting.  Otherwise, just print what
    % files would have been deleted.
    % If do_wet_run is missing or empty, default is to *not* do the deleting.
    %
    % PLEASE DON'T USE THIS ON ANYTHING BUT TEST EXPERIMENT FOLDERS.
    % IT DELETES ANY FILE WHOSE NAME DOES NOT MATCH A WHITELIST.
    
    if ~exist('do_delete_tracking_files', 'var') || isempty(do_delete_tracking_files) ,
        do_delete_tracking_files = false ;
    end
    if ~exist('do_wet_run', 'var') || isempty(do_wet_run) ,
        do_wet_run = false ;
    end
    
    folder_path_from_experiment_index = find_experiment_folders(root_folder_name) ;
    cellfun(@(experiment_folder_name)(clean_experiment_folder(experiment_folder_name, do_delete_tracking_files, do_wet_run)), ...
            folder_path_from_experiment_index) ;
end

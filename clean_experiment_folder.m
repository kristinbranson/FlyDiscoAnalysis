function clean_experiment_folder(experiment_folder_name, do_delete_tracking_files, do_wet_run)
    % Delete the pipeline output files from an experiment folder.
    % If do_delete_tracking_files is true, delete the tracking files also.
    % If do_delete_tracking_files is missing or empty, default is to *not* delete
    % tracking files.
    % PLEASE DON'T USE THIS ON ANYTHING BUT TEST EXPERIMENT FOLDERS.
    % IT DELETES ANY FILE WHOSE NAME DOES NOT MATCH A WHITELIST.
    
    if ~exist('do_delete_tracking_files', 'var') || isempty(do_delete_tracking_files) ,
        do_delete_tracking_files = false ;
    end
    if ~exist('do_wet_run', 'var') || isempty(do_wet_run) ,
        do_wet_run = false ;
    end
    
    file_name_regexp_from_whitelist_index = ...
        {'^cameraSettings.*\.json$', ...
         '^expFigure\.jpg$', ...
         '^expTimeStamp\.txt$', ...
         '^metaData\.xml$', ...
         '^metaData\.xml\.bak$', ...
         '^expNotes\.xml$', ...
         '^movie\.ufmf$', ...
         '^protocol\.mat$', ...
         '^protocol\.xlsx$', ...
         '^protocol\.csv$', ...
         '^flytracker-calibration.mat$', ...
         '^stamp_log_(.*)\.txt$', ...
         '^.*BARCODE.*\.csv$', ...} ;
         '^FlyBowlDataCaptureParams.*\.txt$', ...
         '^Log\.txt$', ...
         '^StimulusTimingLog\.txt$', ...
         '^SUCCESS$', ...
         '^QuickStats\.png$', ...
         '^QuickStats\.txt$'} ;
   
    if ~do_delete_tracking_files ,
        file_name_regexp_from_whitelist_index = [ file_name_regexp_from_whitelist_index, {'^movie-bg\.mat$', '^movie-feat\.mat$', '^movie-track\.mat$'} ] ;
    end
    
    leaf_name_from_file_index = simple_dir(experiment_folder_name) ;
    does_match_whitelist_from_file_index = does_match_some_regexpi(leaf_name_from_file_index, file_name_regexp_from_whitelist_index) ;
    do_delete_from_file_index = ~does_match_whitelist_from_file_index ;
    leaf_file_name_from_to_delete_index = leaf_name_from_file_index(do_delete_from_file_index) ;
    for to_delete_index = 1 : length(leaf_file_name_from_to_delete_index) ,
        leaf_file_name = leaf_file_name_from_to_delete_index{to_delete_index} ;
        file_name = fullfile(experiment_folder_name, leaf_file_name) ;
        if do_wet_run, 
            system_from_list_with_error_handling({'rm', '-rf', file_name}) ;
        else
            fprintf('rm -rf %s\n', file_name) ;
        end
    end
end



function result = does_match_regexpi(strs, re)
    result = cellfun(@(x)(~isempty(x)), regexpi(strs, re, 'once')) ;
end



function result = does_match_some_regexpi(strs, res)
    result = false(size(strs)) ;
    for i = 1 : length(res) ,
        re = res{i} ;        
        result = result | does_match_regexpi(strs, re) ;
    end
end

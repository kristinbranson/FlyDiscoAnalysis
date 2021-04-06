function [n_files_verified, n_dirs_verified, n_file_bytes_verified] = ...
        remote_verify_helper(source_user, source_host, source_parent_path, dest_parent_path, n_files_verified, n_dirs_verified, n_file_bytes_verified, spinner) 
    % get a list of all files and dirs in the source, dest dirs
    [name_from_source_entry_index, ...
     size_from_source_entry_index, ...
     is_file_from_source_entry_index, ...
     is_folder_from_source_entry_index, ...
     ~, ...
     mod_datetime_from_source_entry_index] = ...
        list_remote_dir(source_user, source_host, source_parent_path) ;
    [name_from_dest_entry_index, is_folder_from_dest_entry_index, size_from_dest_entry_index, mod_datetime_from_dest_entry_index] = ...
        simple_dir(dest_parent_path) ;
    
    % separate source file, dir names (ignore links)
    name_from_source_file_index = name_from_source_entry_index(is_file_from_source_entry_index) ;
    size_from_source_file_index = size_from_source_entry_index(is_file_from_source_entry_index) ;
    mod_datetime_from_source_file_index = mod_datetime_from_source_entry_index(is_file_from_source_entry_index) ;
    name_from_source_folder_index = name_from_source_entry_index(is_folder_from_source_entry_index) ;
        
    % separate dest file, dir names
    name_from_dest_folder_index = name_from_dest_entry_index(is_folder_from_dest_entry_index) ;    
    is_file_from_dest_entry_index = ~is_folder_from_dest_entry_index ;
    name_from_dest_file_index = name_from_dest_entry_index(is_file_from_dest_entry_index) ;    
    size_from_dest_file_index = size_from_dest_entry_index(is_file_from_dest_entry_index) ;
    mod_datetime_from_dest_file_index = mod_datetime_from_dest_entry_index(is_file_from_dest_entry_index) ;
    
    % verify that all files in source are also in dest
    are_all_source_files_in_dest = isempty(setdiff(name_from_source_file_index, name_from_dest_file_index)) ;
    if ~are_all_source_files_in_dest ,
        error('The local files in %s do not include all the remote files in %s', dest_parent_path, source_parent_path) ;
    end
    
    % scan the source files, compute the md5sum, compare to that for the dest file
    for source_file_index = 1 : length(name_from_source_file_index) ,
        file_name = name_from_source_file_index{source_file_index} ;
        source_file_path = fullfile(source_parent_path, file_name) ;
        dest_file_path = fullfile(dest_parent_path, file_name) ;        
        source_file_size = size_from_source_file_index(source_file_index) ;
        dest_file_index = find(strcmp(file_name, name_from_dest_file_index)) ;        
        dest_file_size = size_from_dest_file_index(dest_file_index) ;
        spinner.spin() ;
        if dest_file_size == source_file_size ,
            source_mod_datetime = mod_datetime_from_source_file_index(source_file_index) ;
            dest_mod_datetime = mod_datetime_from_dest_file_index(dest_file_index) ;
            if ( source_mod_datetime < dest_mod_datetime ) ,
                % Compare hashes
                source_hex_digest = compute_md5_on_remote(source_user, source_host, source_file_path) ;
                dest_hex_digest = compute_md5_on_local(dest_file_path) ;
                is_file_verified = isequal(source_hex_digest, dest_hex_digest) ;
            else 
                is_file_verified = false ;
            end
        else 
            is_file_verified = false ;
        end
        
        if is_file_verified ,
            n_files_verified = n_files_verified + 1 ;
            n_file_bytes_verified = n_file_bytes_verified + source_file_size ;
        else
            error("There is a problem with destination file %s: It's missing, or the wrong size, or the hashes don't match.", ...
                  dest_file_path) ;
        end
    end

    % Verify that all source folders are in destination
    are_all_source_files_in_dest = isempty(setdiff(name_from_source_folder_index, name_from_dest_folder_index)) ;
    if are_all_source_files_in_dest ,
        n_dirs_verified = n_dirs_verified + length(name_from_source_folder_index) ;
    else
        error('The local folder names in %s do not include all the remote folder names in %s', dest_parent_path, source_parent_path) ;
    end
    
    % for each source folder, recurse
    for i = 1 : length(name_from_source_folder_index) ,
        folder_name = name_from_source_folder_index{i} ;
        source_folder_path = fullfile(source_parent_path, folder_name) ;
        dest_folder_path = fullfile(dest_parent_path, folder_name) ;                
        [n_files_verified, n_dirs_verified, n_file_bytes_verified] = ...
             remote_verify_helper(source_user, ...
                                  source_host, ...
                                  source_folder_path, ...
                                  dest_folder_path, ...
                                  n_files_verified, ...
                                  n_dirs_verified, ... 
                                  n_file_bytes_verified, ...
                                  spinner) ;
    end
end

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
    
    % scan the source files, make sure they're all present in dest, with matching
    % hashes
    for source_file_index = 1 : length(name_from_source_file_index) ,
        file_name = name_from_source_file_index{source_file_index} ;
        source_file_path = fullfile(source_parent_path, file_name) ;
        dest_file_path = fullfile(dest_parent_path, file_name) ;        
        source_file_size = size_from_source_file_index(source_file_index) ;
        dest_file_index = find(strcmp(file_name, name_from_dest_file_index)) ;        
        if isempty(dest_file_index) ,
            error("There is a problem with destination file %s: It's missing.", dest_file_path) ;
        elseif isscalar(dest_file_index) ,
            dest_file_size = size_from_dest_file_index(dest_file_index) ;
            spinner.spin() ;
            if dest_file_size == source_file_size ,
                source_mod_datetime = mod_datetime_from_source_file_index(source_file_index) ;
                dest_mod_datetime = mod_datetime_from_dest_file_index(dest_file_index) ;
                if ( source_mod_datetime < dest_mod_datetime ) ,
                    % Compare hashes
                    source_hex_digest = compute_md5_on_remote(source_user, source_host, source_file_path) ;
                    dest_hex_digest = compute_md5_on_local(dest_file_path) ;
                    if ~isequal(source_hex_digest, dest_hex_digest) ,
                        error("There is a problem with destination file %s: Its hash is %s, but the source file hash is %s.", ...
                            dest_file_path, dest_hex_digest, source_hex_digest) ;
                    end
                else
                    error("There is a problem with destination file %s: Its modification time (%s )is before that of the source file (%s).", ...
                        dest_file_path, char(dest_mod_datetime), char(source_mod_datetime) ) ;
                end
            else
                error("There is a problem with destination file %s: Its size (%d bytes) differs from that of the source file (%d bytes).", ...
                    dest_file_path, dest_file_size, source_file_size) ;
            end
        else
            error('Something has gone horribly wrong in remote_verify_helper().  There seem to be two files with the same name (%s) in destination folder %s', ...
                  file_name, dest_parent_path) ;
        end
        
        % If we get here, destination file is present, it was modified after the source,
        % and the hashes match
        n_files_verified = n_files_verified + 1 ;
        n_file_bytes_verified = n_file_bytes_verified + source_file_size ;
    end

    % Verify that all source folders are in destination
    for source_folder_index = 1 : length(name_from_source_folder_index) ,
        folder_name = name_from_source_folder_index{source_folder_index} ;
        if ~any(strcmp(folder_name, name_from_dest_folder_index)) ,
            dest_folder_path = fullfile(dest_parent_path, folder_name) ;
            error("There is a problem with destination folder %s: It's missing.", dest_folder_path) ;
        end
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

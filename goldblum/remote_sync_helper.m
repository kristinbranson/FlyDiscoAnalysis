function [n_copied, n_failed, n_dir_failed, n_verified, n_dir_failed_to_list, time_spent_copying] = ...
        remote_sync_helper(source_user, ...
                           source_host, ...
                           source_parent_path, ...
                           dest_parent_path, ...
                           n_copied, ...
                           n_failed, ...
                           n_dir_failed, ...
                           n_verified, ...
                           n_dir_failed_to_list, ...
                           time_spent_copying, ...
                           spinner) 
    % get a list of all files and dirs in the source, dest dirs
    try
        [source_entries, source_entry_sizes, is_source_entry_a_file, is_source_entry_a_dir, ~, source_entry_mod_times] = ...
            list_remote_dir(source_user, source_host, source_parent_path) ;
    catch me ,
        % if we can't list the dir, warn but continue
        if isequal(me.identifier, 'list_remote_dir:failed') ,
            spinner.print("Warning: can't list path %s on host %s as user %s", source_path, source_host, source_user) ;
            spinner.print("%s", me.getReport()) ;
            n_dir_failed_to_list = n_dir_failed_to_list + 1 ;
            return
        else
            rethrow(me) ;
        end
    end
    [dest_entries, is_dest_entry_a_folder, dest_file_size, dest_file_mtimes] = simple_dir(dest_parent_path) ;
    
    % separate source file, dir names (ignore links)
    source_file_names = source_entries(is_source_entry_a_file) ;
    source_file_sizes = source_entry_sizes(is_source_entry_a_file) ;
    source_file_mod_times = source_entry_mod_times(is_source_entry_a_file) ;
    source_folder_names = source_entries(is_source_entry_a_dir) ;
        
    % separate dest file, dir names
    dest_folder_names = dest_entries(is_dest_entry_a_folder) ;    
    dest_file_names = dest_entries(~is_dest_entry_a_folder) ;    
    size_from_dest_file_index = dest_file_size(~is_dest_entry_a_folder) ;
    dest_file_mtimes = dest_file_mtimes(~is_dest_entry_a_folder) ;
    
    % scan the source files, copy any that aren't in dest,
    % or that aren't up-to-date
    source_file_count = length(source_file_names) ;
    for i = 1 : source_file_count ,
        file_name = source_file_names{i} ;
        source_file_path = fullfile(source_parent_path, file_name) ;
        dest_file_path = fullfile(dest_parent_path, file_name) ;
        spinner.spin() ;
        %print("  %s" % source_file)
        dest_file_index = find(strcmp(file_name, dest_file_names)) ;        
        if ~isempty(dest_file_index) ,
            source_file_size = source_file_sizes(i) ;
            dest_file_size = size_from_dest_file_index(dest_file_index) ;
            if dest_file_size == source_file_size ,
                source_mod_time = source_file_mod_times(i) ;
                dest_mod_time = dest_file_mtimes(i) ;
                if source_mod_time < dest_mod_time ,
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
        else
            is_file_verified = false ;
        end

        if is_file_verified ,
            n_verified = n_verified + 1 ;
        else 
            try 
                time_spent_copying = time_spent_copying + copy_file_from_remote(source_user, ...
                                                                                source_host, ... 
                                                                                source_file_path, ...
                                                                                dest_file_path) ;
                n_copied = n_copied + 1 ;
            catch me ,
                % orginally catch IOError
                % if we can't copy the file, warn but continue
                if isequal(me.identifier, 'copy_file_from_remote:failed') ,
                    spinner.print('Warning: can''t copy %s',  source_file_path) ;
                    spinner.print("%s", me.getReport()) ;
                    n_failed = n_failed + 1 ;
                else
                    rethrow(me) ;
                end
            end
        end
    end

%     % scan dest dirs, delete any that that aren't in source dirs
%     dest_dirs_as_set = set(dest_dirs)
%     source_dirs_as_set = set(source_dirs)
%     dest_dirs_not_in_source = dest_dirs_as_set - source_dirs_as_set 
%     for dest_dir in dest_dirs_not_in_source :
%         shutil.rmtree(os.path.join(dest_path,dest_dir),1)
    
    % scan source dirs, create any that that aren't in dest dirs
    source_folder_names_not_in_dest = setdiff(source_folder_names, dest_folder_names) ;
    for i = 1 : length(source_folder_names_not_in_dest) ,
        source_folder_name = source_folder_names_not_in_dest{i} ;
        dest_folder_path = fullfile(dest_parent_path, source_folder_name) ;
        [did_succeed, message] = mkdir(dest_folder_path) ;
        if ~did_succeed , 
            % if we can't make the dir, warn but continue
            spinner.print("Warning: can't make directory %s, error message was: ", dest_folder_path, message) ;
            n_dir_failed = n_dir_failed + 1 ;
        end
    end
    
    % need to re-generate dest dirs, because we may have failed to 
    % create some of them
    [dest_entries, is_dest_entry_a_folder] = simple_dir(dest_parent_path) ;
    dest_folder_names = dest_entries(is_dest_entry_a_folder) ;
    folder_names_to_recurse_into = intersect(source_folder_names, dest_folder_names) ;
        
    % for each dir in both source_dirs and dest_folder_names, recurse
    for i = 1 : length(folder_names_to_recurse_into) ,
        folder_name = folder_names_to_recurse_into{i} ;
        [n_copied, n_failed, n_dir_failed, n_verified, n_dir_failed_to_list, time_spent_copying] = ...
             remote_sync_helper(source_user, ...
                                source_host, ...
                                fullfile(source_parent_path, folder_name), ...
                                fullfile(dest_parent_path, folder_name), ...
                                n_copied, ...
                                n_failed, ...
                                n_dir_failed, ...
                                n_verified, ...
                                n_dir_failed_to_list,  ...
                                time_spent_copying, ...
                                spinner) ;
    end
end

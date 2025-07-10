function check_for_shadowing_of_builtin_functions(library_root_folder_path)
    %library_root_folder_path = fileparts(mfilename('fullpath')) ;
    path_from_m_file_index = find_and_list(library_root_folder_path, @is_m_file) ;

    cellfun(@compute_does_shadow_a_builtin, path_from_m_file_index) ;
end

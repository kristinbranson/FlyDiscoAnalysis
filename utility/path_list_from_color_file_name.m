function color_list_file_paths = path_list_from_color_file_name(color_list_glob_file_name, analysis_protocol_folder_path, expdir)
% Convert a single glob (e.g. 'perframe/apt*.mat') that is relative to the expdir to a list of
% matching (absolute) file paths.  The "color" is to communicate that this
% works for a whitelist glob file or a blacklist glob file.

% Convert the file name of the color list file into an absolute path
color_list_glob_file_path = fullfile(analysis_protocol_folder_path, color_list_glob_file_name) ;

% Read the contents of the file into a cellstr
color_list_name_globs = read_file_into_cellstring(color_list_glob_file_path) ;

% Convert the expdir-relative globs into absolute globs
color_list_path_globs = cellfun(@(name_glob)(fullfile(expdir, name_glob)), color_list_name_globs, 'UniformOutput', false) ;

% Get the file of absolute file paths matching any of the globs
color_list_file_paths = file_paths_from_globs(color_list_path_globs) ;
end



function color_list_file_paths = file_paths_from_globs(color_list_path_globs)
% Given a list of absolute globs, return all the absolute file paths that
% match at least one of the globs.  The result will not have repeats in it.
list_of_lists = cellfun(@simple_dir_with_full_path, color_list_path_globs, 'UniformOutput', false) ;
color_list_file_paths = unique(flatten_row_cell_array(list_of_lists)) ;
end



function result = simple_dir_with_full_path(template)
s_from_raw_index = dir(template) ;
name_from_raw_index = {s_from_raw_index.name} ;
is_bs_from_raw_index = ismember(name_from_raw_index, {'.', '..'}) ;
s_from_index = s_from_raw_index(~is_bs_from_raw_index) ;
result = arrayfun(@(s)(fullfile(s.folder, s.name)), s_from_index, 'UniformOutput', false) ;
end

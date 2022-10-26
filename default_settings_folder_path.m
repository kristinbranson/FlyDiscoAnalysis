function result = default_settings_folder_path()
    % Get the path to the default settings folder, based on the location of this
    % script file.  Nowadays this is normally a link to an external settings folder.
    this_script_path = mfilename('fullpath') ;
    fda_path = fileparts(this_script_path) ;
    result = fullfile(fda_path, 'settings') ;
end

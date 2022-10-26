function result = internal_settings_folder_path()
    % Get the path to the settings folder path that is 'included' in the source repo,
    % based on the location of this script file.  Make sure the folder exists before
    % returning.
    % Nowadays this is normally a link to an external settings folder, so it's not
    % really internal.
    this_script_path = mfilename('fullpath') ;
    fda_path = fileparts(this_script_path) ;
    result = fullfile(fda_path, 'settings') ;
    % This is often overidden in client code, so don't check if it exists
end

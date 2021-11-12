function result = internal_settings_folder_path()
    % Get the path to the settings folder path that is included in the source repo,
    % based on the location of this script file.  Make sure the folder exists before
    % returning.
    this_script_path = mfilename('fullpath') ;
    fda_path = fileparts(this_script_path) ;
    result = fullfile(fda_path, 'settings') ;
    if exist(result, 'dir') ,
        % all is well
    else
        error('Internal settings path %s does not exist', result) ;        
    end
end

function result = is_experiment_folder_given_contents(file_names)
    % Given a cellstring representing all the files in a folder (and excluding
    % folder names), determine whether the folder is a a FlyDisco experiment
    % folder.  Returns a logical scalar that is true iff the folder is an
    % experiment folder.

    % We check for three files.  If two or more are present, we consider it an
    % experiment folder
    lowercase_file_names = lower(file_names) ;
    has_movie_file = ( ismember('movie.ufmf', lowercase_file_names) || ismember('movie.avi', lowercase_file_names) ) ;
    point_count = ...
        double(has_movie_file) + ...
        double(ismember('metadata.xml', lowercase_file_names)) + ...
        double(ismember('ABORTED', file_names)) ;
    result = ( point_count >= 2) ;
end

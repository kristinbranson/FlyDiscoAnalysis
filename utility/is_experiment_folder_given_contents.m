function result = is_experiment_folder_given_contents(file_names)

    % specifying experiment directory, not contents
    if ischar(file_names),
      expdir = file_names;
      file_names = dir(expdir);
      file_names = {file_names.name};
    end
    lowercase_file_names = lower(file_names) ;
    % We check for three files.  If two or more are present, we consider it an
    % experiment folder
    has_movie_file = ( ismember('movie.ufmf', lowercase_file_names) || ismember('movie.avi', lowercase_file_names) ) ;
    point_count = ...
        double(has_movie_file) + ...
        double(ismember('metadata.xml', lowercase_file_names)) + ...
        double(ismember('ABORTED', file_names)) ;
    result = ( point_count >= 2) ;
end

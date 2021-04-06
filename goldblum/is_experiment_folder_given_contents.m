function result = is_experiment_folder_given_contents(file_names)
    lowercase_file_names = lower(file_names) ;
    result = ( ismember('movie.ufmf', lowercase_file_names) || ismember('movie.avi', lowercase_file_names) ) && ...
             ismember('metadata.xml', lowercase_file_names) ;
end

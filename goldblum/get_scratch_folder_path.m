function result = get_scratch_folder_path()  
  if exist('/scratch', 'dir') ,
    % If this folder exists, return it, since presumably we're on a cluster node
    result = fullfile('/scratch', get_user_name()) ;
  else
    % Otherwise, return the usual Matlab temp folder name
    result = tempdir() ;
  end
end

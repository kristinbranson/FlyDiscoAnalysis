function turn_on_goldblum(hr, min)
    % Install the goldblum job in crontab.
    % Can give optional hr, min args, which specify the time to run, in 24-hour
    % clock format.  I.e. turn_on_goldblum(23,11) sets it to run once a day at 11:11
    % PM.  Default time is 10:00 PM if no args are given.
    
    if ~exist('hr', 'var') || isempty(hr) ,
        hr = 22 ;
    else
        hr = round(hr) ;
        if hr<0 || hr>23 ,
            error('hr much be an integer between 0 and 23, inclusive') ;
        end
    end
    if ~exist('min', 'var') || isempty(min) ,
        min = 0 ;
    else
        min = round(min) ;
        if min<0 || min>23 ,
            error('min much be an integer between 0 and 59, inclusive') ;
        end
    end
    
    user_name = get_user_name() ;
    index_of_ell_in_lab = regexp(user_name, 'lab$', 'once') ;
    if isempty(index_of_ell_in_lab) ,
        error('No per-lab configuration specified, and username "%s" does not seem to be a shared lab user', user_name) ;
    end
    pi_last_name = user_name(1:index_of_ell_in_lab-1) ;
    configuration_function_name = sprintf('%s_configuration', pi_last_name) ;
    configuration = feval(configuration_function_name) ;
    
    destination_folder_path = configuration.destination_folder ;
    escaped_destination_folder_path = escape_string_for_bash(destination_folder_path) ;
    this_folder_path = fileparts(mfilename('fullpath')) ;
    fly_disco_analysis_folder_path = fileparts(this_folder_path) ;
    escaped_fly_disco_analysis_folder_path = escape_string_for_bash(fly_disco_analysis_folder_path) ;    
    
    goldblum_log_file_path = fullfile(destination_folder_path, 'goldblum.log') ;
    escaped_goldblum_log_file_path = escape_string_for_bash(goldblum_log_file_path) ;
    
    home_folder_path = getenv('HOME') ;
    bash_profile_path = fullfile(home_folder_path, '.bash_profile') ;
    escaped_bash_profile_path = escape_string_for_bash(bash_profile_path) ;
    
    core_command_line = ...
        sprintf(['. /misc/lsf/conf/profile.lsf ; ' ...
                 '. %s ; ' ...
                 'cd %s ; ' ...
                 'bsub -n1 -P %s -o %s -e %s /misc/local/matlab-2019a/bin/matlab -nodisplay -batch ''modpath; goldblum(true, true);'''], ...
                escaped_bash_profile_path, ...
                escaped_fly_disco_analysis_folder_path, ...
                pi_last_name, ...
                escaped_goldblum_log_file_path, ...
                escaped_goldblum_log_file_path)  %#ok<NOPRT>
    escaped_core_command_line = escape_string_for_bash(core_command_line) ;
    
    hash_goldblum = '#GOLDBLUM' ;
    escaped_hash_goldblum = escape_string_for_bash(hash_goldblum) ;
        
    command_line = sprintf('{ crontab -l | grep --invert-match %s; echo "%02d %02d * * *     flock --nonblock %s --command %s   #GOLDBLUM"; } | crontab', ...
                           min, ...
                           hr, ...
                           escaped_hash_goldblum, ...
                           escaped_destination_folder_path, ...
                           escaped_core_command_line)      %#ok<NOPRT> % Run at whenever every day                       
    
    % Clear out any pre-existing #GOLDBLUM crontab lines
    turn_off_goldblum()

    % Execute the command to turn on goldblum
    system_with_error_handling(command_line) ;    
end

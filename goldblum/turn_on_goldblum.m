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
            error('hr must be an integer between 0 and 23, inclusive') ;
        end
    end
    if ~exist('min', 'var') || isempty(min) ,
        min = 0 ;
    else
        min = round(min) ;
        if min<0 || min>59 ,
            error('min must be an integer between 0 and 59, inclusive') ;
        end
    end
    
    user_name = get_user_name() ;
    configuration_function_name = sprintf('%s_configuration', user_name) ;
    configuration = feval(configuration_function_name) ;
    cluster_billing_account_name = configuration.cluster_billing_account_name ;
    
    destination_folder_path = configuration.destination_folder ;
    escaped_destination_folder_path = escape_string_for_bash(destination_folder_path) ;
    this_folder_path = fileparts(mfilename('fullpath')) ;
    fly_disco_analysis_folder_path = fileparts(this_folder_path) ;
    escaped_fly_disco_analysis_folder_path = escape_string_for_bash(fly_disco_analysis_folder_path) ;    
    
    goldblum_logs_folder_path = fullfile(destination_folder_path, 'goldblum-logs') ;
    escaped_goldblum_logs_folder_path = escape_string_for_bash(goldblum_logs_folder_path) ;
    
    home_folder_path = getenv('HOME') ;
    bash_profile_path = fullfile(home_folder_path, '.bash_profile') ;
    escaped_bash_profile_path = escape_string_for_bash(bash_profile_path) ;
    
    launcher_script_path = fullfile(this_folder_path, 'goldblum_launcher.sh') ;
    escaped_launcher_script_path = escape_string_for_bash(launcher_script_path) ;
    
%     escaped_bash_profile_path=${1}
%     escaped_fly_disco_analysis_folder_path=${2}
%     pi_last_name=${3}
%     escaped_goldblum_logs_folder_path=${4}
%     date_as_string=`date +%Y-%m-%d`
%     goldblum_log_file_name="goldblum-${date_as_string}.log"
%     goldblum_log_file_path="${escaped_goldblum_logs_folder_path}/${goldblum_log_file_name}" 

    core_command_line = ...
        sprintf('%s %s %s %s %s', ...
                escaped_launcher_script_path, ...
                escaped_bash_profile_path, ...
                escaped_fly_disco_analysis_folder_path, ...
                cluster_billing_account_name, ...
                escaped_goldblum_logs_folder_path)  %#ok<NOPRT>

%     core_command_line = ...
%         sprintf(['. /misc/lsf/conf/profile.lsf ; ' ...
%                  '. %s ; ' ...
%                  'cd %s ; ' ...
%                  'bsub -n1 -P %s -o %s -e %s /misc/local/matlab-2019a/bin/matlab -nodisplay -batch ''modpath; goldblum(true, true);'''], ...
%                 escaped_bash_profile_path, ...
%                 escaped_fly_disco_analysis_folder_path, ...
%                 pi_last_name, ...
%                 escaped_goldblum_log_folder_path, ...
%                 escaped_goldblum_log_folder_path)  %#ok<NOPRT>
    escaped_core_command_line = escape_string_for_bash(core_command_line) ;
    
    hash_goldblum = '#GOLDBLUM' ;
    escaped_hash_goldblum = escape_string_for_bash(hash_goldblum) ;
        
    command_line = sprintf('{ crontab -l | grep --invert-match %s; echo "%02d %02d * * *     flock --nonblock %s --command %s   #GOLDBLUM"; } | crontab', ...
                           escaped_hash_goldblum, ...
                           min, ...
                           hr, ...
                           escaped_destination_folder_path, ...
                           escaped_core_command_line)      %#ok<NOPRT> % Run at whenever every day                       
    
    % Clear out any pre-existing #GOLDBLUM crontab lines
    turn_off_goldblum()

    % Execute the command to turn on goldblum
    system_with_error_handling(command_line) ;    
end

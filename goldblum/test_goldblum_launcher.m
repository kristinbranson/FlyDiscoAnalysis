function test_goldblum_launcher()    
    pi_last_name = 'branson' ;
    configuration_function_name = sprintf('%s_configuration', pi_last_name) ;
    configuration = feval(configuration_function_name) ;
    
    destination_folder_path = configuration.destination_folder ;
    this_folder_path = fileparts(mfilename('fullpath')) ;
    fly_disco_analysis_folder_path = fileparts(this_folder_path) ;
    
    goldblum_logs_folder_path = fullfile(destination_folder_path, 'goldblum-logs') ;
    
    home_folder_path = getenv('HOME') ;
    bash_profile_path = fullfile(home_folder_path, '.bash_profile') ;
    
    launcher_script_path = fullfile(this_folder_path, 'goldblum_launcher.sh') ;
    
    % Execute the command to turn on goldblum
    stdout = ...
        system_from_list_with_error_handling({launcher_script_path, ...
                                             bash_profile_path, ...
                                             fly_disco_analysis_folder_path, ...
                                             pi_last_name, ...
                                             goldblum_logs_folder_path}) 
end

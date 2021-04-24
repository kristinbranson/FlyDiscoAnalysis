function turn_off_goldblum()
    command_line = 'crontab -l | grep --invert-match ''#GOLDBLUM'' | crontab' ;   % Remove line containing #GOLDBLUM                       
    system_with_error_handling(command_line) ;    
end

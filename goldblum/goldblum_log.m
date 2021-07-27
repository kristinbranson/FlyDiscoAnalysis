function goldblum_log(varargin)
    % Write stuff to a log file, or just to stdout if no log file set up.   
    
    persistent log_file_object
    
    if nargin>=1 && isscalar(varargin{1}) && isnan(varargin{1}) ,
        % this means the call to goldblum_log is a "command"
        command = varargin{2} ;
        if strcmp(command, 'open') ,        
            file_name = varargin{3} ;
            log_file_object = file_object(file_name, 'wt') ;
        elseif strcmp(command, 'close') ,        
            log_file_object = [] ;  % will close the file if no other references to it
        else
            error('goldblum_log:unknown_command', 'Unknown command "%s"', command) ;            
        end
    else
        % Normal calls to this function fprintf to the log, or stdout if no log
        if isempty(log_file_object) || ~isvalid(log_file_object) ,
            % Just send to stdout
            fprintf(1, varargin{:}) ;            
        else
            log_file_object.fprintf(varargin{:}) ;
        end            
    end        
end

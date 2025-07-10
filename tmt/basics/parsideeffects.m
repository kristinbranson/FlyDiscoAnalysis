function parsideeffects(f, x, do_run_in_debug_mode, varargin)
    % Evaluates f(x(i)) for each element of x.  f() is assumed to return nothing.
    % That is, it is run strictly for its side effects.
    % If do_run_in_debug_mode is true, a (serial) for loop is used.
    % Otherwise, a parfor loop is used.
    if ~exist('do_run_in_debug_mode', 'var') || isempty(do_run_in_debug_mode) ,
        do_run_in_debug_mode = false ;
    end
    
    n = length(x) ;
    if do_run_in_debug_mode ,
        for i = 1 : n ,
            feval(f, x(i), varargin{:}) ;
        end
    else
        parfor i = 1 : n ,
            feval(f, x(i), varargin{:}) ;  %#ok<PFBNS>
        end
    end
end

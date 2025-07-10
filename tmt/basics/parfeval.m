function y = pararrayfun(f, x, y1, do_run_in_debug_mode, varargin)
    % Like arrayfun, but uses parfor.
    % If do_run_in_debug_mode is true, runs serially with a for loop
    if ~exist('y1', 'var') || isempty(y1) ,
        y1 = 0 ;
    end
    if ~exist('do_run_in_debug_mode', 'var') || isempty(do_run_in_debug_mode) ,
        do_run_in_debug_mode = false ;
    end
    
    n = length(x) ;
    y = repmat(y1, size(x)) ;
    if do_run_in_debug_mode ,
        for i = 1 : n ,
            y(i) = feval(f, x(i), varargin{:}) ;
        end
    else
        parfor i = 1 : n ,
            y(i) = feval(f, x(i), varargin{:}) ;  %#ok<PFBNS>
        end
    end
end

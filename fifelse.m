function y = fifelse(varargin)
    % Like fif(), but allows for multiple <gate, value> pairs.
    % The number of arguments must be odd, and the last argument is treated as a fallback, similar to
    % the terminal 'else' clause in an if...elseif...else block.
    
    nargin = length(varargin) ;
    if mod(nargin,2) == 0 ,
        error('The number of arguments to %s() must be odd.  The last argument is the fallback value.') ;
    else
        y = varargin{nargin} ;
    end
    for i = 1 : 2 : nargin-1 ,
        if i < nargin , 
            test = varargin{i} ;
            y_true = varargin{i+1} ;
            y(test) = y_true(test) ;
        end 
    end
end

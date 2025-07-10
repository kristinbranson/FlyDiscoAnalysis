function varargout = lapwing(varargin)
    % create the controller instance
    rc = lapwing.Controller(varargin{:}) ;

    if nargout >= 1 ,
        varargout{1} = rc ;
    end
end

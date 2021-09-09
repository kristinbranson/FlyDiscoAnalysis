classdef path_object
    properties
        list_
        is_relative_
    end
    
    methods
        function self = path_object(varargin)
            if nargin==1 ,
                path = varargin{1} ;
                if isempty(path) ,
                    error('path cannot be empty') ;
                end
                if isequal(path(1), '/') ,          
                    self.is_relative_ = false ;
                    if length(path)==1 ,
                        normed_path = '' ;
                    else
                        normed_path = path ;
                    end
                else
                    self.is_relative_ = true ;
                    normed_path = horzcat('/', path) ;
                end
                self.list_ = list_from_normed_path(normed_path) ;
            else
                self.list_ = varargin{1} ;
                self.is_relative_ = varargin{2} ;
            end
        end
        
        function result = list(self) 
            result = self.list_ ;
        end
        
        function result = is_relative(self) 
            result = self.is_relative_ ;
        end
        
        function result = is_absolute(self) 
            result = ~self.is_relative_ ;
        end        
        
        function result = to_char(self)
            normed_path = normed_path_from_list(self.list_) ;
            if self.is_relative_ ,
                result = normed_path(2:end) ;
            else
                result = normed_path ;
            end
        end

        function result = cat(self, other)
            if is_absolute(other) ,
                error('2nd argument can''t be an absolute path') ;
            end
            result = path_object(horzcat(list(self), list(other)), ...             
                                 is_relative(self)) ;
        end
        
        function result = relpath(self, start)
            if ischar(start) ,
                start_path_object = path_object(start) ;
            else
                start_path_object = start ;
            end
            if is_absolute(start_path_object) ~= is_absolute(self),
                error('start and self must both be absolute or both relative') ;
            end
            result = path_object(relpath_helper(list(self), list(start_path_object)), ...
                                 true) ;  % result is always relative            
        end
    end
end



function result = list_from_normed_path(normed_path) 
    if isequal(normed_path, '/') ,
        result = cell(1,0) ;
    else        
        [parent, name] = fileparts2(normed_path) ;
        result = horzcat(list_from_normed_path(parent), {name}) ;
    end
end



function result = normed_path_from_list(list)
    if isempty(list) ,
        result = '' ;
    else
        head = list{1} ;
        tail = list(2:end) ;
        result = horzcat('/', head, normed_path_from_list(tail)) ;
    end
end



function result = relpath_helper(list, start_list)
    if isempty(start_list) ,
        result = list ;
    else
        if isempty(list) ,
            error('start_list is not a prefix of list') ;
        else
            list_head = list{1} ;
            start_list_head = start_list{1} ;
            if isequal(list_head, start_list_head) ,
                result = relpath_helper(list(2:end), start_list(2:end)) ;
            else
                error('start_list is not a prefix of list') ;
            end
        end
    end
end

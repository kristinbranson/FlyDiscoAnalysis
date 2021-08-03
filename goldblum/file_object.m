classdef file_object < handle
    % Simple (almost stupid) file object that fcloses() the fid when no longer referenced.
    properties (Access = protected)
        fid = []
    end
    methods
        function self = file_object(name, mode) 
            [fid, msg] = fopen(name, mode) ;
            if fid>0 ,
                self.fid = fid ;
            else
                error('file_object:unable_to_open', msg) ;
            end
        end
        
        function fprintf(self, varargin)
            if ~isempty(self.fid) && self.fid>0 ,
                fprintf(self.fid, varargin{:}) ;
            else
                error('file_object:file_identifier_is_invalid', 'Can''t fprintf() because the file identifier is invalid') ;
            end
        end
        
        function delete(self)
            if ~isempty(self.fid) && self.fid>0 ,                
                fclose(self.fid) ;
                self.fid = [] ;
            end
        end
    end
end

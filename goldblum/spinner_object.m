classdef spinner_object < handle
    % Class for indicting progress
    properties (SetAccess = private)
        cursors_
        cursor_count_
        cursor_index_
        is_first_call_
        is_mute_
    end
    
    methods
        function self = spinner_object(varargin)
            self.cursors_ = '|/-\\' ;
            self.cursor_count_ = length(self.cursors_) ;
            self.cursor_index_ = 1 ;
            self.is_first_call_ = true ;
            if nargin>=1 && strcmp(varargin{1}, 'mute') ,
                self.is_mute_ = true ;
            else
                self.is_mute_ = false ;
            end
        end

        function spin(self)
            if ~self.is_mute_ ,
                if self.is_first_call_ ,
                    cursor = self.cursors_(self.cursor_index_) ;
                    fprintf('%s', cursor) ;
                    self.is_first_call_ = false ;
                else
                    fprintf('\b')
                    self.cursor_index_ = (self.cursor_index_ + 1) ;
                    if self.cursor_index_ > 4 ,
                        self.cursor_index_ = 1 ;
                    end
                    cursor = self.cursors_(self.cursor_index_) ;
                    fprintf('%s', cursor)
                end
            end
        end

        function print(self, varargin)
            % Want things printed during spinning to look nice
            if ~self.is_mute ,
                fprintf('\b\n') ;  % Delete cursor, then newline
                fprintf(varargin{:}) ;  % print whatever
                cursor = self.cursors_(self.cursor_index_) ;  % get the same cursor back
                fprintf('%s', cursor) ;  % write it again on its own line
            end
        end    

        function stop(self)
            if ~self.is_mute_ ,
                fprintf('\bdone.\n')
            end
        end
    end
end

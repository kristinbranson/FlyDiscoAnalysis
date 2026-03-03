classdef Stack_file_in_memory < handle
    
    properties (SetAccess = private)
        stack_
        data_type_
        n_row_
        n_col_
        n_frame_     
    end  % properties
        
    methods
        function self = Stack_file_in_memory(stack)
            self.stack_ = stack ;
            self.data_type_ = class(stack) ;
            [self.n_row_, self.n_col_, self.n_frame_] = size(stack) ;
        end  % constructor method

        function delete(self)  %#ok<INUSD>
        end
        
        function result = data_type(self)
            result = self.data_type_ ;
        end
        
        function result = n_row(self)
            result =self.n_row_ ;
        end
        
        function result = n_col(self)
            result = self.n_col_ ;
        end
        
        function result = n_frame(self)
            result = self.n_frame_ ;
        end

        function frame = get_frame(self, i_frame)
            frame = self.stack_(:,:,i_frame) ;
        end
        
        function frame = get_next(self, i_last_frame_read)
            % Gets frame *after* i_last_frame_read
            i_frame = i_last_frame_read + 1 ;
            frame = self.get_frame(i_frame) ;
        end        
        
        function to_start(self)  %#ok<MANU>
            % nothing to do
        end
        
    end  % methods
    
end  % classdef

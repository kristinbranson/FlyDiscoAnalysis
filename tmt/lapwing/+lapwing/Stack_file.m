classdef Stack_file < handle
    
    properties
        file_type_  % we never read this currently, but probably nice to have around
        typed_file_        
        i_frame_  % last frame gotten, 1-based index
    end  % properties
        
    methods
        function self = Stack_file(file_name_or_stack)
            if ischar(file_name_or_stack) ,
                file_name = file_name_or_stack ;
                [~,~,ext] = fileparts(file_name) ;
                switch ext ,
                    case '.tif' ,
                        is_file_an_imagej_jumbo_tif = lapwing.is_file_an_imagej_jumbo_tif(file_name) ;
                        if  is_file_an_imagej_jumbo_tif ,
                            self.file_type_ = 'imagej_jumbo_tif' ;                        
                            self.typed_file_ = lapwing.Stack_file_imagej_jumbo_tif(file_name) ;
                        else
                            self.file_type_ = 'tif' ;
                            self.typed_file_ = lapwing.Stack_file_normal_tif(file_name) ;
                        end
                    case '.mj2' ,
                        self.file_type_ = 'mj2' ;
                        self.typed_file_ = lapwing.Stack_file_mj2(file_name) ;
                    case '.h5' ,
                        self.file_type_ = 'h5' ;
                        self.typed_file_ = lapwing.Stack_file_h5(file_name) ;
                    otherwise
                        error('Stack_file:UnableToLoad','Unable to load that file type');
                end
            elseif isnumeric(file_name_or_stack) || islogical(file_name_or_stack) ,
                self.file_type_ = 'in-memory' ;
                self.typed_file_ = lapwing.Stack_file_in_memory(file_name_or_stack) ;                
            else
                error('Stack_file:UnableToLoad', 'Unable to load a thing like that') ;
            end
            self.i_frame_ = 0 ;
        end  % constructor method
        
        function delete(self)
            self.typed_file_ = [] ;
        end
        
        function result = data_type(self)
            result = self.typed_file_.data_type() ;
        end
        
        function result = n_row(self)
            result = self.typed_file_.n_row() ;
        end
        
        function result = n_col(self)
            result = self.typed_file_.n_col() ;
        end
        
        function result = n_frame(self)
            result = self.typed_file_.n_frame() ;
        end

        function result = rate(self)
            result = self.typed_file_.rate() ;
        end        
        
        function result = dt(self)  % s, 1/rate
            rate = self.rate() ;
            if isempty(rate) ,
                result=[];
            else
                result=1/rate;
            end
        end
        
        function frame=get_frame(self, i)
            frame = self.typed_file_.get_frame(i) ;
            self.i_frame_ = i ;
        end
        
        function frame=get_next(self)
            i_frame = self.i_frame_ ;
            frame = self.typed_file_.get_next(i_frame) ;
            self.i_frame_ = i_frame + 1 ;
        end
        
        function frame=get_prev(self)
            frame=self.get_frame(self.i_frame_);
        end
        
        function pred=at_end(self)
            pred=(self.i_frame_==self.n_frame());
        end
        
        function pred=at_start(self)
            pred=(self.i_frame_==0);
        end
        
        function to_start(self)
            self.typed_file_.to_start() ;
            self.i_frame_=0;
        end
    end  % methods
    
end  % classdef

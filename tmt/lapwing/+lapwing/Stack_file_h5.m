classdef Stack_file_h5 < handle
    
    properties
        file_path_
        dataset_path_
        data_type_
        n_row_
        n_col_
        n_frame_     
        has_leading_singleton_dimension_
    end  % properties
        
    methods
        function self = Stack_file_h5(file_name)
            file_path = lapwing.absolute_filename(file_name) ;
            info = h5info(file_path, '/') ;
            dataset_count = length(info.Datasets) ;
            if dataset_count ~= 1 ,
                error('Stack_file:UnableToLoad', ...
                      'Stack_file_h5 only supports files with exactly one dataset in the root group');
            end                
            dataset_name = info.Datasets.Name ;
            dataset_path = ['/' dataset_name] ;
            dataset_info = h5info(file_path, dataset_path) ;
            h5_data_type = dataset_info.Datatype.Type ;
            if isequal(h5_data_type, 'H5T_STD_U16LE') ,
                self.data_type_ = 'uint16' ;
            elseif isequal(h5_data_type, 'H5T_IEEE_F32LE') ,
                self.data_type_ = 'single' ;
            else
                error('Stack_file:UnsupportedPixelType', ...
                      'Stack_file_h5 only supports 16-bit grayscale stacks.');                
            end
            raw_shape = dataset_info.Dataspace.Size ;
            if length(raw_shape)>3 ,
                if raw_shape(1) > 1 , 
                    error('Stack_file:UnableToLoad', ...
                          'Stack_file_h5 only supports 3D data sets');
                end
                shape_ijk = raw_shape(2:end) ;  % assumed xyz order                
                self.has_leading_singleton_dimension_ = true ;
            else
                shape_ijk = raw_shape ;
                self.has_leading_singleton_dimension_ = false ;
            end
            if length(shape_ijk) ~= 3 ,
                error('Stack_file:UnableToLoad', ...
                      'Stack_file_h5 only supports 3D data sets');
            end
            
            % Set the state variables
            self.file_path_ = file_path ;
            self.dataset_path_ = dataset_path ;
            self.n_row_ = shape_ijk(2) ;
            self.n_col_ = shape_ijk(1) ;
            self.n_frame_ = shape_ijk(3) ;
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
            if self.has_leading_singleton_dimension_ ,
                raw_frame = h5read(self.file_path_, self.dataset_path_, [1 1 1 i_frame], [1 self.n_col_ self.n_row_ 1]) ;
                frame = reshape(raw_frame, [self.n_col_ self.n_row_])' ;
            else
                raw_frame = h5read(self.file_path_, self.dataset_path_, [1 1 i_frame], [self.n_col_ self.n_row_ 1]) ;
                frame = raw_frame' ;  % want n_y x n_x
            end
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

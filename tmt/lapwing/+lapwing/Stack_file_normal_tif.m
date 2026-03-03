classdef Stack_file_normal_tif < handle
    
    properties
        tiff_object_  % the file-type specific file object
        original_libtiff_warning_state_  % used to restore libtiff warning state after we close a tiff file        
        original_libtiff_error_as_warning_state_
        n_frame_
        data_type_
        n_row_
        n_col_
    end  % properties
        
    methods
        function self = Stack_file_normal_tif(file_name)
            info=imfinfo(file_name);
            self.n_frame_=length(info);
            self.original_libtiff_warning_state_ = warning('query', 'MATLAB:imagesci:tiffmexutils:libtiffWarning') ;
            warning('off','MATLAB:imagesci:tiffmexutils:libtiffWarning');
            self.original_libtiff_error_as_warning_state_ = warning('query', 'MATLAB:imagesci:tiffmexutils:libtiffErrorAsWarning') ;
            warning('off','MATLAB:imagesci:tiffmexutils:libtiffErrorAsWarning');            
            self.tiff_object_=Tiff(file_name,'r');
            frame=self.tiff_object_.read();
            if ndims(frame)>2  %#ok
                error('Stack_file:UnsupportedPixelType', ...
                      'Stack_file_normal_tif only supports 8- and 16-bit grayscale videos.');
            end
            if isa(frame,'uint8')
                self.data_type_='uint8';
            elseif isa(frame,'uint16')
                self.data_type_='uint16';
            else
                error('Stack_file:UnsupportedPixelType', ...
                      'Stack_file_normal_tif only supports 8- and 16-bit grayscale videos.');
            end
            [self.n_row_,self.n_col_]=size(frame);
            self.tiff_object_.close();
            self.tiff_object_ = [] ;
            self.tiff_object_=Tiff(file_name,'r');
        end  % constructor method

        function delete(self)
            self.tiff_object_.close();
            warning(self.original_libtiff_warning_state_);
            warning(self.original_libtiff_error_as_warning_state_);
            self.tiff_object_ = [] ;
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

        function frame=get_frame(self,i_frame)
            self.tiff_object_.setDirectory(i_frame);
            frame=self.tiff_object_.read();
        end
        
        function frame=get_next(self, i_last_frame_read) %#ok<INUSD>
            % Gets frame *after* i_last_frame_read
            frame=self.tiff_object_.read();  % if i_last_frame_read is right, this should work fine
        end        
        
        function to_start(self)
            if self.n_frame_>0 ,
                self.tiff_object_.setDirectory(1);
            end
        end
        
    end  % methods
    
end  % classdef

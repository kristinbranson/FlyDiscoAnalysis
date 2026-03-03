classdef Stack_file_imagej_jumbo_tif < handle
    
    properties
        header_length_  % in bytes
        frame_stride_  % in bytes, frame i (one-based) starts at file offset header_length_ + (i-1) * frame_stride_
        frame_size_
        n_frame_
        %data_type_
        n_row_
        n_col_
        bytes_per_pel_
        pels_per_frame_
        fid_
        next_frame_offset_
        pixel_class_name_
        byte_order_
    end  % properties
        
    methods
        function self = Stack_file_imagej_jumbo_tif(file_name)
            [tiff_struct, tiff_header] = lapwing.tiffread31_header(file_name) ;
            self.fid_ = tiff_struct.file ;
            tiff_tags = lapwing.tiffread31_readtags(tiff_struct, tiff_header, 1) ;
            [is_imagej_jumbo_tiff, is_imagej_hyperstack, stack_dims] = ...
                lapwing.IsImageJBigStack(tiff_tags, numel(tiff_header));
            if ~is_imagej_jumbo_tiff ,
                error('File does not appear to be an ImageJ "jumbo" TIFF file') ;                
            end
            if is_imagej_hyperstack ,
                error('Unable to open a multi-channel ImageJ "jumbo" TIFF file') ;
            end
            self.n_frame_ = stack_dims(3) ;
            data_type = tiff_header.bits ;
            self.n_row_ = tiff_header.width ;
            self.n_col_ = tiff_header.height ;            
            self.header_length_ = tiff_header.StripOffsets ;
            self.pels_per_frame_ = self.n_row * self.n_col ;
            self.bytes_per_pel_ = ceil(data_type/8) ;    
            if self.bytes_per_pel_==2 ,
                self.pixel_class_name_ = 'uint16' ;
            else
                self.pixel_class_name_ = 'uint8' ;
            end                            
            bytes_per_frame = self.pels_per_frame_ * self.bytes_per_pel_ ;
            self.frame_stride_ = bytes_per_frame ;
            self.frame_size_ = bytes_per_frame ;
            self.byte_order_ = tiff_struct.ByteOrder ;
            self.next_frame_offset_ = self.header_length_ ;
            fseek(self.fid_, self.next_frame_offset_, 'bof') ;            
        end  % constructor method
        
        function delete(self)
            if ~isempty(self.fid_) ,
                fclose(self.fid_) ;
            end
            self.fid_ = [] ;
        end
        
        function result = data_type(self)
            result = self.pixel_class_name_ ;
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
            % - Seek file to beginning of chunk
            offset = self.header_length_ + (i_frame-1) * self.frame_stride_ ;
            fseek(self.fid_, offset, 'bof') ;
            serial_frame = fread(self.fid_, self.pels_per_frame_, [self.pixel_class_name_ '=>' self.pixel_class_name_], 0, self.byte_order_) ;
            frame = reshape(serial_frame, [self.n_row_ self.n_col_])' ;
            self.next_frame_offset_ = offset + self.frame_stride_ ;
        end
        
        function frame=get_next(self, i_last_frame_read) %#ok<INUSD>
            % Could potentially make faster by keeping a file offset around
            serial_frame = fread(self.fid_, self.pels_per_frame_, [self.pixel_class_name_ '=>' self.pixel_class_name_], 0, self.byte_order_) ;
            frame = reshape(serial_frame, [self.n_row_ self.n_col_])' ;
            self.next_frame_offset_ = self.next_frame_offset_ + self.frame_stride_ ;
        end        
        
        function to_start(self)
            self.next_frame_offset_ = self.header_length_ ;
            fseek(self.fid_, self.next_frame_offset_, 'bof') ;            
        end
    end  % methods
    
end  % classdef

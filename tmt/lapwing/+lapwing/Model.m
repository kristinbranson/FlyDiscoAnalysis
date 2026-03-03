classdef Model < handle
    
    properties (Access = private)
        cmap_name_
        cmap_
        fps_
        mode_
        file_name_  % the name of the file currently open
        file_  % the handle of a VideoFile object, the current file (or empty)
        % colorbar_min and colorbar_max are constrained to be integers
        colorbar_max_string_
        colorbar_min_string_
        colorbar_min_  % the colorbar min, derived from cb_min_string,
        colorbar_max_  % the colorbar max, derived from cb_min_string,        
        was_stop_button_hit_
        z_index_
    end
    
    properties (Dependent=true)
        file_name  % the name of the file currently open
        file  % the handle of a VideoFile object, the current file (or empty)
        % colorbar_min and colorbar_max are constrained to be integers
        colorbar_max_string
        colorbar_min_string
        colorbar_min  % the colorbar min, derived from cb_min_string,
        colorbar_max  % the colorbar max, derived from cb_min_string,        
        fps  % Hz, sampling rate for playback & export, possibly different from
        % that in file
        n_rows
        n_cols
        z_slice_count  % number of time samples
        %tl  % 2x1 matrix holding min, max time
        %n_rois;
        %t  % a complete timeline for all z_slices
        is_a_file_open  % true iff a stack is currently open
        z_slice
        indexed_z_slice
        cmap_name
        cmap
        mode
        was_stop_button_hit
        z_index
    end
    
    methods
        function self = Model()
            self.file_name_ = '' ;
            self.file_ = [] ;
            self.z_index_ = [] ;
            self.fps_ = 20 ;  % for playback
            self.cmap_name_ = 'gray' ; 
            
            colorbar_min_string='0';
            colorbar_max_string='255';
            colorbar_min=str2double(colorbar_min_string);
            colorbar_max=str2double(colorbar_max_string);
            
            self.colorbar_max_string_ = colorbar_max_string ;
            self.colorbar_min_string_ = colorbar_min_string ;
            self.colorbar_max_ = colorbar_max ;
            self.colorbar_min_ = colorbar_min ;            
            
            self.mode_ = 'zoom' ;
            self.was_stop_button_hit_ = false ;
        end  % function
        
        function delete(self)
            if ~isempty(self.file_) ,
                self.file_ = [] ;  % should close file handles
            end
        end            
        
        function result = get.fps(self)
            result = self.fps_ ;
        end
        
        
        function result = get.n_rows(self)
            if isempty(self.file_)
                result = [] ;
            else
                result = self.file_.n_row ;
            end
        end
        
        
        function result = get.n_cols(self)
            if isempty(self.file_)
                result = [] ;
            else
                result = self.file_.n_col ;
            end
        end
        
        
        function result = get.z_slice_count(self)
            if isempty(self.file_) ,
                result = [] ;
            else
                result = self.file_.n_frame ;
            end
        end
        
        
        function set.fps(self, new_value)
            if isnumeric(new_value) && isreal(new_value) && isscalar(new_value) && new_value>0 ,
                self.fps_ = double(new_value) ;
            end
        end       
        
        
        function open_file_given_file_name_or_stack(self, file_name_or_stack)
            % filename is a filename, can be relative or absolute
            
            % break up the file name
            %[~,base_name,ext]=fileparts(filename);
            %filename_local=[base_name ext];
            
            % load the optical data
            file = lapwing.Stack_file(file_name_or_stack) ;
                        
            % set the model
            self.file_ = file ;
            if ischar(file_name_or_stack) ,
                self.file_name_ = file_name_or_stack ;
            else
                self.file_name_ = '(in-memory stack)' ;
            end
            
            % determine the colorbar bounds
            [data_min,data_max] = lapwing.pixel_data_type_min_max(self.file_.data_type) ;
            self.colorbar_min_string_ = sprintf('%d',data_min) ;
            self.colorbar_max_string_ = sprintf('%d',data_max) ;
            self.colorbar_min_ = str2double(self.colorbar_min_string) ;
            self.colorbar_max_ = str2double(self.colorbar_max_string) ;

            % reset the z_slice index
            self.z_index_ = 1 ;

            % set the mode to zoom
            self.mode_ = 'zoom' ;                        
        end  % method
        
        
        function close_file(self)
            if ~isempty(self.file_) ,
                self.file_ = [] ;  % should close file handles
            end
            self.file_name_ = '' ;
        end  % method
        
        
        function result = get.is_a_file_open(self)
            result=~isempty(self.file_);
        end  % method

        
        function result = get.z_slice(self)
            % Get the current indexed_z_slice, based on model, z_index,
            % colorbar_min, and colorbar_max.
            z_index = self.z_index_ ;
            result = self.file_.get_frame(z_index) ;
        end
        
        
        function indexed_z_slice = get.indexed_z_slice(self)
            % Get the current indexed_z_slice, based on model, z_index,
            % colorbar_min, and colorbar_max.
            z_slice = double(self.z_slice)  ;
            cb_min = self.colorbar_min ;
            cb_max = self.colorbar_max ;
            indexed_z_slice = uint8(round(255*(z_slice-cb_min)/(cb_max-cb_min))) ;
        end
        
        
        function set_colorbar_bounds_from_strings(self, new_cb_min_string, new_cb_max_string)            
            % Set the view colorbar bounds given max and min values in strings.  No
            % checking is done to make sure the string values are sane.  Note that in
            % the view, although both string and numerical representations of the
            % colorbar bounds are maintained, the string ones are the more fundamental,
            % and the numerical ones are derived from them.  At the moment, this
            % doesn't matter much, since the bounds are constrained to always be
            % integral.  But if we support movies with floating-point pels at some
            % point, this will be important, since the user can explicitly set the
            % bounds, and we want to keep hold of _exactly_ what the user typed.

            % convert to numbers
            new_cb_min = str2double(new_cb_min_string) ;
            new_cb_max = str2double(new_cb_max_string) ;
            
            % if new values are kosher, change colorbar bounds
            if ~isempty(new_cb_min) && ~isempty(new_cb_max) && ...
                    isfinite(new_cb_min) && isfinite(new_cb_max) && ...
                    (new_cb_max>new_cb_min) ,           
                % store the strings
                self.colorbar_min_string_ = new_cb_min_string ;
                self.colorbar_max_string_ = new_cb_max_string ;

                % store the numbers
                self.colorbar_min_ = new_cb_min ;
                self.colorbar_max_ = new_cb_max ;           
            end
        end

        
        function set_colorbar_bounds_from_numbers(self, new_cb_min, new_cb_max)            
            % Set the view colorbar bounds given max and min values as numbers.  This
            % calls self.set_colorbar_bounds_from_strings().  It is assumed that cb_min
            % and cb_max are doubles that happen to be integral.
            
            % if new values are kosher, change colorbar bounds
            if ~isempty(new_cb_min) && ~isempty(new_cb_max) && ...
                    isfinite(new_cb_min) && isfinite(new_cb_max) && ...
                    (new_cb_max>new_cb_min) ,           
                % convert to strings
                new_cb_min_string = sprintf('%d',new_cb_min) ;
                new_cb_max_string = sprintf('%d',new_cb_max) ;            

                % store the strings
                self.colorbar_min_string_ = new_cb_min_string ;
                self.colorbar_max_string_ = new_cb_max_string ;

                % store the numbers
                self.colorbar_min_ = new_cb_min ;
                self.colorbar_max_ = new_cb_max ;           
            end
        end
        
        
        function result = get.cmap_name(self) 
            result = self.cmap_name_ ;
        end
        
        
        function set.cmap_name(self, new_cmap_name)
            % set the chosen cmap_name
            self.cmap_name_ = new_cmap_name ;

            % set the colormap
            self.cmap_ = lapwing.cmap_from_name(new_cmap_name) ;
        end        
        
        
        function result = get.cmap(self) 
            result = self.cmap_ ;
        end
        
        
        function brighten(self)            
            cmap = self.cmap_ ;
            self.cmap_ = brighten(cmap, 0.1) ;
        end
        
        
        function revert_gamma(self)
            self.cmap_ = lapwing.cmap_from_name(self.cmap_name) ;
        end
        
        
        function darken(self)
            cmap = self.cmap_ ;
            self.cmap_ = brighten(cmap, -0.1) ;
        end

        
        function set.z_index(self, new_value)
            % Change the current z_slice to the given z_slice index
            z_slice_count = self.z_slice_count ;
            if (new_value>=1) && (new_value<=z_slice_count) ,
                self.z_index_ = new_value ;
            end
        end
        
        
        function set_colorbar_bounds(self, method, bounds)            
            switch(method)
                case 'pixel_data_type_min_max'
                    [d_min, d_max] = lapwing.pixel_data_type_min_max(self.file_.data_type) ;
                    self.set_colorbar_bounds_from_numbers(d_min,d_max);
                case 'min_max'
                    z_slice = double(self.z_slice) ;
                    d_min = min(min(z_slice)) ;
                    d_max = max(max(z_slice)) ;
                    self.set_colorbar_bounds_from_numbers(d_min,d_max);
                case 'five_95'
                    z_slice=double(self.z_slice);
                    d=lapwing.quantile_mine(z_slice(:),[0.05 0.95]');
                    d_05=d(1);
                    d_95=d(2);
                    d_05=floor(d_05);  % want int, want to span range
                    d_95=ceil(d_95);  % want int, want to span range
                    self.set_colorbar_bounds_from_numbers(d_05,d_95);
                case 'abs_max'
                    z_slice=double(self.z_slice);
                    cb_max=max(max(abs(z_slice)));
                    self.set_colorbar_bounds_from_numbers(-cb_max,+cb_max);
                case 'ninety_symmetric'
                    % need to fix this, since what it does now is useless
                    z_slice=abs(double(self.z_slice));
                    d_90=lapwing.quantile_mine(z_slice(:),0.9);                    
                    cb_max=ceil(d_90);
                    cb_min=-cb_max;
                    self.set_colorbar_bounds_from_numbers(cb_min,cb_max);
                case 'manual'
                    % break out the returned cell array
                    new_cb_max_string=bounds{1};
                    new_cb_min_string=bounds{2};
                    % convert all these strings to real numbers
                    new_cb_min = floor(str2double(new_cb_min_string)) ; 
                    new_cb_max = ceil(str2double(new_cb_max_string)) ;
                    % change colorbar bounds
                    self.set_colorbar_bounds_from_numbers(new_cb_min, new_cb_max) ;
            end
        end  % method
        
        
        function set.mode(self, new_value)
            if isequal(new_value, 'zoom') ,
                self.mode_ = 'zoom' ;
            end
        end
        
        
        function result = get.mode(self)
            result = self.mode_ ;
        end
        
        
        function result = get.file_name(self)
            result = self.file_name_ ;
        end
        
        
        function result = get.colorbar_max_string(self)
            result = self.colorbar_max_string_ ;
        end
        
        
        function result = get.colorbar_min_string(self)
            result = self.colorbar_min_string_ ;
        end

        
        function result = get.colorbar_max(self)
            result = self.colorbar_max_ ;
        end
        
        
        function result = get.colorbar_min(self)
            result = self.colorbar_min_ ;
        end
        
        
        function result = get.was_stop_button_hit(self)
            result = self.was_stop_button_hit_ ;
        end
        

        function set.was_stop_button_hit(self, new_value)
            if (islogical(new_value) || isnumeric(new_value)) && isscalar(new_value) ,
                self.was_stop_button_hit_ = logical(new_value) ;
            end
        end
        
        
        function result = get.z_index(self)
            result = self.z_index_ ;
        end
        
    end  % methods
    
end  % classdef

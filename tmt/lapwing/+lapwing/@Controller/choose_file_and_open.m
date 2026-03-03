function choose_file_and_open(self)
    
    % throw up the dialog box
    [filename,pathname]= ...
        uigetfile({'*.tif','Multi-image TIFF files (*.tif)'; ...
        '*.h5','HDF5 files (*.h5)'; ...
        '*.mj2','Motion JPEG 2000 files (*.mj2)'}, ...
        'Load stack from file...');
    if isnumeric(filename) || isnumeric(pathname)
        % this happens if user hits Cancel
        return;
    end
    
    % Open the file
    file_name_full = fullfile(pathname, filename) ;
    self.open_file_given_file_name_or_stack(file_name_full) ;    
end

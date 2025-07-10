function open_file_given_file_name_or_stack(self, file_name_or_stack)
    
    % filename is a filename, can be relative or absolute
    
    % break up the file name
    if ischar(file_name_or_stack) ,
        [~, base_name, ext] = fileparts(file_name_or_stack) ;
        filename_local = [base_name ext] ;
    else
        filename_local = 'in-memory stack' ;
    end
    
    % load the optical data
    self.hourglass()
    try
        self.model.open_file_given_file_name_or_stack(file_name_or_stack) ;
    catch err
        self.unhourglass();
        if strcmp(err.identifier,'MATLAB:imagesci:imfinfo:whatFormat')
            errordlg(sprintf('Unable to determine file format of %s', ...
                filename_local),...
                'File Error');
            return;
        elseif strcmp(err.identifier,'MATLAB:load:notBinaryFile')
            errordlg(sprintf('%s does not seem to be a binary file.', ...
                filename_local),...
                'File Error');
            return;
        elseif isempty(err.identifier)
            errordlg(sprintf('Unable to read %s, and error lacks identifier.', ...
                filename_local),...
                'File Error');
            return;
        elseif strcmp(err.identifier,'Stack_file:UnableToLoad')
            errordlg(sprintf('Error opening %s.',filename_local),...
                'File Error');
            return;
        elseif strcmp(err.identifier,'Stack_file:UnsupportedPixelType')
            errordlg(sprintf('Error opening %s: Unsupported pixel type.', ...
                filename_local),...
                'File Error');
            return;
        else
            rethrow(err);
        end
    end
    
    % Update the view HG objects to match the model & view state
    self.update();
    
    % OK, we're done.
    self.unhourglass()    
end

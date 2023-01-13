function APT_Labeler_printInfo_without_adding_APT_to_path(label_file_name)

% Get the path, set up a cleaner to set it back to that when we exit this function
original_path = path() ;
cleaner = onCleanup(@()(path(original_path))) ;

% Want to bracket this stuff in the log
fprintf('*** Start of APT Labeler.printInfo() output ***\n') ;

% Where does this script live?
this_script_path = mfilename('fullpath') ;
fly_disco_analysis_folder_path = fileparts(this_script_path) ;
apt_folder_path = fullfile(fly_disco_analysis_folder_path, 'APT') ;

% Set up the path for APT stuff
addpath(apt_folder_path) ;
APT.setpathsmart() ;

% Create a Labeler object using the label file, just to can printInfo() method
labeler = Labeler('isgui', false, 'projfile', label_file_name) ;
labeler.printInfo() ;

% Not sure if I need to do this, but can't hurt
delete(labeler)

% Want to bracket this stuff in the log
fprintf('*** End of APT Labeler.printInfo() output ***\n') ;

end

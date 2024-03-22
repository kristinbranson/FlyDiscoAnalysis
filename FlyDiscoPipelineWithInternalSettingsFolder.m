function FlyDiscoPipelineWithInternalSettingsFolder(expdir, varargin)
%FlyDiscoPipelineWithInternalSettingsFolder  Runs the FlyDisco pipeline, using the internal settings folder.
%   FlyDiscoPipelineWithInternalSettingsFolder(expdir) runs the FlyDisco
%   pipeline on the experiment folder expdir, with the 'settingsdir' argument
%   set to the path to the internal settings folder.  It should be possible to
%   call this function with only the expdir argument and get the pipeline to
%   run on the expdir, but you can also pass key-value pairs that will be
%   passed through to FlyDiscoPipeline().  See help for FlyDiscoPipeline() for
%   more details.

this_script_path = mfilename('fullpath') ;
source_folder_path = fileparts(this_script_path) ;
internal_settings_folder_path = fullfile(source_folder_path, 'settings-internal') ;
FlyDiscoPipeline(expdir, varargin{:}, 'settingsdir', internal_settings_folder_path) ;

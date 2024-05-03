function FlyDiscoComputeAptPerFrameFeatures(expdir, varargin)
% This function implements the computeaptperframefeatures stage.

% varargin will be a sequence of key-value pairs, including at least the keys
% 'settingsdir', 'analysis_protocol', and 'forcecompute'.  It will also
% contain any additional stage-specific key-value pairs passed as the last
% argument to FlyDiscoPipelineStage().

% Parse the optional arguments
[settingsdir, analysis_protocol,forcecompute] = ...
  myparse(varargin,...
          'settingsdir', default_settings_folder_path(), ...
          'analysis_protocol', 'current', ...
          'forcecompute', false) ;

% Do more stuff here...

% If something goes wrong, throw an error.  Any values returned from this
% function are ignored.

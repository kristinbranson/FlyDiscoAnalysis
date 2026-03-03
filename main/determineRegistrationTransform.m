function registration_data = determineRegistrationTransform(expdir, analysis_protocol_folder_path, dataloc_params, registration_params)
% Determine the registration transform.  Along the way, detect the
% registration mark(s), even though they are currently not used for
% registration.

method = registration_params.method ;
if strcmpi(method, 'circle') ,
  registration_data = determineRegistrationTransformForCircleMethod(expdir, analysis_protocol_folder_path, dataloc_params, registration_params) ;
elseif strcmpi(method, 'flytracker') ,
  registration_data = determineRegistrationTransformForFlyTrackerMethod(expdir, analysis_protocol_folder_path, dataloc_params, registration_params) ;
else
  error('Unsupported registration method: %s', method) ;
end

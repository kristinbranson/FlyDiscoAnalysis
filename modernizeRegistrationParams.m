function registration_params = modernizeRegistrationParams(raw_registration_params)

registration_params = raw_registration_params ;
if ~isfield(registration_params, 'plotYAxisPointsUp') ,
  registration_params.plotYAxisPointsUp = 1 ;
    % This makes registration/ctrax-result/apt-results frames appear upside-down
    % (actually flipped) relative to how they're shown by e.g. playfmf(), but it
    % was the norm historically.
end

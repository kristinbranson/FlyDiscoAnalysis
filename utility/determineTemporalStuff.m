function [dotemporalregout, dotemporaltruncationout, recordLengthIdeal] = determineTemporalStuff(registration_params, dotemporalreg, dotemporaltruncation)
% Determine the final settings for temporal registration and truncation.
% This is a pure function.

if isfield(registration_params,'doTemporalRegistration'),
  dotemporalregout = registration_params.doTemporalRegistration;
else
  dotemporalregout = dotemporalreg ;
end
if isfield(registration_params,'doTemporalTruncation')
  if registration_params.doTemporalTruncation > 0
    dotemporaltruncationout = true;
    recordLengthIdeal = registration_params.doTemporalTruncation ;
  else
    dotemporaltruncationout = false;
    recordLengthIdeal = [] ;
  end
else
  dotemporaltruncationout = dotemporaltruncation ;
  recordLengthIdeal = [] ;
end

end  % function

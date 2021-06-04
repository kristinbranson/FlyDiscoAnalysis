% read per-frame params
function ReadPerFrameParams(obj)

obj.perframe_params = ReadParams(obj.dataloc_params.perframeparamsfile);
if ~isfield(obj.perframe_params,'isflytracker'),
  warning('isflytracker not set in perframe parameters, using default value false');
  obj.perframe_params.isflytracker = false;
end
if ~isfield(obj.perframe_params,'fakectrax'),
  warning('fakectrax not set in perframe parameters, using default value true');
  obj.perframe_params.fakectrax = true;
end
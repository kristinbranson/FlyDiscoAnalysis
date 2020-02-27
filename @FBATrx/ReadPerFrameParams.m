% read per-frame params
function ReadPerFrameParams(obj)

obj.perframe_params = ReadParams(obj.dataloc_params.perframeparamsfile);


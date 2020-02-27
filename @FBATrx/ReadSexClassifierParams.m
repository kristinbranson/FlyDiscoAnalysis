% read sex classifier params
function ReadSexClassifierParams(obj)

obj.sexclassifier_params = ReadParams(obj.dataloc_params.sexclassifiertxtfile);


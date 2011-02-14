function ReadAllParams(obj,varargin)

% all extra arguments are parameters
for i = 1:2:nargin-1,
  obj.(varargin{i}) = varargin{i+1};
end

%% read data location parameters
datalocparamsfile = fullfile(obj.settingsdir,obj.analysis_protocol,obj.datalocparamsfilestr);
if ~exist(datalocparamsfile,'file'),
  error('Data location parameters file %s does not exist',datalocparamsfile);
end
dataloc_params = ReadParams(datalocparamsfile);

% save
obj.dataloc_params = dataloc_params;
obj.dataloc_params.datalocparamsfile = datalocparamsfile;

%% read landmark params

obj.dataloc_params.landmarkparamsfile = ...
  fullfile(obj.settingsdir,obj.analysis_protocol,obj.dataloc_params.landmarkparamsfilestr);
obj.ReadLandmarkParams();

%% read per-frame params

obj.dataloc_params.perframeparamsfile = ...
  fullfile(obj.settingsdir,obj.analysis_protocol,obj.dataloc_params.perframeparamsfilestr);
obj.ReadPerFrameParams();

%% read sex classifier

obj.dataloc_params.sexclassifiertxtfile = ...
  fullfile(obj.settingsdir,obj.analysis_protocol,obj.dataloc_params.sexclassifiertxtfilestr);
obj.ReadSexClassifierParams();

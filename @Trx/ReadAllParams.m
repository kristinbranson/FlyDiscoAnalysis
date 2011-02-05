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
old_dataloc_params = obj.dataloc_params;
obj.dataloc_params = dataloc_params;
obj.dataloc_params.datalocparamsfile = datalocparamsfile;

%% reset global locations of data, and reload

if ~isempty(old_dataloc_params) && ...
    (~strcmp(old_dataloc_params.rootreaddir,dataloc_params.rootreaddir) || ...
    ~strcmp(old_dataloc_params.rootwritedir,dataloc_params.rootwritedir)),
  expdir_bases = obj.expdir_bases;
  for i = 1:numel(expdir_bases),
    obj.RemoveExpDir(expdir_bases{i});
  end
  for i = 1:numel(expdir_bases),
    obj.AddExpDir(expdir_bases{i});
  end
end

%% read landmark params

obj.dataloc_params.landmarkparamsfile = ...
  fullfile(obj.settingsdir,obj.analysis_protocol,obj.dataloc_params.landmarkparamsfilestr);
obj.ReadLandmarkParams();

%% read per-frame params

obj.dataloc_params.perframeparamsfile = ...
  fullfile(obj.settingsdir,obj.analysis_protocol,obj.dataloc_params.perframeparamsfilestr);
obj.ReadPerFrameParams();

%% read sex classification params

obj.dataloc_params.sexclassifierparamsfile = ...
  fullfile(obj.settingsdir,obj.analysis_protocol,obj.dataloc_params.sexclassifierparamsfilestr);
obj.ReadSexClassifierParams();

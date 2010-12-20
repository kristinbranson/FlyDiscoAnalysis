% metadata = obj.getMetaDataField(fns,['expdir',expdirs],['n',ns])
% Get metadata field(s) for the specified experiment(s). If multiple
% fields are specified or multiple experiments are specified, the
% result will be a cell of size nexpdirs x nfields.
% Either the expdirs or ns must be specified. If both are specified,
% then ns is used.
function metadata = getMetaDataField(obj,fns,varargin)

if nargin < 3,
  error('Either expdir or n must be input');
end

[expdirs,ns] = myparse(varargin,'expdir',{},'n',[]);

singlefn = ischar(fns);
if singlefn,
  fns = {fns};
end

singleexp = false;
if isempty(ns),
  if ischar(expdirs),
    singleexp = true;
    expdirs = {expdirs};
  end
  ns = obj.exp2n(expdirs);
else
  singleexp = numel(ns) == 1;
end

metadata = cell(numel(ns),numel(fns));
for i = 1:numel(ns),
  n = ns(i);
  for j = 1:numel(fns),
    fn = fns{j};
    switch fn,
      case {'assay','exp_datetime','experimenter','protocol'},
        metadata{i,j} = obj.metadata{n}.attribute.(fn);
      case {'humidity','temperature'},
        metadata{i,j} = obj.metadata{n}.getChildrenByName('environment').attribute.(fn);
      case {'flag_aborted','flag_redo','flag_review','note_behavioral','note_technical'},
        metadata{i,j} = obj.metadata{n}.getChildrenByName('flag_aborted').content;
      case {'bowl','camera','computer','harddrive','plate','rig'},
        metadata{i,j} = obj.metadata{n}.getNodeByUniquePath({'experiment','session','apparatus'}).attribute.(fn);
      case {'count','cross_date','effector','gender','hours_starved','line'},
        metadata{i,j} = obj.metadata{n}.getNodeByUniquePath({'experiment','session','flies'}).attribute.(fn);
      case 'genotype',
        metadata{i,j} = obj.metadata{n}.getNodeByUniquePath({'experiment','session','flies','genotype'}).content;
      case {'datetime_sorting','datetime_starvation','handler_sorting','handler_starvation',...
          'hours_sorted','seconds_fliesloaded','seconds_shiftflytemp'},
        metadata{i,j} = obj.metadata{n}.getNodeByUniquePath({'experiment','session','flies','handling'}).attribute.(fn);
      case 'handling_protocol',
        metadata{i,j} = obj.metadata{n}.getNodeByUniquePath({'experiment','session','flies','handling'}).attribute.protocol;
      case 'incubator',
        metadata{i,j} = obj.metadata{n}.getNodeByUniquePath({'experiment','session','flies','rearing'}).attribute.(fn);
      case 'rearing_protocol',
        metadata{i,j} = obj.metadata{n}.getNodeByUniquePath({'experiment','session','flies','rearing'}).attribute.protocol;
    end
    
    % numeric data
    if ischar(metadata{i,j}) && ismember(fn,{'humidity','temperature','count','hours_starved',...
        'hours_sorted','seconds_fliesloaded','seconds_shiftflytemp'}),
      metadata{i,j} = str2double(metadata{i,j});
    end
    
  end
  
end

if singleexp && singlefn,
  metadata = metadata{1,1};
end
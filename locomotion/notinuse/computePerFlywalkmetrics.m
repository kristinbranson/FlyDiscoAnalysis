function [perfly_walkmetrics] = computePerFlywalkmetrics(perwalk_metrics)


for fly = 1:numel(perwalk_metrics)
    walkfeaturestruct = perwalk_metrics

% % combine across walks
perflywalkfeatures = struct;

% create perfly data
% combine fields within each fly
flds = fields(walkfeaturestruct);

for fld = 1:numel(flds)
    % for non-structure fields just concatenate
    if ~isstruct(walkfeaturestruct(1).(flds{fld}))

        perflywalkfeatures.(flds{fld})= [walkfeaturestruct.(flds{fld})];


    % combine perframe data
    elseif any(strcmp(pfflist_first,flds{fld} )) || any(strcmp(pfflist_none,flds{fld}))
        % isstruct(walkfeaturestruct(1).(flds{fld})) & 
        fieldname = flds{fld};
        perflyperframefeatures = computePerFlywalkperframefeatures(walkfeaturestruct,fieldname);
        perflywalkfeatures.(flds{fld}) = perflyperframefeatures.(flds{fld});

    elseif startsWith(flds{fld}, 'phase')
        fieldname = flds{fld};

        perflyphasefeatures = computePerFlyphasefeatures(walkfeaturestruct,fieldname);
        perflywalkfeatures.(flds{fld}) = perflyphasefeatures.(flds{fld});

    end
end

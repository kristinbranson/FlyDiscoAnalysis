function walkmetrics = computeWalkMetrics(obj,digital_signal)


% params


perwalk_metrics = struct;
perfly_metrics = struct;

nflies = numel(obj.limbBoutData);

% compute walk metrics for each fly 
for fly = 1:nflies
    currfly_digitalindicator = digital_signal{fly};
    [walkfeaturestruct,perflywalkfeatures] = compute_WalkFeatures(obj,fly,currfly_digitalindicator);
    if fly == 1
        perwalk_metrics = walkfeaturestruct;
        perfly_metrics = perflywalkfeatures;
    else 
     perwalk_metrics = [perwalk_metrics, walkfeaturestruct];
     perfly_metrics = [perfly_metrics,perflywalkfeatures];
    end
end

% combine across flies
perexp_metrics = struct;

% list of perframe features to compute stats over bouts that are 
%  'first' = 1st derivative (d)  
pfflist_first = {'velmag_ctr','absdv_ctr','absdu_ctr','absdtheta', ...
    'left_vel','right_vel','forward_vel','backward_vel','right_dtheta','left_dtheta'};   
% 'none' = no derivative
pfflist_none = {'CoM_stability'};
% 'second' = second derivative (dd)

flds = fields(perfly_metrics);

for fld = 1:numel(flds)
    % for non-structure fields just concatenate
    if ~isstruct(perfly_metrics(1).(flds{fld}))

        perexp_metrics.(flds{fld})= [perfly_metrics.(flds{fld})];

    % combine perframe data
    elseif any(strcmp(pfflist_first,flds{fld} )) || any(strcmp(pfflist_none,flds{fld}))
        % isstruct(walkfeaturestruct(1).(flds{fld})) & 
        fieldname = flds{fld};
        perexpperframefeatures = computePerExpwalkperframefeatures(perwalk_metrics,perfly_metrics,fieldname);
        perexp_metrics.(flds{fld}) = perexpperframefeatures.(flds{fld});

    elseif startsWith(flds{fld}, 'phase')
        fieldname = flds{fld};

        perexpphasefeatures = computePerExpphasefeatures(perwalk_metrics,perfly_metrics,fieldname);
        perexp_metrics.(flds{fld}) = perexpphasefeatures.(flds{fld});

    end
end


walkmetrics = struct;
walkmetrics.perwalk = perwalk_metrics;
walkmetrics.perfly = perfly_metrics;
walkmetrics.perexp = perexp_metrics;




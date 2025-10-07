function [walkfeaturestruct,perflywalkfeatures] = compute_WalkFeatures(obj,fly,digital_signal)

% compute_WalkFeatures - Extract walking features from fly limb tracking data
%
% Computes various kinematic and phase-based features for individual walking 
% bouts, including perframe statistics, phase relationships between limbs, and
% gait coordination measures.
%
% Syntax:
%   [walkfeaturestruct,perflywalkfeatures] = compute_WalkFeatures(obj,fly,digital_signal,condition_name)
%
% Input:
%   obj              - LimbBoutAnalyzer object
%   fly              - Integer index of the fly to analyze
%   digital_signal   - Logical array [1 x nframes] indicating frames to exclude from analysis
%                      (e.g., optogenetic stimulation periods)
%   condition_name   - String describing the experimental condition (currently unused in function)
%
% Output:
%   walkfeaturestruct    - Structure array with one element per valid walking bout containing:
%       .fly                      - Fly index
%       .walk_t0                  - Start frame of walking bout
%       .walk_t1                  - End frame of walking bout
%       .velmag_ctr               - Center velocity magnitude statistics during bout
%       .absdv_ctr                - Absolute sideways velocity change statistics
%       .absdu_ctr                - Absolute forward velocity change statistics
%       .absdtheta                - Absolute angular velocity change statistics
%       .left_vel                 - Left velocity statistics
%       .right_vel                - Right velocity statistics
%       .forward_vel              - Forward velocity statistics
%       .backward_vel             - Backward velocity statistics
%       .right_dtheta             - Right turning rate statistics
%       .left_dtheta              - Left turning rate statistics
%       .CoM_stability            - Center of mass stability statistics
%       .phaselag                 - Phase lag data relative to swing onset (method 1)
%       .phasediff_interp         - Phase differences using linear interpolation (method 2)
%       .phasediff_hilbert        - Phase differences using Hilbert transform on bout (method 3)
%       .phasediff_hilbert_global - Phase differences using pre-computed global Hilbert (method 4)
%
%   perflywalkfeatures   - Structure containing aggregated features across all walking bouts:
%       .fly                      - Array of fly indices for each bout
%       .walk_t0                  - Array of start frames for each bout
%       .walk_t1                  - Array of end frames for each bout
%       .(perframe_features)      - Aggregated statistics for each perframe feature
                                    % mean_* are the mean of means per walk bout
                                    % mean, std, n are computed on data (all walk bout frames)                                     
%       .(phase_features)         - Aggregated phase relationship data across bouts
%       
    
%
% Notes:
%   - Only walking bouts with at least 3 stance and 3 swing peaks per limb are included
%   - Phase relationships computed using 4 different methods for comparison
%   - Perframe features are computed as first derivatives unless in pfflist_none
%
% See also: LimbBoutAnalyzer, computePhaseLag, computeContinuousPhaseDiff_linearinterp, 
%           computeContinuousPhaseDiff_hilbert, detect_bouts


% summary plots
debug = 0;

% params for findpeaks
minPeakProminence = 1;

% list of perframe features to compute stats over bouts that are 
%  'first' = 1st derivative (d)  
pfflist_first = {'velmag_ctr','absdv_ctr','absdu_ctr','absdtheta', ...
    'left_vel','right_vel','forward_vel','backward_vel','right_dtheta','left_dtheta'};   
% 'none' = no derivative
pfflist_none = {'CoM_stability'};
% 'second' = second derivative (dd)



currflyboutdata = obj.limbBoutData(fly);
currfly_tips_pos_body = obj.tips_pos_body{fly};
nlimb = size(currfly_tips_pos_body,1);
walkfeaturestruct = struct;

% for hilbert global
currfly_tips_pos_body_Y = squeeze(currfly_tips_pos_body(:,2,:));
norm_ytips_global = zscore(currfly_tips_pos_body_Y')';
phases = nan(size(norm_ytips_global));
for limb = 1:nlimb
    phases(limb,:) = angle(hilbert(norm_ytips_global(limb,:)));
end



% compute frames for analysis
walk_digital = obj.walking_scores{fly};
% digital_signal = obj.digitalindicator{fly};
allowed_walks = walk_digital & ~digital_signal;


% loop over walks
[walk_t0s,walk_t1s] = detect_bouts(allowed_walks);
assert(numel(walk_t0s) == numel(walk_t1s))


% counter of processed walks
ct = 0;

for w = 1:numel(walk_t0s)    
    if isempty(walk_t0s) 
        continue
    end

    % general computations
    walk_t0 = walk_t0s(w);
    walk_t1 = walk_t1s(w);
    currwalk_tips_pos_body_Y = squeeze(currfly_tips_pos_body(:,2,walk_t0:walk_t1));
    % zscore tip data for peaks and hilbert
    norm_ytips = zscore(currwalk_tips_pos_body_Y');
    norm_ytips = norm_ytips';
    currwalk_norm_ytips_global = norm_ytips_global(:,walk_t0:walk_t1);
    currwalk_phases = phases(:,walk_t0:walk_t1);

    % find peaks in the Y pos signal
    if size(norm_ytips,2) < 3
        continue
    end
       
    for limb = 1:nlimb
        [ypeaks_t, loct] = findpeaks(norm_ytips(limb,:),'MinPeakProminence', minPeakProminence);
        loctall{limb} = [ypeaks_t;loct];
        numloct(limb) = numel(loct);
        [ypeaks_b, locb] = findpeaks((1-norm_ytips(limb,:)),'MinPeakProminence', minPeakProminence);
        locball{limb} = [ypeaks_b;locb];
        numlocb(limb) = numel(locb);
    end
    % filter out walks without enough cycles
    if all(numlocb >= 3) && all(numloct >= 3)
        ct = ct+1;
        walkfeaturestruct(ct).fly = fly;
        walkfeaturestruct(ct).walk_t0 = walk_t0;
        walkfeaturestruct(ct).walk_t1 = walk_t1;

        % compute perframe feature per bout  (speed, etc)
        % loop over list of perframe features
        for ifns = 1:numel(pfflist_first)
            fn = pfflist_first{ifns};
            datastruct = compute_StatsofPreframeFeatureDuringBouts(fly,fn,obj.trx,walk_t0,walk_t1,'first');
            walkfeaturestruct(ct).(fn) = datastruct;
        end
        for ifns = 1:numel(pfflist_none)
            fn = pfflist_none{ifns};
            datastruct = compute_StatsofPreframeFeatureDuringBouts(fly,fn,obj.trx,walk_t0,walk_t1,'none');
            walkfeaturestruct(ct).(fn) = datastruct;
        end


        % % method 1 - find phase lags relative to swing onset
        phaseoffsetdata = computePhaseLag(currflyboutdata, walk_t0, walk_t1);
        walkfeaturestruct(ct).phaselag = phaseoffsetdata;

        % % method 2 - find phase diff from peaks and linear interpolation Yang '23
        phasediff_interp = computeContinuousPhaseDiff_linearinterp(norm_ytips,loctall,locball,currwalk_tips_pos_body_Y,w,false);
        walkfeaturestruct(ct).phasediff_interp = phasediff_interp;

        % method 3 - compute hilbert for walk bout
        phasediff_hilbert = computeContinuousPhaseDiff_hilbert(norm_ytips,loctall,locball,currwalk_tips_pos_body_Y,w,debug);
        walkfeaturestruct(ct).phasediff_hilbert = phasediff_hilbert;

        % method 4 - use precomputed hilbert for whole movie - i don't
        % think zscore works as well over the whole movie so local phases
        % look better. 

        phasediff_hilbert_global = computeContinuousPhaseDiff_hilbert_global(currwalk_phases,currwalk_norm_ytips_global,loctall,locball, currwalk_tips_pos_body_Y,w,debug);
        walkfeaturestruct(ct).phasediff_hilbert_global = phasediff_hilbert_global;

        % compute TCS from methor 1 or 2, for each tripod
        % TO DO

    end

end


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





end
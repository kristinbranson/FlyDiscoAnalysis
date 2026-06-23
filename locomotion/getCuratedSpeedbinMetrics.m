function keys = getCuratedSpeedbinMetrics()
% GETCURATEDSPEEDBINMETRICS  LED-stripped keys of the curated metrics that
% receive speed-binned (__slow/__med/__fast) variants in locostatsperexp.
%
% Each key is 'feature__state__qualifier' WITHOUT the LED token, e.g.
% 'length_BL__step__pair1'. The combine* methods strip the LED token from a
% candidate field name and test membership here to decide whether to also
% emit binned variants. Only these metrics are binned; everything else keeps
% its single unbinned field (additive design).
%
% 70 base features: the VNC clustering 49-list MINUS step_direction (pairwise
% circular cancellation; parked for the turning analysis), PLUS the new
% pipeline metrics (duty_factor, TCS, gait_class fractions, P2A lags,
% nfeet_ground) and 8 phase-variability features. See project memory
% project_speed_binning_plan and the curated_metrics doc in
% locomotion_analysis/Alice.

pairs = {'pair1','pair2','pair3'};
add = @(f,s,q) cellfun(@(x) sprintf('%s__%s__%s',f,s,x), q, 'UniformOutput', false);

keys = {};

% --- step kinematics (per-step), x pair1/2/3 ---
stepfeats = {'AEP_BLx','AEP_BLy','PEP_BLx','PEP_BLy', ...
             'amplitude_BL','distance_BL','length_BL', ...
             'instataeous_frequency_steps','duty_factor'};
for i = 1:numel(stepfeats)
    keys = [keys, add(stepfeats{i},'step',pairs)]; %#ok<AGROW>
end

% --- bout durations + swing tip speed (per-bout), x pair1/2/3 ---
keys = [keys, add('durations_time','stance',pairs)];
keys = [keys, add('durations_time','swing',pairs)];
keys = [keys, add('mean_tips_speed_bodyref','swing',pairs)];

% --- whole-body (per-step, __all) ---
wholebody = {'velmag_ctr','absdtheta','absdu_ctr','absdv_ctr','CoM_stability'};
for i = 1:numel(wholebody)
    keys = [keys, add(wholebody{i},'step',{'all'})]; %#ok<AGROW>
end

% --- phase mean (per-frame): signed (circular) + abs (linear) ---
signed_groups = {'tripods_4','ipsi_P2A_4','ipsi_ant_2','ipsi_post_2'};
abs_groups    = {'abscontra_L2R_3','absRF_LF','absRM_LM','absRH_LH'};
keys = [keys, add('phasediff_hilbert','walk',signed_groups)];
keys = [keys, add('absphasediff_hilbert','walk',abs_groups)];

% --- phase variability (per-frame): circ_std (signed) / linear std (abs) ---
% (these fields are added to combinePhaseMetrics in a later piece)
keys = [keys, add('phasevar_hilbert','walk',signed_groups)];
keys = [keys, add('absphasevar_hilbert','walk',abs_groups)];

% --- gait class fractions (per-frame) ---
keys = [keys, add('gait_class','walk', ...
    {'tripod_frac','tetrapod_frac','grounded_frac','airborne_frac','other_frac'})];

% --- nfeet on ground (per-frame) ---
keys = [keys, add('nfeet_ground','walk',{'all'})];

% --- coordination scalars (per-event) ---
keys = [keys, add('TCS','walk',{'all'})];
keys = [keys, add('Pliftoff2Aliftoff_lag','walk',{'all','H_to_M','M_to_F'})];
keys = [keys, add('Ptouchdown2Aliftoff_lag','walk',{'all','H_to_M','M_to_F'})];

% sanity: 70 unique keys
assert(numel(keys) == 70, 'expected 70 curated keys, got %d', numel(keys));
assert(numel(unique(keys)) == numel(keys), 'duplicate curated keys');
end

function [boutfeatures] = computeStepFeatures(fly,trx,aptdata,tip_pos_body,legtip_landmarknums,limb,step_t0s,step_t1s,stance_t0s,stance_t1s,currfly_timestamps,stance_durations_time)
% tip_pos_body = 6 x 2 x T, data from tips_pos_body for 1 fly
% tip_pos_body are in ctrax format - need to use offset from trx
% assume start and end indices are stance only
% assume start and end indices are in movie reference frame
% compute stance_durations_time if not passed in
if nargin < 12 || isempty(stance_durations_time)
    [~, stance_durations_time] = computeBoutDurations(stance_t0s, stance_t1s, currfly_timestamps);
end

boutfeatures = struct;

boutfeatures.start_indices = step_t0s';
boutfeatures.end_indices = step_t1s';

% Body reference metrics
% offcurrfly = trx(fly).off;

% compute average body lengths in pixels (alt - compute from APT data)
meanbodylength = mean(trx(fly).a.*4);

% *step duration* in seconds 
% duration_s = currfly_timestamps(step_t1s)' - currfly_timestamps(step_t0s)';
% boutfeatures.duration_s = duration_s';

%compute the same way as swing and stance
[durations_frames,durations_time] = computeBoutDurations(step_t0s,step_t1s,currfly_timestamps);
boutfeatures.durations_frames = durations_frames;
boutfeatures.durations_time = durations_time; %(milliseconds)

% *duty factor* stance fraction of step cycle (Mendes 2013, Wosnitza 2013)
% duty_factor = stance_duration / step_duration; swing duty = 1 - duty_factor
% stance_durations_time already computed upstream in computeboutmetrics2.
% After restriction to walking bouts, step/stance arrays may lose alignment,
% so match each step to its stance by step_t0 == stance_t0.
nsteps_df = numel(step_t0s);
duty_factor = nan(size(durations_time));
for si = 1:nsteps_df
    midx = find(stance_t0s == step_t0s(si), 1);
    if ~isempty(midx) && durations_time(si) ~= 0
        duty_factor(si) = stance_durations_time(midx) / durations_time(si);
    end
end
boutfeatures.duty_factor = duty_factor;

% *step frequency* number of steps within a second
boutfeatures.instataeous_frequency_steps = 1./(boutfeatures.durations_time./1000);

% *step frequnecy* total
time_stepping = sum(durations_time)./1000;
nsteps = numel(durations_time);
boutfeatures.overall_frequency_steps = nsteps/time_stepping;


%%%% metrics bases on stance start and stance end; only computed for pairs within continguous walking bout
assert(numel(stance_t0s) == numel(stance_t1s));

% AEP - x,y at touch down (anterior extreme position) = stance_t0s, 2 x T
%   stance_t0s is INCLUSIVE (first in-contact frame); see limbSwingStance.m / detect_bouts.m
% PEP - x,y at lift off (posterior extreme position) = stance_t1s, 2 x T
%   velocity[t] = |pos[t+1]-pos[t]| (forward difference), so gc[t]=1 means the foot is
%   planted over the interval [t, t+1]. The foot's last in-contact interval is
%   [stance_t1s-1, stance_t1s] and it first moves over [stance_t1s, stance_t1s+1], so the
%   lift-off position is pos[stance_t1s]. (stance_t1s is the EXCLUSIVE bout end = first
%   swing index.) NB: pos[stance_t1s-1] is one frame too early (was the prior definition).

%TO DO check for flip in data AEP for limb 1 is -,+ instead of +,+

% *AEP* anterior extreme position
AEP = nan(2,numel(stance_t0s));
AEP(:,:) = tip_pos_body(limb,:,stance_t0s);
boutfeatures.AEP = AEP;
boutfeatures.AEP_BL = AEP./meanbodylength;

% *PEP* posterior extrene position
PEP = nan(2,numel(stance_t1s));
PEP(:,:) = tip_pos_body(limb,:,stance_t1s);
boutfeatures.PEP  = PEP;
boutfeatures.PEP_BL = PEP./meanbodylength;


% *step amplitude* (distance between PEP to AEP, Wosnitza et al. 2012, J Exp Biol) 
boutfeatures.amplitude_px = sqrt(sum((PEP' - AEP').^2, 2));
boutfeatures.amplitude_BL = boutfeatures.amplitude_px./meanbodylength;
boutfeatures.amplitude_px = boutfeatures.amplitude_px';
boutfeatures.amplitude_BL = boutfeatures.amplitude_BL';

% *step direction* Yang 2024 (angle of vector from AEP to PEP) 
boutfeatures.step_direction = atan2(PEP(2,:)-AEP(2,:),PEP(1,:)-AEP(1,:));

%%%% metrics based on step = AEP(1) to AEP(2); only computed for pairs within continguous walking bout
assert(numel(step_t0s) == numel(step_t1s));

nsteps = numel(step_t0s);

% *step distance* total distance leg travels AEP to AEP in body ref (Pratt '24)
distance_px = nan(1,nsteps);
for i = 1:nsteps
    curr_step = squeeze(tip_pos_body(limb,:,step_t0s(i):step_t1s(i)));
    dcurr_step = diff(curr_step,1,2)';
    curr_step_distances = hypot(dcurr_step(:,1), dcurr_step(:,2));
    distance_px(i) = sum(curr_step_distances);    
end
boutfeatures.distance_px = distance_px;
boutfeatures.distance_BL = boutfeatures.distance_px./meanbodylength;

% *step speed* step distance divide duration (Pratt '24) this seems weird -
% includes the 'stance' frames
boutfeatures.speed_pxpers = boutfeatures.distance_px./boutfeatures.durations_time./1000;
boutfeatures.speed_BLpers = boutfeatures.distance_BL./boutfeatures.durations_time./1000;


% Global reference metrics

% *step length* total distance leg travels AEP to AEP in global ref (Pratt '24)
aptdata = aptdata.pTrk;
currfly_aptlegdata = aptdata{fly}(legtip_landmarknums,:,:);
limb_data = squeeze(currfly_aptlegdata(limb,:,:));
length_px=nan(1,nsteps);
for i = 1:nsteps
    %only leg tips
    curr_step = limb_data(:,step_t0s(i):step_t1s(i));
    dcurr_step = diff(curr_step,1,2)';
    curr_step_distances = hypot(dcurr_step(:,1), dcurr_step(:,2));
    length_px(i) = sum(curr_step_distances);
end
boutfeatures.length_px= length_px;
boutfeatures.length_BL = length_px./meanbodylength;
boutfeatures.length_BL = boutfeatures.length_BL;

% add mean body velocity

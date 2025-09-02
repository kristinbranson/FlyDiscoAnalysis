function [boutfeatures] = computeStepFeatures(fly,trx,aptdata,tip_pos_body,legtip_landmarknums,limb,step_t0s,step_t1s,stance_t0s,stance_t1s,timestamps)
% tip_pos_body = 6 x 2 x T, data from tips_pos_body for 1 fly
% tip_pos_body are in ctrax format - need to use offset from trx
% assume start and end indices are stance only
% assume start and end indices are in movie reference frame
boutfeatures = struct; 

boutfeatures.start_indices = step_t0s';
boutfeatures.end_indices = step_t1s';

% Body reference metrics
% offcurrfly = trx(fly).off;

% compute average body lengths in pixels
meanbodylength = mean(trx(fly).a.*4);

% *step duration* in seconds
duration_s = timestamps(step_t1s)' - timestamps(step_t0s)';
boutfeatures.duration_s = duration_s';

% *step frequency* number of steps within a second
boutfeatures.instataeous_frequency_steps = 1./boutfeatures.duration_s;

% *step frequnecy* total
time_stepping = sum(duration_s);
nsteps = numel(duration_s);
boutfeatures.overall_frequency_steps = nsteps/time_stepping;


%%%% metrics bases on stance start and stance end; only computed for pairs within continguous walking bout
assert(numel(stance_t0s) == numel(stance_t1s));

% AEP - x,y at touch down - first frame of stance,t0s, 2 x T
% PEP - x,y at lift off - last frame of stance,t1s, 2 x T

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


% *step amplitude* (distance between PEP to AEP Wosnitza 2013 ) 
boutfeatures.amplitude_px = sqrt(sum((PEP' - AEP').^2, 2));
boutfeatures.amplitude_BL = boutfeatures.amplitude_px./meanbodylength;
boutfeatures.amplitude_px = boutfeatures.amplitude_px';
boutfeatures.amplitude_BL = boutfeatures.amplitude_BL';


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
boutfeatures.speed_pxpers = boutfeatures.distance_px./boutfeatures.duration_s;
boutfeatures.speed_BLpers = boutfeatures.distance_BL./boutfeatures.duration_s;


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

function [AEP, AEP_BL, PEP, PEP_BL, amplitude_px,amplitude_BL,distance_px,distance_BL,....
    length_px,length_BL,duration_s,frequency_steppers,speed_pxpers,speed_BLpers] = computeBoutSpatialFeatures(fly,trx,aptdata,tip_pos_body,limb,step_t0s,step_t1s,stance_t0s,stance_t1s,timestamps)
% tip_pos_body = 6 x 2 x T, data from tips_pos_body for 1 fly
% tip_pos_body are in ctrax format - need to use offset from trx
% assume start and end indices are stance only
% assume start and end indices are in movie reference frame

% Body reference metrics
offcurrfly = trx(fly).off;

% compute average body lengths
meanbodylength = mean(trx(fly).a.*4);


% AEP - x,y at touch down - first frame of stance, 2 x T
% PEP - x,y at lift off - last frame of stance (end-indices-1), 2 x T


%%%% metrics bases on stance start and stance end; only computed for pairs within continguous walking bout
assert(numel(stance_t0s) == numel(stance_t1s));

%TO DO check for flip in data AEP for limb 1 is -,+ instead of +,+

% *AEP* anterior extreme position 
AEP = nan(2,numel(stance_t0s));
AEP(:,:) = tip_pos_body(limb,:,stance_t0s+offcurrfly);
AEP_BL = AEP./meanbodylength;

% *PEP* posterior extrene position 
PEP = nan(2,numel(stance_t1s));
PEP(:,:) = tip_pos_body(limb,:,stance_t1s-1+offcurrfly);
PEP_BL = PEP./meanbodylength;


% *step amplitude* (distance between PEP to AEP Wosnitza 2013 ) 
amplitude_px = sqrt(sum((PEP' - AEP').^2, 2));
amplitude_BL = step_amplitude./meanbodylength;

%%%% metrics based on step = AEP(1) to AEP(2); only computed for pairs within continguous walking bout
assert(numel(step_t0s) == numel(step_t1s));

nsteps = numel(step_t0s);

% *step distance* total distance leg travels AEP to AEP in body ref (Pratt '24)
distance_px = nan(nsteps,1);
for i = 1:nsteps
    curr_step = squeeze(tip_pos_body(limb,:,step_t0s(i)+offcurrfly:step_t1s(i)+offcurrfly));
    dcurr_step = diff(curr_step');
    curr_step_distances = sqrt(dcurr_step(:,1).^2 + dcurr_step(:,2).^2);
    distance_px(i) = sum(curr_step_distances)';    
end
distance_BL = distance_px./meanbodylength;

% *step duration* in seconds
duration_s = timestamps(step_t1s)' - timestamps(step_t0s)';

% *step frequency* number of steps within a second
frequency_steppers = 1./duration_s;

% *step speed* step distance divide duration (Pratt '24) this seems weird -
% includes the 'stance' frames
speed_pxpers = distance_px./duration_s;
speed_BLpers = distance_BL./duration_s;


% Global reference metrics

% *step length* total distance leg travels AEP to AEP in global ref (Pratt '24)

length_px=nan(1,nsteps);
for i = 1:nsteps-1
    curr_step = squeeze(aptdata{fly}(limb,:,step_t0s(i)+offcurrfly:step_t1s(i)+offcurrfly));
    dcurr_step = diff(curr_step');
    curr_step_distances = sqrt(dcurr_step(:,1).^2 + dcurr_step(:,2).^2);
    length_px(i) = sum(curr_step_distances);
end
length_BL = length_px./meanbodylength;


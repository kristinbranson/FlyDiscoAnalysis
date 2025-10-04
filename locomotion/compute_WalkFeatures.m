function [walkfeaturestruct] = compute_WalkFeatures(obj,fly,digital_signal,condition_name)
% params
minPeakProminence = 1;
% summary plots
debug = 1;

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
digital_signal = obj.digitalindicator{fly};
allowed_walks = walk_digital & ~digital_signal;


% loop over walks
[walk_t0s,walk_t1s] = detect_bouts(allowed_walks);
assert(numel(walk_t0s) == numel(walk_t1s))

if isempty(walk_t0s)
    return
end

ct = 0;

for w = 1:20%:numel(walk_t0s)

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
        walkfeaturestruct(ct).walk_t0 = walk_t0;
        walkfeaturestruct(ct).walk_t1 = walk_t1;

        % compute perframe feature per bout  (speed, etc)


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


    end

end




end
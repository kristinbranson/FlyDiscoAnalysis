function [walkfeaturestruct] = compute_WalkFeatures(obj,fly,statename,condition_name)


currflyboutdata = obj.limbBoutData(fly);
currfly_tips_pos_body = obj.tips_pos_body{fly};
nlimb = size(currfly_tips_pos_body,1);
walkfeaturestruct = struct;
% loop over walks
[walk_t0s,walk_t1s] = detect_bouts(obj.walking_scores{fly});
assert(numel(walk_t0s) == numel(walk_t1s))

if isempty(walk_t0s)
    return
end


for w = 1:20;%numel(walk_t0s)

   % compute perframe feature per bout  (speed, etc)

   % % method 1 - find phase lags relative to swing onset
    walk_t0 = walk_t0s(w);
    walk_t1 = walk_t1s(w); 

    phaseoffsetdata = computePhaseLag(currflyboutdata, walk_t0, walk_t1);
    walkfeaturestruct(w).phaselag = phaseoffsetdata;   
    
    
    % zscore tip data for peaks and hilbert

    currwalk_tips_pos_body_Y = squeeze(currfly_tips_pos_body(:,2,walk_t0:walk_t1));
    
    [norm_ytips,mu,sigma] = zscore(currwalk_tips_pos_body_Y');
    norm_ytips = norm_ytips';


    % % method 2 - find phase diff from peaks and linear interpolation Yang '23

    phasediff_interp = computeContinuousPhaseDiff_linearinterp(norm_ytips,currwalk_tips_pos_body_Y,w);
    walkfeaturestruct(w).phasediff_interp = phasediff_interp;
    
    % method 3 - compute hilbert for walk bout

    phasediff_hilbert = computeContinuousPhaseDiff_hilbert(norm_ytips,currwalk_tips_pos_body_Y,w);
    walkfeaturestruct(w).phasediff_hilbert = phasediff_hilbert;
 
    % method 4 - use precomputed hilbert for whole movie

    % compute TCS from methor 1 or 2, for each tripod

    


end


 

end
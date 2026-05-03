function data = compute_gait_class(groundcontact)
% compute_gait_class - Per-frame hexapod gait classification.
%
% For each frame, classify the 6-leg ground-contact pattern against
% canonical templates (tripod / tetrapod / grounded / airborne / other).
% See gait_pattern_constants.m for definitions.
%
% Input:
%   groundcontact - {1 x nflies} cell. Each entry is a 6 x nframes
%                   logical/numeric matrix of stance state (1 = stance,
%                   0 = swing). Row order: APT keypoints 12-17 =
%                   (RF, RM, RH, LH, LM, LF), matching tip_pos_body.
%
% Output:
%   data - {1 x nflies} cell. Each entry is a 1 x nframes uint8 array
%          of class codes (1 tripod, 2 tetrapod, 3 grounded, 4 airborne,
%          5 other). Stored as cell-of-rows to match the per-frame
%          feature convention used by other locomotion .mat files
%          (e.g. nfeet_ground.mat).

k = gait_pattern_constants();
nflies = numel(groundcontact);
data = cell(1, nflies);

for fly = 1:nflies
    gc = groundcontact{fly};
    if isempty(gc)
        data{fly} = uint8([]);
        continue;
    end
    % Reorder rows into Mendes (L1 L2 L3 R1 R2 R3) layout, then
    % transpose so each row is one frame's pattern.
    pat = double(gc(k.legorder_to_mendes, :))';   % nframes x 6
    nframes = size(pat, 1);

    classes = repmat(uint8(k.code.other), 1, nframes);

    is_tripod   = ismember(pat, k.tripod,   'rows');
    is_tetrapod = ismember(pat, k.tetrapod, 'rows');
    is_grounded = ismember(pat, k.grounded, 'rows');
    is_airborne = ismember(pat, k.airborne, 'rows');

    classes(is_tripod)   = k.code.tripod;
    classes(is_tetrapod) = k.code.tetrapod;
    classes(is_grounded) = k.code.grounded;
    classes(is_airborne) = k.code.airborne;

    data{fly} = classes;
end

end

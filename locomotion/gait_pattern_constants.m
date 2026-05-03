function k = gait_pattern_constants()
% gait_pattern_constants - Canonical hexapod gait stance patterns.
%
% Per-frame gait classification: each frame's 6-leg ground-contact
% pattern is matched against canonical templates (tripod, tetrapod,
% all-grounded, all-airborne); anything else is "other".
%
% Canonical templates: Wilson 1966 (Annu Rev Entomol 11:103-22,
% originating framework) / Collins & Stewart 1993 (Biol Cybern 68:287-98,
% CPG-symmetry derivation of the same templates).
% Per-frame matching approach: Mendes 2013 (eLife 2:e00231),
% Wosnitza 2012 (J Exp Biol 216(3):480-91), DeAngelis 2019 (eLife 8:e46409).
%
% Class codes (used as integer values in the gait_class perframe feature):
%   1 = tripod
%   2 = tetrapod
%   3 = grounded   (all 6 legs in stance)
%   4 = airborne   (no legs in stance)
%   5 = other
%
% Pattern leg order: Mendes/DeAngelis convention L1 L2 L3 R1 R2 R3
% = (LF, LM, LH, RF, RM, RH). 1 = stance, 0 = swing.
%
% The natural row order of tip_pos_body / groundcontact is APT
% keypoints 12-17 = (RF, RM, RH, LH, LM, LF). The permutation
% k.legorder_to_mendes = [6 5 4 1 2 3] reorders tip_pos_body rows
% into Mendes order. compute_gait_class.m applies it once before
% pattern matching.
%
% Tetrapod patterns: 6 total — 2 cycles of 3 swing-pair phases each
% (left-tetrapod and right-tetrapod, per DeAngelis Table 2). Two of
% these (swing pairs (L1,R3) and (L3,R1)) were missing from Alice's
% earlier compute_gaitclassification.m, so old gait analyses
% under-counted tetrapod and over-counted "other".

k.code.tripod   = 1;
k.code.tetrapod = 2;
k.code.grounded = 3;
k.code.airborne = 4;
k.code.other    = 5;

k.names = {'tripod', 'tetrapod', 'grounded', 'airborne', 'other'};

% Permutation: groundcontact rows (RF RM RH LH LM LF) -> Mendes (L1 L2 L3 R1 R2 R3)
k.legorder_to_mendes = [6 5 4 1 2 3];

% Tripod: 3 legs in swing forming a "triangle of support".
% Tripod A swing = {L1, L3, R2} -> stance pattern [0 1 0 1 0 1]
% Tripod B swing = {L2, R1, R3} -> stance pattern [1 0 1 0 1 0]
k.tripod = [0 1 0 1 0 1;   % swing (L1, L3, R2)
            1 0 1 0 1 0];  % swing (L2, R1, R3)

% Tetrapod: 2 legs in swing simultaneously. 6 canonical patterns total.
% Right-tetrapod cycle (3 phases):
%   swing (L2, R1), (L3, R2), (L1, R3)
% Left-tetrapod cycle (3 phases):
%   swing (L1, R2), (L2, R3), (L3, R1)
k.tetrapod = [1 0 1 0 1 1;   % swing (L2, R1)
              1 1 0 1 0 1;   % swing (L3, R2)
              0 1 1 1 1 0;   % swing (L1, R3)
              0 1 1 1 0 1;   % swing (L1, R2)
              1 0 1 1 1 0;   % swing (L2, R3)
              1 1 0 0 1 1];  % swing (L3, R1)

k.grounded = [1 1 1 1 1 1];
k.airborne = [0 0 0 0 0 0];

end

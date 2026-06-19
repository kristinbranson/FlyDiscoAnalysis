%% test_find_bout_overlap.m
% Unit test for the bout-type-aware overlap filter find_bout_overlap.m.
%
% Convention (limbSwingStanceStep.m / detect_bouts.m):
%   swing/stance end_indices are EXCLUSIVE (end = first frame AFTER the bout).
%   step end_indices are INCLUSIVE-ish (end = next touchdown, a real frame).
% find_bout_overlap keeps a bout only if every real frame of the bout is ON.
% The end_exclusive flag selects start:end-1 (swing/stance) vs start:end (step).
% Returned end indices are ALWAYS the original inputs, unchanged.

addpath('/groups/branson/home/robiea/Code_versioned/FlyDiscoAnalysis'); modpath;

fprintf('=== test_find_bout_overlap ===\n\n');
npass = 0; nfail = 0;

% ON window: frames 1..8 ON, 9..12 OFF.  Last ON frame = 8.
%   frame:  1 2 3 4 5 6 7 8 9 10 11 12
signal = [1 1 1 1 1 1 1 1 0 0 0 0];

%% --- Test 1: swing/stance bout ending exactly at the ON->OFF boundary ----
% A stance bout occupying real frames 5,6,7,8 has exclusive end = 9
% (first frame after).  Its real frames (5:8) are ALL ON -> should be KEPT.
% The OLD code tested signal(5:9), pulled in OFF frame 9, and dropped it.
start_in = 5; end_in = 9;   % exclusive end
[os, oe] = find_bout_overlap(signal, start_in, end_in, true);
if isequal(os(:),5) && isequal(oe(:),9)
    fprintf('Test 1 PASS: boundary swing/stance bout kept; end index returned unchanged (9).\n');
    npass = npass + 1;
else
    fprintf('Test 1 FAIL: got start=%s end=%s, expected 5 / 9.\n', mat2str(os(:)'), mat2str(oe(:)'));
    nfail = nfail + 1;
end

%% --- Test 2: same bout, but actually straddling the boundary -> DROPPED ---
% Real frames 6,7,8,9 (exclusive end = 10): frame 9 is OFF -> must be dropped.
[os, oe] = find_bout_overlap(signal, 6, 10, true);
if isempty(os) && isempty(oe)
    fprintf('Test 2 PASS: straddling swing/stance bout correctly dropped.\n');
    npass = npass + 1;
else
    fprintf('Test 2 FAIL: got start=%s end=%s, expected empty.\n', mat2str(os(:)'), mat2str(oe(:)'));
    nfail = nfail + 1;
end

%% --- Test 3: step bout whose end IS the last ON frame -> KEPT (inclusive) --
% A step ends at the next touchdown (a real frame).  Step real frames 5..8,
% end = 8 (the touchdown).  With end_exclusive=false, test signal(5:8) -> all
% ON -> kept.  (A blanket end-1 would wrongly test 5:7 and mis-handle step.)
[os, oe] = find_bout_overlap(signal, 5, 8, false);
if isequal(os(:),5) && isequal(oe(:),8)
    fprintf('Test 3 PASS: step bout ending on last ON frame kept (inclusive end).\n');
    npass = npass + 1;
else
    fprintf('Test 3 FAIL: got start=%s end=%s, expected 5 / 8.\n', mat2str(os(:)'), mat2str(oe(:)'));
    nfail = nfail + 1;
end

%% --- Test 4: step bout whose touchdown lands on the first OFF frame -> DROP -
% Step real frames 6..9, end = 9 (touchdown on OFF frame) -> dropped.
[os, oe] = find_bout_overlap(signal, 6, 9, false);
if isempty(os) && isempty(oe)
    fprintf('Test 4 PASS: step bout with touchdown on OFF frame correctly dropped.\n');
    npass = npass + 1;
else
    fprintf('Test 4 FAIL: got start=%s end=%s, expected empty.\n', mat2str(os(:)'), mat2str(oe(:)'));
    nfail = nfail + 1;
end

%% --- Test 5: default (no flag) == inclusive (backwards compatible) ---------
% Old callers passed no flag and got start:end behavior; default must match.
[os_def, oe_def] = find_bout_overlap(signal, 5, 8);
[os_exp, oe_exp] = find_bout_overlap(signal, 5, 8, false);
if isequal(os_def, os_exp) && isequal(oe_def, oe_exp) && isequal(os_def(:),5)
    fprintf('Test 5 PASS: default flag reproduces historical start:end behavior.\n');
    npass = npass + 1;
else
    fprintf('Test 5 FAIL: default behavior changed.\n');
    nfail = nfail + 1;
end

%% --- Test 6: multiple bouts, mixed, returned ends preserve convention ------
% Two swing/stance bouts: [2,4) interior (real 2,3 ON -> keep),
% [7,9) boundary (real 7,8 ON -> keep), [8,10) straddle (frame 9 OFF -> drop).
starts = [2; 7; 8];
ends   = [4; 9; 10];
[os, oe] = find_bout_overlap(signal, starts, ends, true);
if isequal(os(:),[2;7]) && isequal(oe(:),[4;9])
    fprintf('Test 6 PASS: mixed swing/stance bouts filtered; exclusive ends preserved.\n');
    npass = npass + 1;
else
    fprintf('Test 6 FAIL: got start=%s end=%s, expected [2 7] / [4 9].\n', ...
        mat2str(os(:)'), mat2str(oe(:)'));
    nfail = nfail + 1;
end

fprintf('\n=== SUMMARY: %d passed, %d failed ===\n', npass, nfail);

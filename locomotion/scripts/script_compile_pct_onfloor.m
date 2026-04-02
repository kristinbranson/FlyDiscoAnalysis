%% script_compile_pct_onfloor
% Compile per-experiment pct_onfloor for all VNC2+VNC3 experiments.
% onfloor = ~(onceiling | nottracking) — union of both classifiers.
%
% Input:
%   expdirs_VNC23_20251011.mat — metadata with experiment paths
%
% Output:
%   pct_onfloor_VNC23.mat containing:
%     expdirs         — cell array of full paths (5525)
%     expnames        — cell array of experiment directory names
%     pct_onfloor     — per-experiment % of frames onfloor (pooled across flies)
%     pct_onceiling   — per-experiment % onceiling only (diagnostic)
%     pct_nottracking — per-experiment % nottracking only (diagnostic)
%     screen_type     — 'VNC2' or 'VNC3'
%     is_control      — logical, true for YNA_K_162984
%     is_poscontrol   — logical, true for EXT_VGLUT-GAL4

outfile = '/groups/branson/bransonlab/flydisco_linelevel_VNC/pct_onfloor_VNC23.mat';

%% Load experiment list
M = load('/groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/expdirs_VNC23_20251011.mat');
expdirs = {M.metadata.file_system_path}';
n = numel(expdirs);
fprintf('Processing %d experiments...\n', n);

%% Preallocate
pct_onfloor = nan(n, 1);
pct_onceiling = nan(n, 1);
pct_nottracking = nan(n, 1);
expnames = cell(n, 1);
screen_type = cell(n, 1);
is_control = false(n, 1);
is_poscontrol = false(n, 1);
issuccess = true(n, 1);

%% Loop over experiments
parfor i = 1:n
    [~, ename] = fileparts(expdirs{i});
    expnames{i} = ename;

    % Screen type
    if startsWith(ename, 'VNC3_'), screen_type{i} = 'VNC3';
    elseif startsWith(ename, 'VNC2_'), screen_type{i} = 'VNC2';
    else, screen_type{i} = 'unknown';
    end

    % Control and positive control lines
    is_control(i) = contains(ename, 'YNA_K_162984');
    is_poscontrol(i) = contains(ename, 'EXT_VGLUT-GAL4');

    try
        % Load both score files
        OC = load(fullfile(expdirs{i}, 'scores_onceiling_resnet_v2.mat'), 'allScores');
        NT = load(fullfile(expdirs{i}, 'scores_nottracking.mat'), 'allScores');

        nflies = numel(OC.allScores.postprocessed);
        total_frames = 0;
        n_onceiling = 0;
        n_nottracking = 0;
        n_onfloor = 0;

        for f = 1:nflies
            oc_pp = OC.allScores.postprocessed{f};
            nt_pp = NT.allScores.postprocessed{f};
            nf = numel(oc_pp);
            total_frames = total_frames + nf;
            n_onceiling = n_onceiling + sum(oc_pp == 1);
            n_nottracking = n_nottracking + sum(nt_pp == 1);
            n_onfloor = n_onfloor + sum(~(oc_pp == 1 | nt_pp == 1));
        end

        pct_onceiling(i) = 100 * n_onceiling / total_frames;
        pct_nottracking(i) = 100 * n_nottracking / total_frames;
        pct_onfloor(i) = 100 * n_onfloor / total_frames;
    catch ME
        fprintf('  Error on %s: %s\n', ename, ME.message);
        issuccess(i) = false;
    end

    if mod(i, 500) == 0, fprintf('  %d/%d\n', i, n); end
end

fprintf('Done. %d/%d succeeded.\n', sum(issuccess), n);

%% Save
save(outfile, 'expdirs', 'expnames', 'pct_onfloor', 'pct_onceiling', ...
    'pct_nottracking', 'screen_type', 'is_control', 'is_poscontrol', ...
    'issuccess', '-v7.3');
fprintf('\nSaved to %s\n', outfile);

%% Load (if not already in memory)
if ~exist('pct_onfloor', 'var')
    outfile = '/groups/branson/bransonlab/flydisco_linelevel_VNC/pct_onfloor_VNC23.mat';
    load(outfile);
    fprintf('Loaded from %s\n', outfile);
end

%% Summary: experiments to exclude (pct_onfloor < 10%)
idx_nonvglut = ~is_poscontrol & issuccess;
idx_vglut = is_poscontrol & issuccess;
n_nonvglut = sum(idx_nonvglut);
n_vglut = sum(idx_vglut);

fprintf('\n=== Exclusion summary (pct_onfloor < 10%% threshold) ===\n');
fprintf('\nNon-VGLUT (%d total):\n', n_nonvglut);
fprintf('  Excluded (pct_onfloor < 10%%): %d/%d (%.1f%%)\n', ...
    sum(pct_onfloor(idx_nonvglut) < 10), n_nonvglut, ...
    100*sum(pct_onfloor(idx_nonvglut) < 10)/n_nonvglut);
fprintf('  Kept (pct_onfloor >= 10%%):    %d/%d (%.1f%%)\n', ...
    sum(pct_onfloor(idx_nonvglut) >= 10), n_nonvglut, ...
    100*sum(pct_onfloor(idx_nonvglut) >= 10)/n_nonvglut);

fprintf('\nVGLUT (%d total):\n', n_vglut);
fprintf('  Excluded (pct_onfloor < 10%%): %d/%d (%.1f%%)\n', ...
    sum(pct_onfloor(idx_vglut) < 10), n_vglut, ...
    100*sum(pct_onfloor(idx_vglut) < 10)/n_vglut);
fprintf('  Kept (pct_onfloor >= 10%%):    %d/%d (%.1f%%)\n', ...
    sum(pct_onfloor(idx_vglut) >= 10), n_vglut, ...
    100*sum(pct_onfloor(idx_vglut) >= 10)/n_vglut);

fprintf('\nAll experiments (%d total):\n', n_nonvglut + n_vglut);
fprintf('  Excluded: %d/%d\n', ...
    sum(pct_onfloor(issuccess) < 10), sum(issuccess));
fprintf('  Kept:     %d/%d\n', ...
    sum(pct_onfloor(issuccess) >= 10), sum(issuccess));

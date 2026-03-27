%% script_compile_walk_onceiling_overlap
% Compile per-fly walk/onceiling overlap for VNC2+VNC3 experiments
% with 25-90% onceiling (excluding VGLUT).
%
% For each fly computes:
%   pct_oc         — % of frames onceiling
%   pct_walk       — % of frames walking (any surface)
%   pct_walk_floor — % of frames walking AND NOT onceiling (salvageable)
%   pct_walk_ceil  — % of frames walking AND onceiling (to discard)
%
% Output:
%   walk_onceiling_overlap.mat containing:
%     exp_idx       — index into allVNC summary for each experiment in this subset
%     exp_names     — experiment names
%     exp_pct_oc    — per-experiment overall % onceiling
%     exp_screen_type — VNC2 or VNC3
%     fly_data      — struct array with per-fly stats

summaryfile = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/allVNC_onceiling_resnet_v2_summary.mat';
outfile = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/walk_onceiling_overlap_vnc23_25to90.mat';

load(summaryfile, 'pct_oc', 'screen_type', 'expnames', 'expdirs');

is_vglut = contains(expnames, 'VGLUT');
is_vnc23 = ismember(screen_type, {'VNC2', 'VNC3'}) & ~is_vglut;
in_range = is_vnc23 & pct_oc >= 25 & pct_oc <= 90;

exp_idx = find(in_range);
n_exp = numel(exp_idx);
fprintf('Compiling %d experiments...\n', n_exp);

exp_names = expnames(exp_idx);
exp_pct_oc = pct_oc(exp_idx);
exp_screen_type = screen_type(exp_idx);

% Pre-allocate fly_data as struct array
fly_data = struct('expname', {}, 'exp_i', {}, 'fly', {}, 'sex', {}, ...
    'nframes', {}, 'pct_oc', {}, 'pct_walk', {}, ...
    'pct_walk_floor', {}, 'pct_walk_ceil', {});

for i = 1:n_exp
    edir = expdirs{exp_idx(i)};
    ename = exp_names{i};

    soc = load(fullfile(edir, 'scores_onceiling_resnet_v2.mat'));
    swk = load(fullfile(edir, 'scores_Walk2.mat'));

    nflies = numel(soc.allScores.scores);

    for f = 1:nflies
        oc = soc.allScores.scores{f} > 0;
        wk = swk.allScores.scores{f} > 0;
        nf = min(numel(oc), numel(wk));
        oc = oc(1:nf);
        wk = wk(1:nf);

        if f <= 5, sex = 'M'; else, sex = 'F'; end

        entry.expname = ename;
        entry.exp_i = i;
        entry.fly = f;
        entry.sex = sex;
        entry.nframes = nf;
        entry.pct_oc = 100 * mean(oc);
        entry.pct_walk = 100 * mean(wk);
        entry.pct_walk_floor = 100 * mean(wk & ~oc);
        entry.pct_walk_ceil = 100 * mean(wk & oc);

        fly_data(end+1) = entry; %#ok<SAGROW>
    end

    if mod(i, 200) == 0, fprintf('  %d/%d\n', i, n_exp); end
end
fprintf('Done. %d flies across %d experiments.\n', numel(fly_data), n_exp);

save(outfile, 'exp_idx', 'exp_names', 'exp_pct_oc', 'exp_screen_type', ...
    'fly_data', '-v7.3');
fprintf('Saved to %s\n', outfile);

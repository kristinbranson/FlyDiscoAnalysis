%% script_explore_walk_classifier_overlap.m
% Compute per-walk frac_notonfloor (union of onceiling | nottracking)
% for experiments with moderate ceiling activity (25-90% median on-ceiling).
% Plot histogram of frac_notonfloor to inform per-walk filtering threshold.

modpath;

%% Config
rootdir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20260225_VNC23controldata4onceilingwalkexplore';
savedir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Locomotion_analysis/claudeplots_20260324_walk_classifier_filtering';
cache_file = fullfile(rootdir, 'walk_classifier_overlap_cache_allexp.mat');

if ~exist(savedir, 'dir'), mkdir(savedir); end

%% Section 1: Find all experiment directories
dd = dir(fullfile(rootdir, 'VNC*'));
dd = dd([dd.isdir]);
sel_expnames = {dd.name};
nexp = numel(sel_expnames);
fprintf('Found %d experiment directories\n', nexp);

%% Section 2: Compute per-walk frac_notonfloor
if exist(cache_file, 'file')
    fprintf('Loading cached walk data from %s\n', cache_file);
    C = load(cache_file);
    all_frac_notonfloor = C.all_frac_notonfloor;
    all_velmag = C.all_velmag;
    all_walk_duration = C.all_walk_duration;
    all_fly = C.all_fly;
    all_exp_idx = C.all_exp_idx;
    all_exp_median_oc = C.all_exp_median_oc;
    fprintf('Loaded %d walks from cache\n', numel(all_frac_notonfloor));
else
    fprintf('Computing per-walk frac_notonfloor for %d experiments...\n', nexp);

    % Pre-allocate with cell arrays, concatenate at end
    frac_notonfloor_cell = cell(nexp, 1);
    velmag_cell = cell(nexp, 1);
    walk_duration_cell = cell(nexp, 1);
    fly_cell = cell(nexp, 1);
    exp_idx_cell = cell(nexp, 1);
    for i = 1:nexp
        if mod(i, 50) == 0
            fprintf('  %d / %d experiments\n', i, nexp);
        end

        expdir = fullfile(rootdir, sel_expnames{i});

        % Load walk struct
        ws_file = fullfile(expdir, 'locomotion_walkstruct.mat');
        if ~exist(ws_file, 'file')
            warning('Missing walkstruct for %s', sel_expnames{i});
            continue;
        end
        W = load(ws_file, 'walk_struct_OFF');
        ws = W.walk_struct_OFF;
        nwalks = numel(ws.fly);
        if nwalks == 0, continue; end

        % Load classifier scores
        oc_file = fullfile(expdir, 'scores_onceiling_resnet_v2.mat');
        nt_file = fullfile(expdir, 'scores_nottracking_apt.mat');
        if ~exist(oc_file, 'file')
            warning('Missing onceiling scores for %s', sel_expnames{i});
            continue;
        end
        if ~exist(nt_file, 'file')
            warning('Missing nottracking scores for %s', sel_expnames{i});
            continue;
        end
        OC = load(oc_file, 'allScores');
        NT = load(nt_file, 'allScores');

        % Compute frac_notonfloor per walk
        frac = nan(1, nwalks);
        for w = 1:nwalks
            fly_id = ws.fly(w);
            t0 = ws.walk_t0(w);
            t1 = ws.walk_t1(w);
            walk_frames = t0:t1;

            pp_oc = OC.allScores.postprocessed{fly_id};
            pp_nt = NT.allScores.postprocessed{fly_id};

            % Frame bounds check
            valid = walk_frames >= 1 & walk_frames <= min(numel(pp_oc), numel(pp_nt));
            if sum(valid) == 0, continue; end

            vf = walk_frames(valid);
            oc_vals = pp_oc(vf);
            nt_vals = pp_nt(vf);
            % Treat NaN as not-flagged (0) for the OR
            oc_vals(isnan(oc_vals)) = 0;
            nt_vals(isnan(nt_vals)) = 0;
            frac(w) = mean(oc_vals > 0 | nt_vals > 0);
        end

        frac_notonfloor_cell{i} = frac;
        velmag_cell{i} = ws.velmag_ctr;
        walk_duration_cell{i} = ws.walk_duration;
        fly_cell{i} = ws.fly;
        exp_idx_cell{i} = repmat(i, 1, nwalks);
    end

    all_frac_notonfloor = [frac_notonfloor_cell{:}];
    all_velmag = [velmag_cell{:}];
    all_walk_duration = [walk_duration_cell{:}];
    all_fly = [fly_cell{:}];
    all_exp_idx = [exp_idx_cell{:}];

    fprintf('Computed frac_notonfloor for %d walks across %d experiments\n', ...
        numel(all_frac_notonfloor), nexp);

    save(cache_file, 'all_frac_notonfloor', 'all_velmag', 'all_walk_duration', ...
        'all_fly', 'all_exp_idx', 'sel_expnames');
    fprintf('Saved cache to %s\n', cache_file);
end

%% Section 3: Plot histogram of frac_notonfloor
valid_walks = ~isnan(all_frac_notonfloor);
frac = all_frac_notonfloor(valid_walks);
nwalks_total = numel(frac);

fprintf('\nSummary:\n');
fprintf('  Total walks: %d\n', nwalks_total);
fprintf('  frac_notonfloor == 0: %d (%.1f%%)\n', sum(frac == 0), 100*mean(frac == 0));
fprintf('  frac_notonfloor > 0: %d (%.1f%%)\n', sum(frac > 0), 100*mean(frac > 0));
fprintf('  frac_notonfloor > 0.5: %d (%.1f%%)\n', sum(frac > 0.5), 100*mean(frac > 0.5));
fprintf('  frac_notonfloor == 1: %d (%.1f%%)\n', sum(frac == 1), 100*mean(frac == 1));

edges = 0:0.02:1;

fig = figure('Position', [100 100 900 500]);

histogram(frac, edges, 'Normalization', 'count', ...
    'FaceColor', [0.3 0.5 0.8], 'EdgeColor', 'w');
set(gca, 'YScale', 'log');
xlabel('frac\_notonfloor (fraction of walk frames: onceiling | nottracking)', 'Interpreter', 'none');
ylabel('Number of walks (log scale)');
title(sprintf('Per-walk frac\\_notonfloor distribution\n%d walks from %d experiments (all VNC/VNC2/VNC3)', ...
    nwalks_total, nexp), 'Interpreter', 'none');
grid on;

% Add text annotation with key percentiles
pcts = [0 0.01 0.05 0.1 0.25 0.5 0.75 1.0];
for p = pcts
    n_above = sum(frac > p);
    pct_above = 100 * n_above / nwalks_total;
    if p == 0
        text_str = sprintf('>%.0f%%: %d walks (%.1f%%)', p*100, n_above, pct_above);
    else
        text_str = sprintf('>%.0f%%: %d (%.1f%%)', p*100, n_above, pct_above);
    end
    fprintf('  %s\n', text_str);
end

saveas(fig, fullfile(savedir, 'histogram_frac_notonfloor.png'));
fprintf('\nSaved histogram to %s\n', fullfile(savedir, 'histogram_frac_notonfloor.png'));

% Linear scale version
fig2 = figure('Position', [100 100 900 500]);
histogram(frac, edges, 'Normalization', 'count', ...
    'FaceColor', [0.3 0.5 0.8], 'EdgeColor', 'w');
xlabel('frac\_notonfloor (fraction of walk frames: onceiling | nottracking)', 'Interpreter', 'none');
ylabel('Number of walks');
title(sprintf('Per-walk frac\\_notonfloor distribution (linear scale)\n%d walks from %d experiments (all VNC/VNC2/VNC3)', ...
    nwalks_total, nexp), 'Interpreter', 'none');
grid on;

saveas(fig2, fullfile(savedir, 'histogram_frac_notonfloor_linear.png'));
fprintf('Saved linear histogram to %s\n', fullfile(savedir, 'histogram_frac_notonfloor_linear.png'));

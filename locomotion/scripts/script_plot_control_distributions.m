%% script_plot_control_distributions.m
% Plot distributions of all 1498 locomotion statistics across control
% experiments (YNA_K_162984) under 4 conditions:
%   1. Raw means
%   2. Ceiling-filtered means
%   3. Normalized means
%   4. Ceiling-filtered + normalized means
%
% Each stat produces a 3-panel figure:
%   Left   = per-experiment histogram (pdf-normalized, lillietest)
%   Middle = per-set histogram (pdf-normalized, lillietest)
%   Right  = bootstrap null distribution (pdf-normalized, bootstrap diagnostics)
%
% Histograms are pdf-normalized so the 3 panels are directly comparable.
% Experiment and set panels use lillietest (linear) or Rayleigh (circular)
% plus Hartigan's dip test for multimodality (linear stats only).
% Bootstrap panel shows diagnostics for p-value reliability:
%   skewness, excess kurtosis, n_ctrl_sets, tail_ratio (observed/expected |z|>3)
%
% Outputs:
%   - PNGs in 4 subdirectories
%   - normality_pvalues.csv with lillietest/Rayleigh + dip test p-values for exp and set panels
%   - bootstrap_diagnostics.csv with skewness, kurtosis, n_sets, tail_ratio per condition

modpath;

%% Config
collected_stats_file = '/groups/branson/bransonlab/flydisco_linelevel_VNC/CollectedVNC23PerFrameStats20260205.mat';
ceiling_data_dir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20260225_VNC23controldata4onceilingwalkexplore';
outdir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Locomotion_analysis/claudeplots_20260304_control_distributions';
cache_file = fullfile(outdir, 'ceiling_filter_cache.mat');

ceiling_fly_thresh = 0.5;   % fraction of frames on ceiling to call a fly "ceiling fly"
ceiling_exp_thresh = 0.5;   % fraction of flies that are ceiling flies to exclude experiment
nbins = 30;
nsamples = 10000;            % bootstrap samples for null distribution
k_sets = 2;                  % median number of sets per line
dip_nboot = 200;             % bootstrap samples for Hartigan's dip significance test

%% Load collected stats
fprintf('Loading collected stats...\n');
S = load(collected_stats_file, 'expstats', 'setstats', 'exp2lineidx', ...
    'set2lineidx', 'idxcontrol', 'statfns', 'iscircstat', 'setidx', 'expdirs');

nstats = numel(S.statfns);
ctrl_exp_idx = find(S.exp2lineidx == S.idxcontrol);
ctrl_set_idx = find(S.set2lineidx == S.idxcontrol);
nctrl_exp = numel(ctrl_exp_idx);
nctrl_set = numel(ctrl_set_idx);
fprintf('Control experiments: %d, control sets: %d\n', nctrl_exp, nctrl_set);

%% Create output directories
subdirs = {'raw', 'ceiling_filtered', 'normalized', 'ceiling_filtered_normalized'};
for i = 1:numel(subdirs)
    d = fullfile(outdir, subdirs{i});
    if ~exist(d, 'dir'), mkdir(d); end
end

%% Compute ceiling filter (with cache)
if exist(cache_file, 'file')
    fprintf('Loading cached ceiling filter from %s\n', cache_file);
    C = load(cache_file, 'exclude_exp', 'frac_ceiling_per_exp', 'n_excluded');
    exclude_exp = C.exclude_exp;
    frac_ceiling_per_exp = C.frac_ceiling_per_exp;
    fprintf('Cached: %d experiments excluded\n', C.n_excluded);
else
    fprintf('Computing ceiling filter for %d control experiments...\n', nctrl_exp);
    exclude_exp = false(1, nctrl_exp);
    frac_ceiling_per_exp = nan(1, nctrl_exp);

    for i = 1:nctrl_exp
        [~, expname] = fileparts(S.expdirs{ctrl_exp_idx(i)});
        score_file = fullfile(ceiling_data_dir, expname, 'scores_onceiling.mat');
        if ~exist(score_file, 'file')
            warning('Missing ceiling scores for %s', expname);
            continue;
        end
        Sc = load(score_file, 'allScores');
        nflies = numel(Sc.allScores.postprocessed);
        n_ceiling_flies = 0;
        for f = 1:nflies
            pp = Sc.allScores.postprocessed{f};
            if mean(pp) > ceiling_fly_thresh
                n_ceiling_flies = n_ceiling_flies + 1;
            end
        end
        frac = n_ceiling_flies / nflies;
        frac_ceiling_per_exp(i) = frac;
        if frac >= ceiling_exp_thresh
            exclude_exp(i) = true;
        end
        if mod(i, 50) == 0
            fprintf('  %d/%d experiments processed\n', i, nctrl_exp);
        end
    end
    n_excluded = sum(exclude_exp);
    fprintf('Ceiling filter: excluding %d / %d experiments\n', n_excluded, nctrl_exp);
    save(cache_file, 'exclude_exp', 'frac_ceiling_per_exp', 'n_excluded');
end

%% Build indices for 4 conditions
% Experiment level
exp_idx_raw = ctrl_exp_idx;                          % all control experiments
exp_idx_ceil = ctrl_exp_idx(~exclude_exp);           % ceiling-filtered

% Set level: for raw conditions, use setstats directly
set_idx_raw = ctrl_set_idx;

% For ceiling-filtered sets, we need to recompute set means from remaining experiments.
% Build mapping: for each control set, which control experiments belong to it?
set_to_exp = cell(1, numel(ctrl_set_idx));
set_to_exp_ceil = cell(1, numel(ctrl_set_idx));
for si = 1:numel(ctrl_set_idx)
    setid = ctrl_set_idx(si);
    % Find control experiments in this set
    in_set = S.setidx(ctrl_exp_idx) == setid;
    set_to_exp{si} = find(in_set);  % indices into ctrl_exp_idx
    set_to_exp_ceil{si} = find(in_set & ~exclude_exp(:));
end
% Which sets still have experiments after ceiling filtering?
set_has_exps = cellfun(@(x) ~isempty(x), set_to_exp_ceil);
fprintf('Sets with experiments after ceiling filter: %d / %d\n', sum(set_has_exps), numel(ctrl_set_idx));

%% Pre-generate bootstrap indices (shared across all stats)
rng('default');
boot_idx = randi(nctrl_set, nsamples, k_sets);          % for raw/normalized
nctrl_set_ceil = sum(set_has_exps);
boot_idx_ceil = randi(nctrl_set_ceil, nsamples, k_sets); % for ceiling-filtered

%% Initialize output tables
% Normality p-values: 4 conditions x 4 values (lillie_exp, dip_exp, lillie_set, dip_set) = 16 columns
pval_table = nan(nstats, 16);

% Bootstrap diagnostics: 4 conditions x 4 metrics (skewness, excess_kurtosis, n_sets, tail_ratio) = 16 columns
boot_diag_table = nan(nstats, 16);

%% Main loop over stats
fprintf('Plotting %d stats x 4 conditions...\n', nstats);
t0 = tic;

for si = 1:nstats
    fn = S.statfns{si};
    is_circ = S.iscircstat(si);

    % Extract experiment-level data for all 4 conditions
    exp_raw = S.expstats.means.(fn)(exp_idx_raw);
    exp_ceil = S.expstats.means.(fn)(exp_idx_ceil);
    exp_norm = S.expstats.normmeans.(fn)(exp_idx_raw);
    exp_ceil_norm = S.expstats.normmeans.(fn)(exp_idx_ceil);

    % Set-level: raw and normalized come directly from setstats
    set_raw = S.setstats.means.(fn)(set_idx_raw);
    set_norm = S.setstats.normmeans.(fn)(set_idx_raw);

    % Set-level ceiling-filtered: recompute from remaining experiments
    set_ceil = nan(1, numel(ctrl_set_idx));
    set_ceil_norm = nan(1, numel(ctrl_set_idx));
    for ki = 1:numel(ctrl_set_idx)
        exp_in_set = set_to_exp_ceil{ki};
        if isempty(exp_in_set), continue; end
        vals = S.expstats.means.(fn)(ctrl_exp_idx(exp_in_set));
        vals_norm = S.expstats.normmeans.(fn)(ctrl_exp_idx(exp_in_set));
        vals = vals(~isnan(vals));
        vals_norm = vals_norm(~isnan(vals_norm));
        if is_circ
            if ~isempty(vals), set_ceil(ki) = circ_mean(vals(:)); end
            if ~isempty(vals_norm), set_ceil_norm(ki) = circ_mean(vals_norm(:)); end
        else
            if ~isempty(vals), set_ceil(ki) = mean(vals); end
            if ~isempty(vals_norm), set_ceil_norm(ki) = mean(vals_norm); end
        end
    end
    set_ceil = set_ceil(set_has_exps);
    set_ceil_norm = set_ceil_norm(set_has_exps);

    % Bootstrap null distributions (mean of k_sets randomly sampled control sets)
    boot_raw = bootstrap_null(set_raw, boot_idx, is_circ);
    boot_ceil = bootstrap_null(set_ceil, boot_idx_ceil, is_circ);
    boot_norm = bootstrap_null(set_norm, boot_idx, is_circ);
    boot_ceil_norm = bootstrap_null(set_ceil_norm, boot_idx_ceil, is_circ);

    % Plot each condition and collect metrics
    cond_labels = {'raw', 'ceiling_filtered', 'normalized', 'ceiling_filtered_normalized'};
    exp_data = {exp_raw, exp_ceil, exp_norm, exp_ceil_norm};
    set_data = {set_raw, set_ceil, set_norm, set_ceil_norm};
    boot_data = {boot_raw, boot_ceil, boot_norm, boot_ceil_norm};

    for ci = 1:4
        [p_exp, p_set, bdiag, dip_exp, dip_set] = plot_one_stat(fn, exp_data{ci}, set_data{ci}, ...
            boot_data{ci}, is_circ, nbins, cond_labels{ci}, ...
            fullfile(outdir, subdirs{ci}), k_sets, dip_nboot);
        pval_table(si, (ci-1)*4 + (1:4)) = [p_exp, dip_exp, p_set, dip_set];
        boot_diag_table(si, (ci-1)*4 + (1:4)) = bdiag;
    end

    if mod(si, 50) == 0
        fprintf('  %d / %d stats done (%.1f min elapsed)\n', si, nstats, toc(t0)/60);
    end
end

%% Write normality p-value CSV (exp and set panels, lillie + dip)
csv_file = fullfile(outdir, 'normality_pvalues.csv');
fprintf('Writing normality p-values to %s\n', csv_file);
fid = fopen(csv_file, 'w');
fprintf(fid, 'statname,iscircular');
for ci = 1:4
    for metric = {'lillie_exp', 'dip_exp', 'lillie_set', 'dip_set'}
        fprintf(fid, ',%s_%s', subdirs{ci}, metric{1});
    end
end
fprintf(fid, '\n');
for si = 1:nstats
    fprintf(fid, '%s,%d', S.statfns{si}, S.iscircstat(si));
    for col = 1:16
        fprintf(fid, ',%.6e', pval_table(si, col));
    end
    fprintf(fid, '\n');
end
fclose(fid);

%% Write bootstrap diagnostics CSV
diag_file = fullfile(outdir, 'bootstrap_diagnostics.csv');
fprintf('Writing bootstrap diagnostics to %s\n', diag_file);
fid = fopen(diag_file, 'w');
fprintf(fid, 'statname,iscircular');
for ci = 1:4
    for metric = {'skewness', 'excess_kurtosis', 'n_ctrl_sets', 'tail_ratio'}
        fprintf(fid, ',%s_%s', subdirs{ci}, metric{1});
    end
end
fprintf(fid, '\n');
for si = 1:nstats
    fprintf(fid, '%s,%d', S.statfns{si}, S.iscircstat(si));
    for col = 1:16
        fprintf(fid, ',%.6e', boot_diag_table(si, col));
    end
    fprintf(fid, '\n');
end
fclose(fid);

fprintf('Done. Total time: %.1f min\n', toc(t0)/60);
fprintf('Output in: %s\n', outdir);


%% Helper functions

function boot_means = bootstrap_null(set_vals, boot_idx, is_circ)
    % Bootstrap null: for each sample, draw k sets and compute their mean.
    % boot_idx is nsamples x k, indexing into set_vals.
    set_vals = set_vals(:)';  % ensure row
    sampled = set_vals(boot_idx);  % nsamples x k
    % Remove rows where any sampled value is NaN
    bad = any(isnan(sampled), 2);
    if is_circ
        boot_means = circ_mean(sampled(~bad, :), [], 2);
    else
        boot_means = mean(sampled(~bad, :), 2);
    end
end

function [p_exp, p_set, bdiag, dip_exp, dip_set] = plot_one_stat(statname, exp_vals, set_vals, ...
        boot_vals, is_circ, nbins, cond_label, outdir, k_sets, dip_nboot)
    % Remove NaNs
    exp_vals = exp_vals(~isnan(exp_vals));
    n_ctrl_sets = numel(set_vals(~isnan(set_vals)));
    set_vals = set_vals(~isnan(set_vals));
    boot_vals = boot_vals(~isnan(boot_vals));

    if is_circ
        [p_exp, p_set, bdiag] = plot_circular_stat(statname, exp_vals, ...
            set_vals, boot_vals, nbins, cond_label, outdir, k_sets, n_ctrl_sets);
        dip_exp = NaN;  % dip test not applicable to circular data
        dip_set = NaN;
    else
        [p_exp, p_set, bdiag, dip_exp, dip_set] = plot_linear_stat(statname, exp_vals, ...
            set_vals, boot_vals, nbins, cond_label, outdir, k_sets, n_ctrl_sets, dip_nboot);
    end
end

function [p_exp, p_set, bdiag, dip_exp, dip_set] = plot_linear_stat(statname, exp_vals, ...
        set_vals, boot_vals, nbins, cond_label, outdir, k_sets, n_ctrl_sets, dip_nboot)
    fig = figure('Visible', 'off', 'Position', [100 100 1400 400]);

    % Compute shared bin edges across all 3 panels
    all_vals = [exp_vals(:); set_vals(:); boot_vals(:)];
    if isempty(all_vals) || all(isnan(all_vals))
        close(fig);
        p_exp = NaN; p_set = NaN; bdiag = [NaN NaN NaN NaN];
        dip_exp = NaN; dip_set = NaN;
        return;
    end
    all_vals = all_vals(~isnan(all_vals));
    lo = min(all_vals);
    hi = max(all_vals);
    if lo == hi
        lo = lo - 1;
        hi = hi + 1;
    end
    edges = linspace(lo, hi, nbins + 1);

    % Compute p-values for exp and set
    [p_exp, pstr_exp] = run_lillietest(exp_vals);
    [p_set, pstr_set] = run_lillietest(set_vals);

    % Compute Hartigan's dip test for exp and set
    [dip_exp, dipstr_exp] = run_diptest(exp_vals, dip_nboot);
    [dip_set, dipstr_set] = run_diptest(set_vals, dip_nboot);

    % Compute bootstrap diagnostics
    bdiag = compute_boot_diagnostics(boot_vals, n_ctrl_sets);
    boot_str = format_boot_diagnostics(bdiag);

    % Left panel: experiments
    subplot(1, 3, 1);
    histogram(exp_vals, edges, 'Normalization', 'pdf', ...
        'FaceColor', [0.3 0.5 0.8], 'EdgeColor', 'w');
    title(sprintf('%s [%s]\nExperiments (n=%d) %s %s', statname, cond_label, ...
        numel(exp_vals), pstr_exp, dipstr_exp), ...
        'Interpreter', 'none', 'FontSize', 8);
    xlabel('Value', 'Interpreter', 'none');
    ylabel('Probability density');

    % Middle panel: sets
    subplot(1, 3, 2);
    histogram(set_vals, edges, 'Normalization', 'pdf', ...
        'FaceColor', [0.8 0.4 0.3], 'EdgeColor', 'w');
    title(sprintf('%s [%s]\nSets (n=%d) %s %s', statname, cond_label, ...
        numel(set_vals), pstr_set, dipstr_set), ...
        'Interpreter', 'none', 'FontSize', 8);
    xlabel('Value', 'Interpreter', 'none');
    ylabel('Probability density');

    % Right panel: bootstrap null
    subplot(1, 3, 3);
    histogram(boot_vals, edges, 'Normalization', 'pdf', ...
        'FaceColor', [0.4 0.7 0.4], 'EdgeColor', 'w');
    title(sprintf('%s [%s]\nBootstrap null (n=%d, k=%d) %s', ...
        statname, cond_label, numel(boot_vals), k_sets, boot_str), ...
        'Interpreter', 'none', 'FontSize', 8);
    xlabel('Value', 'Interpreter', 'none');
    ylabel('Probability density');

    % Save
    savename = fullfile(outdir, [statname, '.png']);
    saveas(fig, savename);
    close(fig);
end

function [p_exp, p_set, bdiag] = plot_circular_stat(statname, exp_vals, ...
        set_vals, boot_vals, nbins, cond_label, outdir, k_sets, n_ctrl_sets)
    fig = figure('Visible', 'off', 'Position', [100 100 1400 450]);

    % Compute p-values (Rayleigh test) for exp and set
    p_exp = run_circ_rtest(exp_vals);
    p_set = run_circ_rtest(set_vals);

    % Bootstrap diagnostics (circular versions)
    bdiag = compute_boot_diagnostics_circ(boot_vals, n_ctrl_sets);
    boot_str = format_boot_diagnostics(bdiag);

    panels = {exp_vals, set_vals, boot_vals};
    colors = {[0.3 0.5 0.8], [0.8 0.4 0.3], [0.4 0.7 0.4]};
    panel_strs = { ...
        sprintf('Experiments (n=%d) Rayleigh p=%.2e', numel(exp_vals), p_exp), ...
        sprintf('Sets (n=%d) Rayleigh p=%.2e', numel(set_vals), p_set), ...
        sprintf('Bootstrap null (n=%d, k=%d) %s', numel(boot_vals), k_sets, boot_str)};

    for pi = 1:3
        ax = subplot(1, 3, pi);
        pos = get(ax, 'Position');
        delete(ax);
        pax = polaraxes('Position', pos);
        if ~isempty(panels{pi})
            polarhistogram(pax, panels{pi}, nbins, 'Normalization', 'pdf', ...
                'FaceColor', colors{pi}, 'EdgeColor', 'w');
        end
        title(pax, sprintf('%s [%s]\n%s', statname, cond_label, panel_strs{pi}), ...
            'Interpreter', 'none', 'FontSize', 8);
    end

    % Save
    savename = fullfile(outdir, [statname, '.png']);
    saveas(fig, savename);
    close(fig);
end

function bdiag = compute_boot_diagnostics(boot_vals, n_ctrl_sets)
    % Compute diagnostics for bootstrap null quality.
    % Returns [skewness, excess_kurtosis, n_ctrl_sets, tail_ratio]
    boot_vals = boot_vals(~isnan(boot_vals));
    if numel(boot_vals) < 4
        bdiag = [NaN, NaN, n_ctrl_sets, NaN];
        return;
    end
    sk = skewness(boot_vals);
    ek = kurtosis(boot_vals) - 3;  % excess kurtosis (normal = 0)
    % Tail ratio: observed fraction with |z| > 3 vs normal expectation (0.0027)
    z = (boot_vals - mean(boot_vals)) / std(boot_vals);
    frac_beyond_3 = mean(abs(z) > 3);
    tail_ratio = frac_beyond_3 / 0.0027;
    bdiag = [sk, ek, n_ctrl_sets, tail_ratio];
end

function bdiag = compute_boot_diagnostics_circ(boot_vals, n_ctrl_sets)
    % Circular version of bootstrap diagnostics.
    % Skewness and kurtosis use circular definitions.
    % Tail ratio uses angular distance from circular mean.
    boot_vals = boot_vals(~isnan(boot_vals));
    if numel(boot_vals) < 4
        bdiag = [NaN, NaN, n_ctrl_sets, NaN];
        return;
    end
    % Circular skewness and kurtosis (from circstat toolbox)
    mu = circ_mean(boot_vals(:));
    % Use angular deviations for skewness/kurtosis approximation
    diffs = circ_dist(boot_vals(:), mu);
    sk = skewness(diffs);
    ek = kurtosis(diffs) - 3;
    % Tail ratio: fraction of angular deviations beyond 3*circ_std
    cstd = circ_std(boot_vals(:));
    frac_beyond_3 = mean(abs(diffs) > 3 * cstd);
    tail_ratio = frac_beyond_3 / 0.0027;
    bdiag = [sk, ek, n_ctrl_sets, tail_ratio];
end

function str = format_boot_diagnostics(bdiag)
    % Format bootstrap diagnostics for title annotation
    str = sprintf('skew=%.2f kurt=%.2f n_{sets}=%d tail=%.1fx', ...
        bdiag(1), bdiag(2), bdiag(3), bdiag(4));
end

function [pval, pval_str] = run_lillietest(vals)
    vals = vals(~isnan(vals));
    if numel(vals) < 4
        pval = NaN;
        pval_str = 'lillie: n<4';
        return;
    end
    if all(vals == vals(1))
        pval = NaN;
        pval_str = 'lillie: constant';
        return;
    end
    try
        [~, pval] = lillietest(vals);
        pval_str = sprintf('lillie p=%.2e', pval);
    catch ME
        pval = NaN;
        pval_str = sprintf('lillie: %s', ME.message);
    end
end

function pval = run_circ_rtest(vals)
    vals = vals(~isnan(vals));
    if numel(vals) < 2
        pval = NaN;
        return;
    end
    try
        [pval, ~] = circ_rtest(vals(:));
    catch
        pval = NaN;
    end
end

function [pval, pval_str] = run_diptest(vals, nboot)
    % Run Hartigan's dip test for multimodality
    vals = vals(~isnan(vals));
    if numel(vals) < 4
        pval = NaN;
        pval_str = 'dip: n<4';
        return;
    end
    if all(vals == vals(1))
        pval = NaN;
        pval_str = 'dip: constant';
        return;
    end
    try
        [~, pval] = HartigansDipSignifTest(vals(:), nboot);
        pval_str = sprintf('dip p=%.2e', pval);
    catch ME
        pval = NaN;
        pval_str = sprintf('dip: %s', ME.message);
    end
end

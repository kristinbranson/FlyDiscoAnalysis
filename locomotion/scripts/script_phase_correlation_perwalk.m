%% Find metrics that separate high vs low phase walks (PER-WALK analysis)
% Three approaches:
%   1. Cohen's d (group comparison by day and median split)
%   2. Pearson correlation (phase as continuous variable)
%   3. AUC (classification discrimination)
%
% Note: Walks within the same fly are not independent - interpret with caution

modpath

%% Configuration
rootdir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20260211_exploringphaseissue';
outputdir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Locomotion_analysis/claudeplots_20260211';

% Experiment directories
expdirs = {
    % Dec 5 (low phase day)
    fullfile(rootdir, 'VNC2_YNA_K_162984_RigA_20231205T114519')
    fullfile(rootdir, 'VNC2_YNA_K_162984_RigB_20231205T114631')
    fullfile(rootdir, 'VNC2_YNA_K_162984_RigC_20231205T114709')
    fullfile(rootdir, 'VNC2_YNA_K_162984_RigD_20231205T114743')
    % Dec 6 (high phase day)
    fullfile(rootdir, 'VNC2_YNA_K_162984_RigA_20231206T125420')
    fullfile(rootdir, 'VNC2_YNA_K_162984_RigB_20231206T125502')
    fullfile(rootdir, 'VNC2_YNA_K_162984_RigC_20231206T125605')
    fullfile(rootdir, 'VNC2_YNA_K_162984_RigD_20231206T125700')
};

day_labels = [repmat({'Dec5'}, 4, 1); repmat({'Dec6'}, 4, 1)];

%% Load all per-walk data (or load from cache)
intermediate_file = fullfile(rootdir, 'all_walk_data_cache.mat');

if exist(intermediate_file, 'file')
    fprintf('Loading cached per-walk data from: %s\n', intermediate_file);
    load(intermediate_file, 'all_walk_data', 'walk_counter');
    fprintf('Loaded %d walks from cache.\n', walk_counter);
else
    fprintf('Loading per-walk data from all experiments...\n');

    all_walk_data = struct();
    walk_counter = 0;
    global_fly_counter = 0;

    for e = 1:numel(expdirs)
        expdir = expdirs{e};
        [~, expname] = fileparts(expdir);
        day = day_labels{e};

        datafile = fullfile(expdir, 'locomotionmetricsswingstanceboutstats.mat');
        if ~exist(datafile, 'file')
            fprintf('  WARNING: %s not found, skipping\n', datafile);
            continue;
        end

        data = load(datafile);
        fprintf('Processing %s (%s)...\n', expname, day);

        % Focus on LED off condition
        if isempty(data.walk_metrics_OFF) || ~isfield(data.walk_metrics_OFF, 'perwalk')
            fprintf('  No perwalk data, skipping\n');
            continue;
        end

        perwalk = data.walk_metrics_OFF.perwalk;
        nwalks_exp = numel(perwalk);
        nflies_exp = numel(data.walk_metrics_OFF.perfly);

        fprintf('  %d walks from %d flies\n', nwalks_exp, nflies_exp);

        % Metrics to extract from perwalk
        % Source: walk_metrics_OFF.perwalk(w).<metric>.mean
        %   These are means across frames within each walk bout
        walk_metrics_to_extract = {'velmag_ctr', 'forward_vel', 'backward_vel', ...
            'left_vel', 'right_vel', 'absdtheta', 'absdv_ctr', 'absdu_ctr', ...
            'CoM_stability', 'left_dtheta', 'right_dtheta'};

        for w = 1:nwalks_exp
            walk_counter = walk_counter + 1;

            % Metadata
            all_walk_data(walk_counter).exp_idx = e;
            all_walk_data(walk_counter).expname = expname;
            all_walk_data(walk_counter).day = day;
            all_walk_data(walk_counter).fly_in_exp = perwalk(w).fly;
            all_walk_data(walk_counter).global_fly_id = global_fly_counter + perwalk(w).fly;
            all_walk_data(walk_counter).walk_t0 = perwalk(w).walk_t0;
            all_walk_data(walk_counter).walk_t1 = perwalk(w).walk_t1;
            all_walk_data(walk_counter).walk_duration = perwalk(w).walk_t1 - perwalk(w).walk_t0;

            % Extract phase difference (target variable)
            % Source: perwalk(w).phasediff_hilbert.ipsi_post_2.mean
            if isfield(perwalk(w), 'phasediff_hilbert') && ...
               isfield(perwalk(w).phasediff_hilbert, 'ipsi_post_2') && ...
               isfield(perwalk(w).phasediff_hilbert.ipsi_post_2, 'mean')
                all_walk_data(walk_counter).phase_ipsi_post_2 = perwalk(w).phasediff_hilbert.ipsi_post_2.mean;
            else
                all_walk_data(walk_counter).phase_ipsi_post_2 = NaN;
            end

            % Extract walk metrics
            for wm = 1:numel(walk_metrics_to_extract)
                wmname = walk_metrics_to_extract{wm};
                if isfield(perwalk(w), wmname) && isfield(perwalk(w).(wmname), 'mean')
                    all_walk_data(walk_counter).(wmname) = perwalk(w).(wmname).mean;
                else
                    all_walk_data(walk_counter).(wmname) = NaN;
                end
            end
        end

        global_fly_counter = global_fly_counter + nflies_exp;
    end

    fprintf('\nTotal walks loaded: %d\n', walk_counter);
    fprintf('Total unique flies: %d\n', global_fly_counter);

    % Save intermediate file
    fprintf('Saving cached per-walk data to: %s\n', intermediate_file);
    save(intermediate_file, 'all_walk_data', 'walk_counter');
end

%% Build data table
fprintf('\nBuilding data table...\n');

% Get all numeric field names (excluding metadata)
all_fields = fieldnames(all_walk_data);
metadata_fields = {'exp_idx', 'expname', 'day', 'fly_in_exp', 'global_fly_id', 'walk_t0', 'walk_t1'};
numeric_fields = setdiff(all_fields, metadata_fields);

% List all metrics being analyzed
fprintf('\nMetrics being analyzed (%d total):\n', numel(numeric_fields));
for i = 1:numel(numeric_fields)
    fprintf('  %d. %s\n', i, numeric_fields{i});
end
fprintf('\n');

% Create arrays
phase_values = [all_walk_data.phase_ipsi_post_2]';
day_values = {all_walk_data.day}';
is_dec5 = strcmp(day_values, 'Dec5');
is_dec6 = strcmp(day_values, 'Dec6');

% Median split
valid_phase = ~isnan(phase_values);
median_phase = median(phase_values(valid_phase));
is_low_phase = phase_values < median_phase;
is_high_phase = phase_values >= median_phase;

fprintf('Median phase: %.4f rad\n', median_phase);
fprintf('Dec5 walks: %d, Dec6 walks: %d\n', sum(is_dec5), sum(is_dec6));
fprintf('Low phase walks: %d, High phase walks: %d\n', sum(is_low_phase & valid_phase), sum(is_high_phase & valid_phase));

%% ========================================
%% ANALYSIS 1: Cohen's d (Group Comparison)
%% ========================================
fprintf('\n========================================\n');
fprintf('ANALYSIS 1: Cohen''s d (Group Comparison)\n');
fprintf('========================================\n');

% --- 1a: By Day ---
fprintf('\n--- 1a: By Day (Dec5=low vs Dec6=high) ---\n');

results_cohend_day = [];
metric_counter = 0;

for m = 1:numel(numeric_fields)
    fname = numeric_fields{m};
    if strcmp(fname, 'phase_ipsi_post_2'), continue; end

    values = [all_walk_data.(fname)]';
    vals_dec5 = values(is_dec5); vals_dec5 = vals_dec5(~isnan(vals_dec5));
    vals_dec6 = values(is_dec6); vals_dec6 = vals_dec6(~isnan(vals_dec6));

    if numel(vals_dec5) < 3 || numel(vals_dec6) < 3, continue; end

    metric_counter = metric_counter + 1;
    n1 = numel(vals_dec5); n2 = numel(vals_dec6);
    mean1 = mean(vals_dec5); mean2 = mean(vals_dec6);
    std1 = std(vals_dec5); std2 = std(vals_dec6);
    pooled_std = sqrt(((n1-1)*std1^2 + (n2-1)*std2^2) / (n1+n2-2));
    cohens_d = (mean2 - mean1) / max(pooled_std, eps);
    [~, pval] = ttest2(vals_dec5, vals_dec6);

    results_cohend_day(metric_counter).metric = fname;
    results_cohend_day(metric_counter).mean_low = mean1;
    results_cohend_day(metric_counter).mean_high = mean2;
    results_cohend_day(metric_counter).cohens_d = cohens_d;
    results_cohend_day(metric_counter).pval = pval;
    results_cohend_day(metric_counter).n_low = n1;
    results_cohend_day(metric_counter).n_high = n2;
end

[~, idx] = sort(abs([results_cohend_day.cohens_d]), 'descend');
results_cohend_day = results_cohend_day(idx);

fprintf('\nTop 15 metrics by |Cohen''s d| (by day):\n');
fprintf('%-35s %10s %10s %10s %10s %8s %8s\n', 'Metric', 'Mean_Dec5', 'Mean_Dec6', 'Cohen_d', 'p-value', 'n_Dec5', 'n_Dec6');
fprintf('%s\n', repmat('-', 1, 95));
for i = 1:min(15, numel(results_cohend_day))
    r = results_cohend_day(i);
    fprintf('%-35s %10.4f %10.4f %10.3f %10.2e %8d %8d\n', r.metric, r.mean_low, r.mean_high, r.cohens_d, r.pval, r.n_low, r.n_high);
end

% --- 1b: By Median Split ---
fprintf('\n--- 1b: By Median Split on Phase ---\n');

results_cohend_median = [];
metric_counter = 0;

for m = 1:numel(numeric_fields)
    fname = numeric_fields{m};
    if strcmp(fname, 'phase_ipsi_post_2'), continue; end

    values = [all_walk_data.(fname)]';
    vals_low = values(is_low_phase & ~isnan(values));
    vals_high = values(is_high_phase & ~isnan(values));

    if numel(vals_low) < 3 || numel(vals_high) < 3, continue; end

    metric_counter = metric_counter + 1;
    n1 = numel(vals_low); n2 = numel(vals_high);
    mean1 = mean(vals_low); mean2 = mean(vals_high);
    std1 = std(vals_low); std2 = std(vals_high);
    pooled_std = sqrt(((n1-1)*std1^2 + (n2-1)*std2^2) / (n1+n2-2));
    cohens_d = (mean2 - mean1) / max(pooled_std, eps);
    [~, pval] = ttest2(vals_low, vals_high);

    results_cohend_median(metric_counter).metric = fname;
    results_cohend_median(metric_counter).mean_low = mean1;
    results_cohend_median(metric_counter).mean_high = mean2;
    results_cohend_median(metric_counter).cohens_d = cohens_d;
    results_cohend_median(metric_counter).pval = pval;
    results_cohend_median(metric_counter).n_low = n1;
    results_cohend_median(metric_counter).n_high = n2;
end

[~, idx] = sort(abs([results_cohend_median.cohens_d]), 'descend');
results_cohend_median = results_cohend_median(idx);

fprintf('\nTop 15 metrics by |Cohen''s d| (median split):\n');
fprintf('%-35s %10s %10s %10s %10s %8s %8s\n', 'Metric', 'Mean_Low', 'Mean_High', 'Cohen_d', 'p-value', 'n_Low', 'n_High');
fprintf('%s\n', repmat('-', 1, 95));
for i = 1:min(15, numel(results_cohend_median))
    r = results_cohend_median(i);
    fprintf('%-35s %10.4f %10.4f %10.3f %10.2e %8d %8d\n', r.metric, r.mean_low, r.mean_high, r.cohens_d, r.pval, r.n_low, r.n_high);
end

%% ========================================
%% ANALYSIS 2: Pearson Correlation
%% ========================================
fprintf('\n========================================\n');
fprintf('ANALYSIS 2: Pearson Correlation with Phase\n');
fprintf('========================================\n');

results_corr = [];
metric_counter = 0;

for m = 1:numel(numeric_fields)
    fname = numeric_fields{m};
    if strcmp(fname, 'phase_ipsi_post_2'), continue; end

    values = [all_walk_data.(fname)]';
    valid = ~isnan(phase_values) & ~isnan(values);

    if sum(valid) < 5, continue; end

    metric_counter = metric_counter + 1;
    [r, pval] = corr(phase_values(valid), values(valid));

    results_corr(metric_counter).metric = fname;
    results_corr(metric_counter).r = r;
    results_corr(metric_counter).pval = pval;
    results_corr(metric_counter).n = sum(valid);
end

[~, idx] = sort(abs([results_corr.r]), 'descend');
results_corr = results_corr(idx);

fprintf('\nTop 15 metrics by |Pearson r|:\n');
fprintf('%-35s %10s %12s %10s\n', 'Metric', 'r', 'p-value', 'n');
fprintf('%s\n', repmat('-', 1, 70));
for i = 1:min(15, numel(results_corr))
    r = results_corr(i);
    fprintf('%-35s %10.4f %12.2e %10d\n', r.metric, r.r, r.pval, r.n);
end

%% ========================================
%% ANALYSIS 3: AUC (ROC Discrimination)
%% ========================================
fprintf('\n========================================\n');
fprintf('ANALYSIS 3: AUC (ROC Discrimination)\n');
fprintf('========================================\n');

% --- 3a: By Day ---
fprintf('\n--- 3a: By Day (Dec5=low vs Dec6=high) ---\n');

results_auc_day = [];
metric_counter = 0;

for m = 1:numel(numeric_fields)
    fname = numeric_fields{m};
    if strcmp(fname, 'phase_ipsi_post_2'), continue; end

    values = [all_walk_data.(fname)]';
    vals_dec5 = values(is_dec5); vals_dec5 = vals_dec5(~isnan(vals_dec5));
    vals_dec6 = values(is_dec6); vals_dec6 = vals_dec6(~isnan(vals_dec6));

    if numel(vals_dec5) < 3 || numel(vals_dec6) < 3, continue; end

    metric_counter = metric_counter + 1;
    n1 = numel(vals_dec5); n2 = numel(vals_dec6);

    % Compute AUC via Mann-Whitney U (vectorized for speed)
    % Count how often dec6 > dec5
    [~, pval, stats] = ranksum(vals_dec5, vals_dec6);
    % U = stats.ranksum - n1*(n1+1)/2 gives U for group 1
    % AUC = P(dec6 > dec5) = 1 - U1/(n1*n2)
    U1 = stats.ranksum - n1*(n1+1)/2;
    auc = 1 - U1/(n1*n2);

    results_auc_day(metric_counter).metric = fname;
    results_auc_day(metric_counter).auc = auc;
    results_auc_day(metric_counter).pval = pval;
    results_auc_day(metric_counter).n1 = n1;
    results_auc_day(metric_counter).n2 = n2;
end

% Sort by distance from 0.5
[~, idx] = sort(abs([results_auc_day.auc] - 0.5), 'descend');
results_auc_day = results_auc_day(idx);

fprintf('\nTop 15 metrics by |AUC - 0.5| (by day):\n');
fprintf('%-35s %10s %12s %10s\n', 'Metric', 'AUC', 'p-value', 'Direction');
fprintf('%s\n', repmat('-', 1, 70));
for i = 1:min(15, numel(results_auc_day))
    r = results_auc_day(i);
    if r.auc > 0.5
        direction = 'Dec6>Dec5';
    else
        direction = 'Dec5>Dec6';
    end
    fprintf('%-35s %10.4f %12.2e %10s\n', r.metric, r.auc, r.pval, direction);
end

% --- 3b: By Median Split ---
fprintf('\n--- 3b: By Median Split on Phase ---\n');

results_auc_median = [];
metric_counter = 0;

for m = 1:numel(numeric_fields)
    fname = numeric_fields{m};
    if strcmp(fname, 'phase_ipsi_post_2'), continue; end

    values = [all_walk_data.(fname)]';
    vals_low = values(is_low_phase & ~isnan(values));
    vals_high = values(is_high_phase & ~isnan(values));

    if numel(vals_low) < 3 || numel(vals_high) < 3, continue; end

    metric_counter = metric_counter + 1;
    n1 = numel(vals_low); n2 = numel(vals_high);

    [~, pval, stats] = ranksum(vals_low, vals_high);
    U1 = stats.ranksum - n1*(n1+1)/2;
    auc = 1 - U1/(n1*n2);

    results_auc_median(metric_counter).metric = fname;
    results_auc_median(metric_counter).auc = auc;
    results_auc_median(metric_counter).pval = pval;
end

[~, idx] = sort(abs([results_auc_median.auc] - 0.5), 'descend');
results_auc_median = results_auc_median(idx);

fprintf('\nTop 15 metrics by |AUC - 0.5| (median split):\n');
fprintf('%-35s %10s %12s %10s\n', 'Metric', 'AUC', 'p-value', 'Direction');
fprintf('%s\n', repmat('-', 1, 70));
for i = 1:min(15, numel(results_auc_median))
    r = results_auc_median(i);
    if r.auc > 0.5
        direction = 'High>Low';
    else
        direction = 'Low>High';
    end
    fprintf('%-35s %10.4f %12.2e %10s\n', r.metric, r.auc, r.pval, direction);
end

%% ========================================
%% SUMMARY: Consensus across analyses
%% ========================================
fprintf('\n========================================\n');
fprintf('SUMMARY: Consensus across analyses\n');
fprintf('========================================\n');

% Get top 10 from each analysis
top_cohend_day = {results_cohend_day(1:min(10, numel(results_cohend_day))).metric};
top_cohend_median = {results_cohend_median(1:min(10, numel(results_cohend_median))).metric};
top_corr = {results_corr(1:min(10, numel(results_corr))).metric};
top_auc_day = {results_auc_day(1:min(10, numel(results_auc_day))).metric};
top_auc_median = {results_auc_median(1:min(10, numel(results_auc_median))).metric};

% Count appearances
all_tops = [top_cohend_day, top_cohend_median, top_corr, top_auc_day, top_auc_median];
unique_metrics = unique(all_tops);
counts = zeros(size(unique_metrics));
for i = 1:numel(unique_metrics)
    counts(i) = sum(strcmp(all_tops, unique_metrics{i}));
end

[counts_sorted, idx] = sort(counts, 'descend');
unique_metrics = unique_metrics(idx);

fprintf('\nMetrics appearing in multiple top-10 lists:\n');
fprintf('%-35s %s\n', 'Metric', 'Count (out of 5 analyses)');
fprintf('%s\n', repmat('-', 1, 50));
for i = 1:numel(unique_metrics)
    if counts_sorted(i) >= 2
        fprintf('%-35s %d\n', unique_metrics{i}, counts_sorted(i));
    end
end

%% Save results
save(fullfile(outputdir, 'phase_correlation_perwalk_results.mat'), ...
    'results_cohend_day', 'results_cohend_median', ...
    'results_corr', ...
    'results_auc_day', 'results_auc_median', ...
    'all_walk_data', 'median_phase', 'numeric_fields');
fprintf('\nResults saved to: %s\n', fullfile(outputdir, 'phase_correlation_perwalk_results.mat'));

%% Scatter plots for top consensus metrics
fprintf('\nGenerating scatter plots for top consensus metrics...\n');

% Get metrics that appear in 3+ analyses
consensus_metrics = unique_metrics(counts_sorted >= 3);
if numel(consensus_metrics) < 3
    consensus_metrics = unique_metrics(1:min(6, numel(unique_metrics)));
end

n_plots = min(6, numel(consensus_metrics));
figure('Position', [100, 100, 1200, 800]);
tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:n_plots
    fname = consensus_metrics{i};
    nexttile;

    values = [all_walk_data.(fname)]';

    hold on;
    % Subsample for plotting if too many points
    dec5_idx = find(is_dec5);
    dec6_idx = find(is_dec6);
    if numel(dec5_idx) > 2000
        dec5_idx = dec5_idx(randperm(numel(dec5_idx), 2000));
    end
    if numel(dec6_idx) > 2000
        dec6_idx = dec6_idx(randperm(numel(dec6_idx), 2000));
    end

    scatter(phase_values(dec5_idx), values(dec5_idx), 10, [0.1, 0.3, 0.8], 'filled', 'MarkerFaceAlpha', 0.3);
    scatter(phase_values(dec6_idx), values(dec6_idx), 10, [0.8, 0.1, 0.1], 'filled', 'MarkerFaceAlpha', 0.3);
    xline(median_phase, 'k--', 'LineWidth', 1.5);

    % Zoom in on the data range (exclude outliers)
    phase_valid = phase_values(~isnan(phase_values));
    xlim([prctile(phase_valid, 1), prctile(phase_valid, 99)]);

    xlabel('Phase diff ipsi\_post\_2 (rad)', 'Interpreter', 'none');
    ylabel(strrep(fname, '_', '\_'), 'Interpreter', 'tex');

    % Add correlation
    valid = ~isnan(phase_values) & ~isnan(values);
    if sum(valid) > 3
        [rho, ~] = corr(phase_values(valid), values(valid));
        title(sprintf('r = %.3f', rho), 'Interpreter', 'none');
    end

    if i == 1
        legend({'Dec5', 'Dec6', 'Median'}, 'Location', 'best');
    end
end

sgtitle('Per-walk: Top metrics separating high vs low phase', 'Interpreter', 'none');
exportgraphics(gcf, fullfile(outputdir, 'phase_correlation_perwalk_scatter.png'), 'Resolution', 150);
fprintf('Scatter plot saved.\n');

fprintf('\nDone.\n');
fprintf('\nNOTE: Walks within flies are not independent. Consider:\n');
fprintf('  - Mixed-effects models for proper inference\n');
fprintf('  - Within-fly correlations to check if effects hold within individuals\n');

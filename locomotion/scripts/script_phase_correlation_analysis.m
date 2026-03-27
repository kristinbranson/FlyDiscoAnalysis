%% Find metrics that separate high vs low phase flies
% Three approaches:
%   1. Cohen's d (group comparison by day and median split)
%   2. Pearson correlation (phase as continuous variable)
%   3. AUC (classification discrimination)

modpath

%% Configuration
rootdir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20260211_exploringphaseissue';
outputdir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Locomotion_analysis/claudeplots_20260211';

% Experiment directories (correct paths)
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

%% Load all per-fly data (or load from cache)
intermediate_file = fullfile(rootdir, 'all_fly_data_cache.mat');

if exist(intermediate_file, 'file')
    fprintf('Loading cached per-fly data from: %s\n', intermediate_file);
    load(intermediate_file, 'all_fly_data', 'fly_counter');
    fprintf('Loaded %d flies from cache.\n', fly_counter);
else
    fprintf('Loading per-fly data from all experiments...\n');

    all_fly_data = struct();
    fly_counter = 0;

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
    if isempty(data.walk_metrics_OFF)
        fprintf('  No walk_metrics_OFF, skipping\n');
        continue;
    end

    perfly = data.walk_metrics_OFF.perfly;
    nflies = numel(perfly);

    for f = 1:nflies
        fly_counter = fly_counter + 1;

        % Metadata
        all_fly_data(fly_counter).exp_idx = e;
        all_fly_data(fly_counter).expname = expname;
        all_fly_data(fly_counter).day = day;
        all_fly_data(fly_counter).fly_in_exp = f;

        % Extract phase difference (target variable)
        if isfield(perfly(f), 'phasediff_hilbert') && ...
           isfield(perfly(f).phasediff_hilbert, 'ipsi_post_2') && ...
           isfield(perfly(f).phasediff_hilbert.ipsi_post_2, 'walk_mean_fly')
            all_fly_data(fly_counter).phase_ipsi_post_2 = perfly(f).phasediff_hilbert.ipsi_post_2.walk_mean_fly;
        else
            all_fly_data(fly_counter).phase_ipsi_post_2 = NaN;
        end

        % Extract other phase metrics
        phase_groups = {'tripods_4', 'ipsi_ant_2', 'conta_L2R_3', 'abscontra_L2R_3'};
        for pg = 1:numel(phase_groups)
            pgname = phase_groups{pg};
            fname = ['phase_' pgname];
            if isfield(perfly(f), 'phasediff_hilbert') && ...
               isfield(perfly(f).phasediff_hilbert, pgname) && ...
               isfield(perfly(f).phasediff_hilbert.(pgname), 'walk_mean_fly')
                all_fly_data(fly_counter).(fname) = perfly(f).phasediff_hilbert.(pgname).walk_mean_fly;
            else
                all_fly_data(fly_counter).(fname) = NaN;
            end
        end

        % Extract walk metrics (scalar per-fly values)
        % Source: walk_metrics_OFF.perfly(f).<metric>.mean
        % These are means across ALL walk frames for that fly (not swing/stance specific)
        %   velmag_ctr - mean speed during walk frames
        %   forward_vel, backward_vel, left_vel, right_vel - directional velocities
        %   absdtheta - absolute angular velocity
        %   absdu_ctr, absdv_ctr - absolute lateral/forward acceleration
        %   CoM_stability - center of mass stability
        %   left_dtheta, right_dtheta - turning rates
        walk_metrics_to_extract = {'velmag_ctr', 'forward_vel', 'backward_vel', ...
            'left_vel', 'right_vel', 'absdtheta', 'absdv_ctr', 'absdu_ctr', ...
            'CoM_stability', 'left_dtheta', 'right_dtheta'};

        for wm = 1:numel(walk_metrics_to_extract)
            wmname = walk_metrics_to_extract{wm};
            if isfield(perfly(f), wmname) && isfield(perfly(f).(wmname), 'frm_mean_fly')
                all_fly_data(fly_counter).(wmname) = perfly(f).(wmname).frm_mean_fly;
            else
                all_fly_data(fly_counter).(wmname) = NaN;
            end
        end

        % Extract number of walks
        % Source: walk_metrics_OFF.perfly(f).phasediff_hilbert.ipsi_post_2.walk_n_fly
        %   n_walks - number of walk bouts for this fly
        if isfield(perfly(f), 'phasediff_hilbert') && ...
           isfield(perfly(f).phasediff_hilbert, 'ipsi_post_2') && ...
           isfield(perfly(f).phasediff_hilbert.ipsi_post_2, 'walk_n_fly')
            all_fly_data(fly_counter).n_walks = perfly(f).phasediff_hilbert.ipsi_post_2.walk_n_fly;
        else
            all_fly_data(fly_counter).n_walks = NaN;
        end
    end

    % Also extract bout metrics if available
    % Source: bout_metrics_OFF.perfly(f).perlimb(limb).swing.durations_time and
    %         bout_metrics_OFF.perfly(f).perlimb(limb).stance.durations_time
    %   swing_dur_* - mean of all swing bout durations for that limb (per fly)
    %   stance_dur_* - mean of all stance bout durations for that limb (per fly)
    if ~isempty(data.bout_metrics_OFF) && isfield(data.bout_metrics_OFF, 'perfly')
        bout_perfly = data.bout_metrics_OFF.perfly;
        limb_names = {'RF', 'RM', 'RH', 'LH', 'LM', 'LF'};

        for f = 1:nflies
            fly_idx = fly_counter - nflies + f;

            if f > numel(bout_perfly) || ~isfield(bout_perfly(f), 'perlimb')
                % No bout data for this fly, set all to NaN
                for limb = 1:6
                    all_fly_data(fly_idx).(['swing_dur_' limb_names{limb}]) = NaN;
                    all_fly_data(fly_idx).(['stance_dur_' limb_names{limb}]) = NaN;
                end
                continue;
            end

            % Extract mean swing/stance durations per limb
            for limb = 1:6
                limb_name = limb_names{limb};

                % Swing duration
                swing_fname = ['swing_dur_' limb_name];
                if limb <= numel(bout_perfly(f).perlimb) && ...
                   isfield(bout_perfly(f).perlimb(limb), 'swing') && ...
                   isfield(bout_perfly(f).perlimb(limb).swing, 'durations_time')
                    durs = bout_perfly(f).perlimb(limb).swing.durations_time;
                    if ~isempty(durs)
                        all_fly_data(fly_idx).(swing_fname) = mean(durs);
                    else
                        all_fly_data(fly_idx).(swing_fname) = NaN;
                    end
                else
                    all_fly_data(fly_idx).(swing_fname) = NaN;
                end

                % Stance duration
                stance_fname = ['stance_dur_' limb_name];
                if limb <= numel(bout_perfly(f).perlimb) && ...
                   isfield(bout_perfly(f).perlimb(limb), 'stance') && ...
                   isfield(bout_perfly(f).perlimb(limb).stance, 'durations_time')
                    durs = bout_perfly(f).perlimb(limb).stance.durations_time;
                    if ~isempty(durs)
                        all_fly_data(fly_idx).(stance_fname) = mean(durs);
                    else
                        all_fly_data(fly_idx).(stance_fname) = NaN;
                    end
                else
                    all_fly_data(fly_idx).(stance_fname) = NaN;
                end
            end
        end
    end
    end

    fprintf('\nTotal flies loaded: %d\n', fly_counter);

    % Save intermediate file
    fprintf('Saving cached per-fly data to: %s\n', intermediate_file);
    save(intermediate_file, 'all_fly_data', 'fly_counter');
end

%% Build data table
fprintf('\nBuilding data table...\n');

% Get all numeric field names (excluding metadata and other phase metrics)
all_fields = fieldnames(all_fly_data);
metadata_fields = {'exp_idx', 'expname', 'day', 'fly_in_exp'};
numeric_fields = setdiff(all_fields, metadata_fields);
% Exclude other phase metrics (we already know they correlate with phase_ipsi_post_2)
phase_metrics_to_exclude = numeric_fields(startsWith(numeric_fields, 'phase_'));
numeric_fields = setdiff(numeric_fields, phase_metrics_to_exclude);
fprintf('Excluding %d phase metrics from analysis.\n', numel(phase_metrics_to_exclude));

% List all metrics being analyzed
fprintf('\nMetrics being analyzed (%d total):\n', numel(numeric_fields));
for i = 1:numel(numeric_fields)
    fprintf('  %d. %s\n', i, numeric_fields{i});
end
fprintf('\n');

% Create arrays for each field
phase_values = [all_fly_data.phase_ipsi_post_2]';
day_values = {all_fly_data.day}';
is_dec5 = strcmp(day_values, 'Dec5');
is_dec6 = strcmp(day_values, 'Dec6');

% Median split
valid_phase = ~isnan(phase_values);
median_phase = median(phase_values(valid_phase));
is_low_phase = phase_values < median_phase;
is_high_phase = phase_values >= median_phase;

fprintf('Median phase: %.4f rad\n', median_phase);
fprintf('Dec5 flies: %d, Dec6 flies: %d\n', sum(is_dec5), sum(is_dec6));
fprintf('Low phase flies: %d, High phase flies: %d\n', sum(is_low_phase & valid_phase), sum(is_high_phase & valid_phase));

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

    values = [all_fly_data.(fname)]';
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
end

[~, idx] = sort(abs([results_cohend_day.cohens_d]), 'descend');
results_cohend_day = results_cohend_day(idx);

fprintf('\nTop 15 metrics by |Cohen''s d| (by day):\n');
fprintf('%-35s %10s %10s %10s %10s\n', 'Metric', 'Mean_Dec5', 'Mean_Dec6', 'Cohen_d', 'p-value');
fprintf('%s\n', repmat('-', 1, 80));
for i = 1:min(15, numel(results_cohend_day))
    r = results_cohend_day(i);
    fprintf('%-35s %10.4f %10.4f %10.3f %10.4f\n', r.metric, r.mean_low, r.mean_high, r.cohens_d, r.pval);
end

% --- 1b: By Median Split ---
fprintf('\n--- 1b: By Median Split on Phase ---\n');

results_cohend_median = [];
metric_counter = 0;

for m = 1:numel(numeric_fields)
    fname = numeric_fields{m};
    if strcmp(fname, 'phase_ipsi_post_2'), continue; end

    values = [all_fly_data.(fname)]';
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
end

[~, idx] = sort(abs([results_cohend_median.cohens_d]), 'descend');
results_cohend_median = results_cohend_median(idx);

fprintf('\nTop 15 metrics by |Cohen''s d| (median split):\n');
fprintf('%-35s %10s %10s %10s %10s\n', 'Metric', 'Mean_Low', 'Mean_High', 'Cohen_d', 'p-value');
fprintf('%s\n', repmat('-', 1, 80));
for i = 1:min(15, numel(results_cohend_median))
    r = results_cohend_median(i);
    fprintf('%-35s %10.4f %10.4f %10.3f %10.4f\n', r.metric, r.mean_low, r.mean_high, r.cohens_d, r.pval);
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

    values = [all_fly_data.(fname)]';
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
fprintf('%-35s %10s %10s %10s\n', 'Metric', 'r', 'p-value', 'n');
fprintf('%s\n', repmat('-', 1, 70));
for i = 1:min(15, numel(results_corr))
    r = results_corr(i);
    fprintf('%-35s %10.4f %10.4f %10d\n', r.metric, r.r, r.pval, r.n);
end

%% ========================================
%% ANALYSIS 3: AUC (ROC Discrimination)
%% ========================================
fprintf('\n========================================\n');
fprintf('ANALYSIS 3: AUC (ROC Discrimination)\n');
fprintf('========================================\n');

% Helper function to compute AUC using Mann-Whitney U statistic
% AUC = U / (n1 * n2), where U is Mann-Whitney U statistic

% --- 3a: By Day ---
fprintf('\n--- 3a: By Day (Dec5=low vs Dec6=high) ---\n');

results_auc_day = [];
metric_counter = 0;

for m = 1:numel(numeric_fields)
    fname = numeric_fields{m};
    if strcmp(fname, 'phase_ipsi_post_2'), continue; end

    values = [all_fly_data.(fname)]';
    vals_dec5 = values(is_dec5); vals_dec5 = vals_dec5(~isnan(vals_dec5));
    vals_dec6 = values(is_dec6); vals_dec6 = vals_dec6(~isnan(vals_dec6));

    if numel(vals_dec5) < 3 || numel(vals_dec6) < 3, continue; end

    metric_counter = metric_counter + 1;
    n1 = numel(vals_dec5); n2 = numel(vals_dec6);

    % Compute AUC via Mann-Whitney U
    % For each pair, count how often dec6 > dec5
    U = 0;
    for i = 1:n1
        for j = 1:n2
            if vals_dec6(j) > vals_dec5(i)
                U = U + 1;
            elseif vals_dec6(j) == vals_dec5(i)
                U = U + 0.5;
            end
        end
    end
    auc = U / (n1 * n2);

    % Wilcoxon rank-sum test for p-value
    [pval, ~] = ranksum(vals_dec5, vals_dec6);

    results_auc_day(metric_counter).metric = fname;
    results_auc_day(metric_counter).auc = auc;
    results_auc_day(metric_counter).pval = pval;
    results_auc_day(metric_counter).n1 = n1;
    results_auc_day(metric_counter).n2 = n2;
end

% Sort by distance from 0.5 (how discriminative)
[~, idx] = sort(abs([results_auc_day.auc] - 0.5), 'descend');
results_auc_day = results_auc_day(idx);

fprintf('\nTop 15 metrics by |AUC - 0.5| (by day):\n');
fprintf('%-35s %10s %10s %10s\n', 'Metric', 'AUC', 'p-value', 'Direction');
fprintf('%s\n', repmat('-', 1, 70));
for i = 1:min(15, numel(results_auc_day))
    r = results_auc_day(i);
    if r.auc > 0.5
        direction = 'Dec6>Dec5';
    else
        direction = 'Dec5>Dec6';
    end
    fprintf('%-35s %10.4f %10.4f %10s\n', r.metric, r.auc, r.pval, direction);
end

% --- 3b: By Median Split ---
fprintf('\n--- 3b: By Median Split on Phase ---\n');

results_auc_median = [];
metric_counter = 0;

for m = 1:numel(numeric_fields)
    fname = numeric_fields{m};
    if strcmp(fname, 'phase_ipsi_post_2'), continue; end

    values = [all_fly_data.(fname)]';
    vals_low = values(is_low_phase & ~isnan(values));
    vals_high = values(is_high_phase & ~isnan(values));

    if numel(vals_low) < 3 || numel(vals_high) < 3, continue; end

    metric_counter = metric_counter + 1;
    n1 = numel(vals_low); n2 = numel(vals_high);

    % Compute AUC via Mann-Whitney U
    U = 0;
    for i = 1:n1
        for j = 1:n2
            if vals_high(j) > vals_low(i)
                U = U + 1;
            elseif vals_high(j) == vals_low(i)
                U = U + 0.5;
            end
        end
    end
    auc = U / (n1 * n2);

    [pval, ~] = ranksum(vals_low, vals_high);

    results_auc_median(metric_counter).metric = fname;
    results_auc_median(metric_counter).auc = auc;
    results_auc_median(metric_counter).pval = pval;
end

[~, idx] = sort(abs([results_auc_median.auc] - 0.5), 'descend');
results_auc_median = results_auc_median(idx);

fprintf('\nTop 15 metrics by |AUC - 0.5| (median split):\n');
fprintf('%-35s %10s %10s %10s\n', 'Metric', 'AUC', 'p-value', 'Direction');
fprintf('%s\n', repmat('-', 1, 70));
for i = 1:min(15, numel(results_auc_median))
    r = results_auc_median(i);
    if r.auc > 0.5
        direction = 'High>Low';
    else
        direction = 'Low>High';
    end
    fprintf('%-35s %10.4f %10.4f %10s\n', r.metric, r.auc, r.pval, direction);
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
save(fullfile(outputdir, 'phase_correlation_results.mat'), ...
    'results_cohend_day', 'results_cohend_median', ...
    'results_corr', ...
    'results_auc_day', 'results_auc_median', ...
    'all_fly_data', 'median_phase', 'numeric_fields');
fprintf('\nResults saved to: %s\n', fullfile(outputdir, 'phase_correlation_results.mat'));

%% Scatter plots for top consensus metrics
fprintf('\nGenerating scatter plots for top consensus metrics...\n');

% Get metrics that appear in 3+ analyses
consensus_metrics = unique_metrics(counts_sorted >= 3);
if isempty(consensus_metrics)
    consensus_metrics = unique_metrics(1:min(6, numel(unique_metrics)));
end

n_plots = min(6, numel(consensus_metrics));
figure('Position', [100, 100, 1200, 800]);
tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:n_plots
    fname = consensus_metrics{i};
    nexttile;

    values = [all_fly_data.(fname)]';

    hold on;
    scatter(phase_values(is_dec5), values(is_dec5), 40, [0.1, 0.3, 0.8], 'filled', 'MarkerFaceAlpha', 0.6);
    scatter(phase_values(is_dec6), values(is_dec6), 40, [0.8, 0.1, 0.1], 'filled', 'MarkerFaceAlpha', 0.6);
    xline(median_phase, 'k--', 'LineWidth', 1.5);

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

sgtitle('Top metrics separating high vs low phase flies', 'Interpreter', 'none');
exportgraphics(gcf, fullfile(outputdir, 'phase_correlation_scatter.png'), 'Resolution', 150);
fprintf('Scatter plot saved.\n');

fprintf('\nDone.\n');

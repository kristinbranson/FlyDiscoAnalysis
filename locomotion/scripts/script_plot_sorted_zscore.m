%% script_plot_sorted_zscore.m
% Plot sorted z-scores (normalized difference from control) colored by
% significance for representative metrics.
% z-score = (line_normmean - controlmean) / controlstd

modpath;

%% load data
S = load('/groups/branson/bransonlab/flydisco_linelevel_VNC/CollectedVNC23PerFrameStats20260402.mat', ...
    'linestats', 'line_names', 'nlines', 'statfns', 'controlmean', 'controlstd', ...
    'qvalue_bigger_adj', 'qvalue_smaller_adj', 'idxcontrol');

plotdir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Locomotion_analysis/plots_20260402_sorteddata';
if ~isfolder(plotdir), mkdir(plotdir); end

fdr_alpha = 0.1;

%% metrics to plot — LEDon
metrics_LEDon = {
    'velmag_ctr__walk__LEDon__all'                          'body speed (walk) LEDon'
    'durations_time__swing__LEDon__all'                     'swing duration LEDon'
    'durations_time__stance__LEDon__all'                    'stance duration LEDon'
    'amplitude_BL__step__LEDon__all'                        'step amplitude LEDon'
    'instataeous_frequency_steps__step__LEDon__all'         'step frequency LEDon'
    'AEPy__step__LEDon__all'                                'foot placement (AEPy) LEDon'
    'mean_tips_speed_bodyref__swing__LEDon__all'            'leg tip speed (swing) LEDon'
    'CoM_stability__walk__LEDon__all'                       'CoM stability LEDon'
    'phasediff_hilbert__walk__LEDon__tripods_4'             'phase diff (tripods) LEDon'
    'absphasediff_hilbert__walk__LEDon__absRM_LM'           'abs phase diff (RM-LM) LEDon'
    'absdtheta__walk__LEDon__all'                           'turning rate LEDon'
};

%% metrics to plot — LEDoff
metrics_LEDoff = {
    'velmag_ctr__walk__LEDoff__all'                          'body speed (walk) LEDoff'
    'durations_time__swing__LEDoff__all'                     'swing duration LEDoff'
    'durations_time__stance__LEDoff__all'                    'stance duration LEDoff'
    'amplitude_BL__step__LEDoff__all'                        'step amplitude LEDoff'
    'instataeous_frequency_steps__step__LEDoff__all'         'step frequency LEDoff'
    'AEPy__step__LEDoff__all'                                'foot placement (AEPy) LEDoff'
    'mean_tips_speed_bodyref__swing__LEDoff__all'            'leg tip speed (swing) LEDoff'
    'CoM_stability__walk__LEDoff__all'                       'CoM stability LEDoff'
    'phasediff_hilbert__walk__LEDoff__tripods_4'             'phase diff (tripods) LEDoff'
    'absphasediff_hilbert__walk__LEDoff__absRM_LM'           'abs phase diff (RM-LM) LEDoff'
    'absdtheta__walk__LEDoff__all'                           'turning rate LEDoff'
};

metrics = [metrics_LEDon; metrics_LEDoff];

%% plot
for mi = 1:size(metrics, 1)
    statfn = metrics{mi, 1};
    label = metrics{mi, 2};
    stati = find(strcmp(S.statfns, statfn));
    if isempty(stati)
        fprintf('WARNING: stat %s not found, skipping\n', statfn);
        continue;
    end

    vals = S.linestats.normmeans.(statfn);
    if S.iscircstat(stati)
        zscore_vals = circ_dist(vals, S.controlmean(stati)) / real(S.controlstd(stati));
    else
        zscore_vals = (vals - S.controlmean(stati)) / S.controlstd(stati);
    end

    q_bigger = S.qvalue_bigger_adj(:, stati);
    q_smaller = S.qvalue_smaller_adj(:, stati);


    % classify each line
    is_sig_bigger = q_bigger < fdr_alpha & ~isnan(q_bigger);
    is_sig_smaller = q_smaller < fdr_alpha & ~isnan(q_smaller);
    is_control = false(S.nlines, 1);
    is_control(S.idxcontrol) = true;
    is_nonsig = ~is_sig_bigger & ~is_sig_smaller & ~is_control & ~isnan(zscore_vals);

    % sort by z-score
    [zs_sorted, si] = sort(zscore_vals);
    is_sig_bigger_sorted = is_sig_bigger(si);
    is_sig_smaller_sorted = is_sig_smaller(si);
    is_control_sorted = is_control(si);
    is_nonsig_sorted = is_nonsig(si);

    n_bigger = sum(is_sig_bigger);
    n_smaller = sum(is_sig_smaller);

    figure('Position', [100 100 1000 500], 'Visible', 'off');
    hold on;

    % plot each category, collecting handles for legend
    h = gobjects(0);
    leg = {};

    idx = find(is_nonsig_sorted);
    if ~isempty(idx)
        h(end+1) = plot(idx, zs_sorted(idx), '.', 'Color', [0.7 0.7 0.7], 'MarkerSize', 8);
        leg{end+1} = 'not significant';
    end

    idx = find(is_sig_smaller_sorted);
    if ~isempty(idx)
        h(end+1) = plot(idx, zs_sorted(idx), 'b.', 'MarkerSize', 10);
        leg{end+1} = sprintf('sig smaller (%d, FDR<%.1f)', n_smaller, fdr_alpha);
    end

    idx = find(is_sig_bigger_sorted);
    if ~isempty(idx)
        h(end+1) = plot(idx, zs_sorted(idx), 'r.', 'MarkerSize', 10);
        leg{end+1} = sprintf('sig bigger (%d, FDR<%.1f)', n_bigger, fdr_alpha);
    end

    idx = find(is_control_sorted);
    if ~isempty(idx)
        h(end+1) = plot(idx, zs_sorted(idx), 'k.', 'MarkerSize', 12);
        leg{end+1} = 'control';
    end

    yline(0, 'k:', 'LineWidth', 0.5);
    % symmetric y-axis around 0
    ymax = max(abs(zs_sorted));
    ylim([-ymax ymax]);
    legend(h, leg, 'Location', 'southeast', 'Interpreter', 'none');
    xlabel('lines (sorted by z-score)', 'Interpreter', 'none');
    ylabel('z-score (normalized difference)', 'Interpreter', 'none');
    title(sprintf('%s', label), 'Interpreter', 'none');

    saveas(gcf, fullfile(plotdir, sprintf('zscore_%s.png', statfn)));
    saveas(gcf, fullfile(plotdir, sprintf('zscore_%s.pdf', statfn)));
    close(gcf);

    fprintf('Plotted %s: range [%.1f, %.1f]\n', label, min(zs_sorted), max(zs_sorted));
end

fprintf('\nPlots saved to %s\n', plotdir);

%% ========================================================================
%% Z-score plots with max curvature elbows and significance coloring
%% ========================================================================

plotdir_elbow = '/groups/branson/home/robiea/Projects_data/FlyDisco/Locomotion_analysis/plots_20260402_zscore_elbow';
if ~isfolder(plotdir_elbow), mkdir(plotdir_elbow); end

for mi = 1:size(metrics, 1)
    statfn = metrics{mi, 1};
    label = metrics{mi, 2};
    stati = find(strcmp(S.statfns, statfn));
    if isempty(stati), continue; end

    vals = S.linestats.normmeans.(statfn);
    if S.iscircstat(stati)
        zscore_vals = circ_dist(vals, S.controlmean(stati)) / real(S.controlstd(stati));
    else
        zscore_vals = (vals - S.controlmean(stati)) / S.controlstd(stati);
    end

    q_bigger = S.qvalue_bigger_adj(:, stati);
    q_smaller = S.qvalue_smaller_adj(:, stati);

    % classify significance
    is_sig_bigger = q_bigger < fdr_alpha & ~isnan(q_bigger);
    is_sig_smaller = q_smaller < fdr_alpha & ~isnan(q_smaller);
    is_control = false(S.nlines, 1);
    is_control(S.idxcontrol) = true;
    is_nonsig = ~is_sig_bigger & ~is_sig_smaller & ~is_control & ~isnan(zscore_vals);

    n_bigger = sum(is_sig_bigger);
    n_smaller = sum(is_sig_smaller);

    % sort
    [zs_sorted, si] = sort(zscore_vals);
    x = (1:numel(zs_sorted))';
    n = numel(zs_sorted);

    % max curvature elbow detection
    zs_smooth = movmean(zs_sorted, 11);
    d2 = diff(zs_smooth, 2);

    z_left = nan; z_right = nan;
    knee_left = nan; knee_right = nan;
    n_below_elbow = 0; n_above_elbow = 0;

    idx_neg = find(zs_sorted(2:end-1) < 0);
    if numel(idx_neg) > 2
        [~, kl] = max(d2(idx_neg));
        knee_left = idx_neg(kl) + 1;
        z_left = zs_sorted(knee_left);
        n_below_elbow = knee_left;
    end

    idx_pos = find(zs_sorted(2:end-1) > 0);
    if numel(idx_pos) > 2
        [~, kr] = max(d2(idx_pos));
        knee_right = idx_pos(kr) + 1;
        z_right = zs_sorted(knee_right);
        n_above_elbow = n - knee_right + 1;
    end

    % plot with significance coloring and elbow lines
    figure('Position', [100 100 1000 500], 'Visible', 'off');
    hold on;

    h = gobjects(0);
    leg = {};

    idx = find(is_nonsig(si));
    if ~isempty(idx)
        h(end+1) = plot(idx, zs_sorted(idx), '.', 'Color', [0.7 0.7 0.7], 'MarkerSize', 8);
        leg{end+1} = 'not significant';
    end

    idx = find(is_sig_smaller(si));
    if ~isempty(idx)
        h(end+1) = plot(idx, zs_sorted(idx), 'b.', 'MarkerSize', 10);
        leg{end+1} = sprintf('sig smaller (%d, FDR<%.1f)', n_smaller, fdr_alpha);
    end

    idx = find(is_sig_bigger(si));
    if ~isempty(idx)
        h(end+1) = plot(idx, zs_sorted(idx), 'r.', 'MarkerSize', 10);
        leg{end+1} = sprintf('sig bigger (%d, FDR<%.1f)', n_bigger, fdr_alpha);
    end

    idx = find(is_control(si));
    if ~isempty(idx)
        h(end+1) = plot(idx, zs_sorted(idx), 'k.', 'MarkerSize', 12);
        leg{end+1} = 'control';
    end

    yline(0, 'k:', 'LineWidth', 0.5);

    % elbow lines with annotations
    if ~isnan(z_left)
        yline(z_left, 'b--', 'LineWidth', 1.5);
        text(n*0.02, z_left - 0.08*max(abs(zs_sorted)), ...
            sprintf('z=%.2f (%d lines)', z_left, n_below_elbow), ...
            'Color', 'b', 'FontSize', 9, 'Interpreter', 'none');
    end
    if ~isnan(z_right)
        yline(z_right, 'r--', 'LineWidth', 1.5);
        text(n*0.02, z_right + 0.08*max(abs(zs_sorted)), ...
            sprintf('z=%.2f (%d lines)', z_right, n_above_elbow), ...
            'Color', 'r', 'FontSize', 9, 'Interpreter', 'none');
    end

    ymax = max(abs(zs_sorted)); ylim([-ymax ymax]);
    legend(h, leg, 'Location', 'southeast', 'Interpreter', 'none');
    xlabel('lines (sorted by z-score)', 'Interpreter', 'none');
    ylabel('z-score (normalized difference)', 'Interpreter', 'none');
    title(sprintf('%s', label), 'Interpreter', 'none');

    saveas(gcf, fullfile(plotdir_elbow, sprintf('zscore_elbow_%s.png', statfn)));
    saveas(gcf, fullfile(plotdir_elbow, sprintf('zscore_elbow_%s.pdf', statfn)));
    close(gcf);

    fprintf('Plotted %s: elbow left z=%.2f (%d), right z=%.2f (%d)\n', ...
        label, z_left, n_below_elbow, z_right, n_above_elbow);
end

fprintf('\nZ-score + elbow plots saved to %s\n', plotdir_elbow);

% %% ========================================================================
% %% Alternative elbow methods (commented out — kept for reference)
% %% ========================================================================
%
% %% Method: Kneedle (max distance from diagonal)
% % Normalize x,y to [0,1] on each side of zero, find max distance from diagonal.
% % Sensitive to search range on asymmetric data. See plots_20260402_elbows/.
%
% %% Method: Two-line fit (minimize total residual)
% % For each candidate split, fit two lines and sum residuals. Pick min error.
% % Results between kneedle and max curvature. Slow (brute force search).

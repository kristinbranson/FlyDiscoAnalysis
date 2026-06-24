%% script_plot_dual_criteria.m
% Plot sorted z-scores with dual-criteria significance coloring:
% - LEDon vs control (original stat significant)
% - LEDon vs LEDoff (delta stat significant)
% - Both (dual-criteria hit)
%
% Colors:
%   light blue/red  = significant for LEDon-vs-control only
%   dark blue/red   = dual-criteria hit (both LEDon-vs-control AND LEDon-vs-LEDoff)
%   gray            = not significant
%   black           = control

modpath;

%% load data
S = load('/groups/branson/bransonlab/flydisco_linelevel_VNC/CollectedVNC23PerFrameStats20260406_pvalues.mat', ...
    'linestats', 'line_names', 'nlines', 'statfns', 'controlmean', 'controlstd', ...
    'qvalue_bigger_adj', 'qvalue_smaller_adj', 'idxcontrol', 'isdeltastat', 'iscircstat');

plotdir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Locomotion_analysis/plots_20260407_dualcriteria';
if ~isfolder(plotdir), mkdir(plotdir); end

fdr_alpha = 0.1;

%% metrics to plot (LEDon stats — delta counterpart found automatically)
metrics = {
    'velmag_ctr__walk__LEDon__all'                          'body speed (walk)'
    'durations_time__swing__LEDon__all'                     'swing duration'
    'durations_time__stance__LEDon__all'                    'stance duration'
    'amplitude_BL__step__LEDon__all'                        'step amplitude'
    'instataeous_frequency_steps__step__LEDon__all'         'step frequency'
    'AEPy__step__LEDon__all'                                'foot placement (AEPy)'
    'mean_tips_speed_bodyref__swing__LEDon__all'            'leg tip speed (swing)'
    'CoM_stability__walk__LEDon__all'                       'CoM stability'
    'phasediff_hilbert__walk__LEDon__tripods_4'             'phase diff (tripods)'
    'absphasediff_hilbert__walk__LEDon__absRM_LM'           'abs phase diff (RM-LM)'
    'absdtheta__walk__LEDon__all'                           'turning rate'
};
%% AR hack for DIG 
plotdir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Locomotion_analysis/plots_20260410_dualcriteria';
if ~isfolder(plotdir), mkdir(plotdir); end
metrics = {'AEPy_middle__walk__LEDon__pair2', 'foot placement (AEPy)'};
%% plot
for mi = 1:size(metrics, 1)
    statfn = metrics{mi, 1};
    label = metrics{mi, 2};
    stati_on = find(strcmp(S.statfns, statfn));

    % find corresponding delta stat
    fn_delta = [statfn '_norm_LEDoff'];
    stati_delta = find(strcmp(S.statfns, fn_delta));
    if isempty(stati_on) || isempty(stati_delta)
        fprintf('WARNING: missing stat pair for %s, skipping\n', statfn);
        continue;
    end

    % z-scores (using LEDon stat values)
    vals = S.linestats.normmeans.(statfn);
    if S.iscircstat(stati_on)
        % circular z-score: angular distance from control mean / circular std
        zscore_vals = circ_dist(vals, S.controlmean(stati_on)) / real(S.controlstd(stati_on));
    else
        zscore_vals = (vals - S.controlmean(stati_on)) / S.controlstd(stati_on);
    end

    % significance for LEDon-vs-control
    sig_bigger_on = S.qvalue_bigger_adj(:, stati_on) < fdr_alpha & ~isnan(S.qvalue_bigger_adj(:, stati_on));
    sig_smaller_on = S.qvalue_smaller_adj(:, stati_on) < fdr_alpha & ~isnan(S.qvalue_smaller_adj(:, stati_on));

    % significance for LEDon-vs-LEDoff (delta)
    sig_bigger_delta = S.qvalue_bigger_adj(:, stati_delta) < fdr_alpha & ~isnan(S.qvalue_bigger_adj(:, stati_delta));
    sig_smaller_delta = S.qvalue_smaller_adj(:, stati_delta) < fdr_alpha & ~isnan(S.qvalue_smaller_adj(:, stati_delta));

    % classify lines
    is_control = false(S.nlines, 1);
    is_control(S.idxcontrol) = true;

    % LEDon-vs-LEDoff significant (either direction)
    is_sig_delta = (sig_smaller_delta | sig_bigger_delta);

    % classify by vs-control significance (same as before)
    is_nonsig = ~sig_bigger_on & ~sig_smaller_on & ~is_control & ~isnan(zscore_vals);

    % sort by z-score
    [zs_sorted, si] = sort(zscore_vals);
    n = numel(zs_sorted);

    n_sig_smaller = sum(sig_smaller_on);
    n_sig_bigger = sum(sig_bigger_on);
    n_sig_delta_total = sum(is_sig_delta & ~is_control);

    % plot
    figure('Position', [100 100 1100 500], 'Visible', 'off');
    hold on;

    h = gobjects(0);
    leg = {};

    % gray dots: not significant vs control
    idx = find(is_nonsig(si));
    if ~isempty(idx)
        h(end+1) = plot(idx, zs_sorted(idx), '.', 'Color', [0.7 0.7 0.7], 'MarkerSize', 8);
        leg{end+1} = 'not significant';
    end

    % blue dots: significant smaller vs control
    idx = find(sig_smaller_on(si));
    if ~isempty(idx)
        h(end+1) = plot(idx, zs_sorted(idx), 'b.', 'MarkerSize', 10);
        leg{end+1} = sprintf('sig smaller (%d)', n_sig_smaller);
    end

    % red dots: significant bigger vs control
    idx = find(sig_bigger_on(si));
    if ~isempty(idx)
        h(end+1) = plot(idx, zs_sorted(idx), 'r.', 'MarkerSize', 10);
        leg{end+1} = sprintf('sig bigger (%d)', n_sig_bigger);
    end

    % black dot: control
    idx = find(is_control(si));
    if ~isempty(idx)
        h(end+1) = plot(idx, zs_sorted(idx), 'k.', 'MarkerSize', 12);
        leg{end+1} = 'control';
    end

    % overlay + markers on all LEDon-vs-LEDoff significant lines
    idx = find(is_sig_delta(si) & ~is_control(si));
    if ~isempty(idx)
        h(end+1) = plot(idx, zs_sorted(idx), '+', 'Color', [0 0 0], 'MarkerSize', 5, 'LineWidth', 0.5);
        leg{end+1} = sprintf('LEDon-vs-LEDoff sig (%d)', n_sig_delta_total);
    end

    yline(0, 'k:', 'LineWidth', 0.5);

    % max curvature elbow detection
    zs_smooth = movmean(zs_sorted, 11);
    d2 = diff(zs_smooth, 2);

    z_left = nan; z_right = nan;
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

    ymax = max(abs(zs_sorted));
    ylim([-ymax ymax]);
    legend(h, leg, 'Location', 'southeast', 'Interpreter', 'none', 'FontSize', 8);
    xlabel('lines (sorted by z-score)', 'Interpreter', 'none');
    ylabel('z-score (normalized difference)', 'Interpreter', 'none');
    title(sprintf('%s', label), 'Interpreter', 'none');

    plotfilename = fullfile(plotdir, sprintf('dual_%s', statfn));
    exportgraphics(gcf, [plotfilename '.png']);
    exportgraphics(gcf, [plotfilename '.pdf']);
    close(gcf);

    fprintf('Plotted %s: sig smaller=%d, sig bigger=%d, LEDon-vs-LEDoff sig=%d, elbow L=%.2f (%d) R=%.2f (%d)\n', ...
        label, n_sig_smaller, n_sig_bigger, n_sig_delta_total, z_left, n_below_elbow, z_right, n_above_elbow);
end

fprintf('\nPlots saved to %s\n', plotdir);

%% ========================================================================
%% Alternative: blue = vs control only, red = vs LEDoff only, purple = both
%% ========================================================================

plotdir2 = '/groups/branson/home/robiea/Projects_data/FlyDisco/Locomotion_analysis/plots_20260407_dualcriteria_v2';
if ~isfolder(plotdir2), mkdir(plotdir2); end

for mi = 1:size(metrics, 1)
    statfn = metrics{mi, 1};
    label = metrics{mi, 2};
    stati_on = find(strcmp(S.statfns, statfn));
    fn_delta = [statfn '_norm_LEDoff'];
    stati_delta = find(strcmp(S.statfns, fn_delta));
    if isempty(stati_on) || isempty(stati_delta), continue; end

    % z-scores
    if S.iscircstat(stati_on)
        zscore_vals = circ_dist(S.linestats.normmeans.(statfn), S.controlmean(stati_on)) / real(S.controlstd(stati_on));
    else
        zscore_vals = (S.linestats.normmeans.(statfn) - S.controlmean(stati_on)) / S.controlstd(stati_on);
    end

    % significance
    sig_on = (S.qvalue_bigger_adj(:, stati_on) < fdr_alpha & ~isnan(S.qvalue_bigger_adj(:, stati_on))) | ...
             (S.qvalue_smaller_adj(:, stati_on) < fdr_alpha & ~isnan(S.qvalue_smaller_adj(:, stati_on)));
    sig_delta = (S.qvalue_bigger_adj(:, stati_delta) < fdr_alpha & ~isnan(S.qvalue_bigger_adj(:, stati_delta))) | ...
                (S.qvalue_smaller_adj(:, stati_delta) < fdr_alpha & ~isnan(S.qvalue_smaller_adj(:, stati_delta)));

    is_control = false(S.nlines, 1);
    is_control(S.idxcontrol) = true;

    is_both = sig_on & sig_delta & ~is_control;
    is_only_control = sig_on & ~sig_delta & ~is_control;
    is_only_delta = ~sig_on & sig_delta & ~is_control;
    is_nonsig = ~sig_on & ~sig_delta & ~is_control & ~isnan(zscore_vals);

    [zs_sorted, si] = sort(zscore_vals);
    n = numel(zs_sorted);

    figure('Position', [100 100 1100 500], 'Visible', 'off');
    hold on;

    h = gobjects(0); leg = {};

    idx = find(is_nonsig(si));
    if ~isempty(idx)
        h(end+1) = plot(idx, zs_sorted(idx), '.', 'Color', [0.8 0.8 0.8], 'MarkerSize', 8);
        leg{end+1} = 'not significant';
    end

    idx = find(is_only_control(si));
    if ~isempty(idx)
        h(end+1) = plot(idx, zs_sorted(idx), '.', 'Color', [0.2 0.4 1.0], 'MarkerSize', 10);
        leg{end+1} = sprintf('vs control only (%d)', sum(is_only_control));
    end

    idx = find(is_only_delta(si));
    if ~isempty(idx)
        h(end+1) = plot(idx, zs_sorted(idx), '.', 'Color', [1.0 0.2 0.2], 'MarkerSize', 10);
        leg{end+1} = sprintf('LEDon-vs-LEDoff only (%d)', sum(is_only_delta));
    end

    idx = find(is_both(si));
    if ~isempty(idx)
        h(end+1) = plot(idx, zs_sorted(idx), '.', 'Color', [0.6 0.0 0.8], 'MarkerSize', 12);
        leg{end+1} = sprintf('both (%d)', sum(is_both));
    end

    idx = find(is_control(si));
    if ~isempty(idx)
        h(end+1) = plot(idx, zs_sorted(idx), 'k.', 'MarkerSize', 12);
        leg{end+1} = 'control';
    end

    yline(0, 'k:', 'LineWidth', 0.5);

    % elbow detection
    zs_smooth = movmean(zs_sorted, 11);
    d2 = diff(zs_smooth, 2);
    z_left = nan; z_right = nan; n_below_elbow = 0; n_above_elbow = 0;
    idx_neg = find(zs_sorted(2:end-1) < 0);
    if numel(idx_neg) > 2
        [~, kl] = max(d2(idx_neg));
        knee_left = idx_neg(kl) + 1;
        z_left = zs_sorted(knee_left); n_below_elbow = knee_left;
    end
    idx_pos = find(zs_sorted(2:end-1) > 0);
    if numel(idx_pos) > 2
        [~, kr] = max(d2(idx_pos));
        knee_right = idx_pos(kr) + 1;
        z_right = zs_sorted(knee_right); n_above_elbow = n - knee_right + 1;
    end
    if ~isnan(z_left)
        yline(z_left, 'b--', 'LineWidth', 1.5);
        text(n*0.02, z_left - 0.08*max(abs(zs_sorted)), sprintf('z=%.2f (%d)', z_left, n_below_elbow), ...
            'Color', 'b', 'FontSize', 9, 'Interpreter', 'none');
    end
    if ~isnan(z_right)
        yline(z_right, 'r--', 'LineWidth', 1.5);
        text(n*0.02, z_right + 0.08*max(abs(zs_sorted)), sprintf('z=%.2f (%d)', z_right, n_above_elbow), ...
            'Color', 'r', 'FontSize', 9, 'Interpreter', 'none');
    end

    ymax = max(abs(zs_sorted)); ylim([-ymax ymax]);
    legend(h, leg, 'Location', 'southeast', 'Interpreter', 'none', 'FontSize', 8);
    xlabel('lines (sorted by z-score)', 'Interpreter', 'none');
    ylabel('z-score (normalized difference)', 'Interpreter', 'none');
    title(sprintf('%s', label), 'Interpreter', 'none');

    plotfilename = fullfile(plotdir2, sprintf('dual_%s', statfn));
    exportgraphics(gcf, [plotfilename '.png']);
    exportgraphics(gcf, [plotfilename '.pdf']);
    close(gcf);

    fprintf('Plotted %s: ctrl-only=%d, delta-only=%d, both=%d\n', ...
        label, sum(is_only_control), sum(is_only_delta), sum(is_both));
end

fprintf('\nPlots saved to %s\n', plotdir2);

%% ========================================================================
%% Cool/warm color scheme: cool = smaller, warm = bigger
%% lighter = single criterion, darker = dual criterion
%% ========================================================================

% smaller (cool): light blue -> blue green -> blue
c_smaller_ctrl = [86, 180, 233]/255;   % light blue: vs control only
c_smaller_led  = [0, 158, 115]/255;    % blue green: vs LEDoff only
c_smaller_both = [0, 114, 178]/255;    % blue: both

% bigger (warm): yellow -> orange -> red
c_bigger_ctrl  = [240, 228, 66]/255;   % yellow: vs control only
c_bigger_led   = [230, 159, 0]/255;    % orange: vs LEDoff only
c_bigger_both  = [213, 94, 0]/255;     % red: both

plotdir3 = '/groups/branson/home/robiea/Projects_data/FlyDisco/Locomotion_analysis/plots_20260407_dualcriteria_v3';
if ~isfolder(plotdir3), mkdir(plotdir3); end

for mi = 1:size(metrics, 1)
    statfn = metrics{mi, 1};
    label = metrics{mi, 2};
    stati_on = find(strcmp(S.statfns, statfn));
    fn_delta = [statfn '_norm_LEDoff'];
    stati_delta = find(strcmp(S.statfns, fn_delta));
    if isempty(stati_on) || isempty(stati_delta), continue; end

    % z-scores (column vector)
    vals = S.linestats.normmeans.(statfn);
    if S.iscircstat(stati_on)
        zscore_vals = circ_dist(vals, S.controlmean(stati_on)) / real(S.controlstd(stati_on));
    else
        zscore_vals = (vals - S.controlmean(stati_on)) / S.controlstd(stati_on);
    end
    % ensure column
    zscore_vals = zscore_vals(:);

    % significance (all column vectors)
    sig_smaller_on = S.qvalue_smaller_adj(:, stati_on) < fdr_alpha & ~isnan(S.qvalue_smaller_adj(:, stati_on));
    sig_bigger_on = S.qvalue_bigger_adj(:, stati_on) < fdr_alpha & ~isnan(S.qvalue_bigger_adj(:, stati_on));
    sig_on = sig_smaller_on | sig_bigger_on;
    sig_delta = (S.qvalue_bigger_adj(:, stati_delta) < fdr_alpha & ~isnan(S.qvalue_bigger_adj(:, stati_delta))) | ...
                (S.qvalue_smaller_adj(:, stati_delta) < fdr_alpha & ~isnan(S.qvalue_smaller_adj(:, stati_delta)));

    is_control = false(S.nlines, 1);
    is_control(S.idxcontrol) = true;

    % 6 categories
    is_smaller_both     = sig_smaller_on & sig_delta & ~is_control;
    is_smaller_ctrlonly = sig_smaller_on & ~sig_delta & ~is_control;
    is_bigger_both      = sig_bigger_on & sig_delta & ~is_control;
    is_bigger_ctrlonly  = sig_bigger_on & ~sig_delta & ~is_control;

    % delta-only: split by z-score sign
    is_delta_only = sig_delta & ~sig_on & ~is_control;
    is_delta_only_neg = is_delta_only & zscore_vals < 0;
    is_delta_only_pos = is_delta_only & zscore_vals >= 0;

    is_nonsig = ~sig_on & ~sig_delta & ~is_control & ~isnan(zscore_vals);

    % sort
    [zs_sorted, si] = sort(zscore_vals);
    n = numel(zs_sorted);

    figure('Position', [100 100 1100 500], 'Visible', 'off');
    hold on;
    h = gobjects(0); leg = {};

    % gray: not significant
    idx = find(is_nonsig(si));
    if ~isempty(idx)
        h(end+1) = plot(idx, zs_sorted(idx), '.', 'Color', [0.8 0.8 0.8], 'MarkerSize', 8);
        leg{end+1} = 'not significant';
    end

    % SMALLER (cool)
    idx = find(is_smaller_ctrlonly(si));
    if ~isempty(idx)
        h(end+1) = plot(idx, zs_sorted(idx), '.', 'Color', c_smaller_ctrl, 'MarkerSize', 10);
        leg{end+1} = sprintf('smaller vs ctrl (%d)', sum(is_smaller_ctrlonly));
    end
    idx = find(is_delta_only_neg(si));
    if ~isempty(idx)
        h(end+1) = plot(idx, zs_sorted(idx), '.', 'Color', c_smaller_led, 'MarkerSize', 10);
        leg{end+1} = sprintf('smaller vs LEDoff (%d)', sum(is_delta_only_neg));
    end
    idx = find(is_smaller_both(si));
    if ~isempty(idx)
        h(end+1) = plot(idx, zs_sorted(idx), '.', 'Color', c_smaller_both, 'MarkerSize', 12);
        leg{end+1} = sprintf('smaller both (%d)', sum(is_smaller_both));
    end

    % BIGGER (warm)
    idx = find(is_bigger_ctrlonly(si));
    if ~isempty(idx)
        h(end+1) = plot(idx, zs_sorted(idx), '.', 'Color', c_bigger_ctrl, 'MarkerSize', 10);
        leg{end+1} = sprintf('bigger vs ctrl (%d)', sum(is_bigger_ctrlonly));
    end
    idx = find(is_delta_only_pos(si));
    if ~isempty(idx)
        h(end+1) = plot(idx, zs_sorted(idx), '.', 'Color', c_bigger_led, 'MarkerSize', 10);
        leg{end+1} = sprintf('bigger vs LEDoff (%d)', sum(is_delta_only_pos));
    end
    idx = find(is_bigger_both(si));
    if ~isempty(idx)
        h(end+1) = plot(idx, zs_sorted(idx), '.', 'Color', c_bigger_both, 'MarkerSize', 12);
        leg{end+1} = sprintf('bigger both (%d)', sum(is_bigger_both));
    end

    % control
    idx = find(is_control(si));
    if ~isempty(idx)
        h(end+1) = plot(idx, zs_sorted(idx), 'k.', 'MarkerSize', 12);
        leg{end+1} = 'control';
    end

    yline(0, 'k:', 'LineWidth', 0.5);

    % elbow detection (max curvature)
    zs_smooth = movmean(zs_sorted, 21);
    d2 = diff(zs_smooth, 2);
    z_left = nan; z_right = nan; n_below_elbow = 0; n_above_elbow = 0;
    idx_neg = find(zs_sorted(2:end-1) < 0);
    if numel(idx_neg) > 2
        [~, kl] = max(d2(idx_neg)); knee_left = idx_neg(kl)+1;
        z_left = zs_sorted(knee_left); n_below_elbow = knee_left;
    end
    idx_pos = find(zs_sorted(2:end-1) > 0);
    if numel(idx_pos) > 2
        [~, kr] = max(d2(idx_pos)); knee_right = idx_pos(kr)+1;
        z_right = zs_sorted(knee_right); n_above_elbow = n - knee_right + 1;
    end
    if ~isnan(z_left)
        yline(z_left, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5);
        text(n*0.02, z_left-0.08*max(abs(zs_sorted)), sprintf('z=%.2f (%d)', z_left, n_below_elbow), ...
            'Color', [0.3 0.3 0.3], 'FontSize', 9, 'Interpreter', 'none');
    end
    if ~isnan(z_right)
        yline(z_right, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5);
        text(n*0.02, z_right+0.08*max(abs(zs_sorted)), sprintf('z=%.2f (%d)', z_right, n_above_elbow), ...
            'Color', [0.3 0.3 0.3], 'FontSize', 9, 'Interpreter', 'none');
    end

    % highlight: both criteria AND beyond elbow (larger dots)
    is_beyond_left = zscore_vals <= z_left;
    is_beyond_right = zscore_vals >= z_right;
    is_highlight_smaller = is_smaller_both & is_beyond_left;
    is_highlight_bigger = is_bigger_both & is_beyond_right;
    n_highlight = sum(is_highlight_smaller) + sum(is_highlight_bigger);

    idx = find(is_highlight_smaller(si));
    if ~isempty(idx)
        h(end+1) = plot(idx, zs_sorted(idx), 'o', 'Color', c_smaller_both, 'MarkerSize', 3,'MarkerFaceColor','k','LineWidth',1);
        leg{end+1} = sprintf('smaller both+elbow (%d)', sum(is_highlight_smaller));
    end
    idx = find(is_highlight_bigger(si));
    if ~isempty(idx)
        h(end+1) = plot(idx, zs_sorted(idx), 'o', 'Color', c_bigger_both, 'MarkerSize', 3,'MarkerFaceColor','k','LineWidth',1);
        leg{end+1} = sprintf('bigger both+elbow (%d)', sum(is_highlight_bigger));
    end

    ymax = max(abs(zs_sorted)); ylim([-ymax ymax]);
    legend(h, leg, 'Location', 'southeast', 'Interpreter', 'none', 'FontSize', 7);
    xlabel('lines (sorted by z-score)', 'Interpreter', 'none');
    ylabel('z-score (normalized difference)', 'Interpreter', 'none');
    title(sprintf('%s', label), 'Interpreter', 'none');

    plotfilename = fullfile(plotdir3, sprintf('dual_%s', statfn));
    exportgraphics(gcf, [plotfilename '.png']);
    exportgraphics(gcf, [plotfilename '.pdf']);
    close(gcf);

    fprintf('Plotted %s: smaller(ctrl=%d,led=%d,both=%d) bigger(ctrl=%d,led=%d,both=%d) highlight=%d\n', ...
        label, sum(is_smaller_ctrlonly), sum(is_delta_only_neg), sum(is_smaller_both), ...
        sum(is_bigger_ctrlonly), sum(is_delta_only_pos), sum(is_bigger_both), n_highlight);
end

fprintf('\nPlots saved to %s\n', plotdir3);

%% ========================================================================
%% Simple version: gray = not dual, dark blue = dual smaller, dark red = dual bigger
%% No elbow lines
%% ========================================================================

c_dual_smaller = [0, 114, 178]/255;
c_dual_bigger  = [213, 94, 0]/255;

plotdir4 = '/groups/branson/home/robiea/Projects_data/FlyDisco/Locomotion_analysis/plots_20260408_dualcriteria_simple';
if ~isfolder(plotdir4), mkdir(plotdir4); end

for mi = 1:size(metrics, 1)
    statfn = metrics{mi, 1}; label = metrics{mi, 2};
    stati_on = find(strcmp(S.statfns, statfn));
    fn_delta = [statfn '_norm_LEDoff'];
    stati_delta = find(strcmp(S.statfns, fn_delta));
    if isempty(stati_on) || isempty(stati_delta), continue; end

    vals = S.linestats.normmeans.(statfn);
    if S.iscircstat(stati_on)
        zscore_vals = circ_dist(vals, S.controlmean(stati_on)) / real(S.controlstd(stati_on));
    else
        zscore_vals = (vals - S.controlmean(stati_on)) / S.controlstd(stati_on);
    end
    zscore_vals = zscore_vals(:);

    sig_smaller_on = S.qvalue_smaller_adj(:, stati_on) < fdr_alpha & ~isnan(S.qvalue_smaller_adj(:, stati_on));
    sig_bigger_on = S.qvalue_bigger_adj(:, stati_on) < fdr_alpha & ~isnan(S.qvalue_bigger_adj(:, stati_on));
    sig_delta = (S.qvalue_bigger_adj(:, stati_delta) < fdr_alpha & ~isnan(S.qvalue_bigger_adj(:, stati_delta))) | ...
                (S.qvalue_smaller_adj(:, stati_delta) < fdr_alpha & ~isnan(S.qvalue_smaller_adj(:, stati_delta)));
    is_control = false(S.nlines, 1); is_control(S.idxcontrol) = true;

    is_dual_smaller = sig_smaller_on & sig_delta & ~is_control;
    is_dual_bigger = sig_bigger_on & sig_delta & ~is_control;
    is_other = ~is_dual_smaller & ~is_dual_bigger & ~is_control & ~isnan(zscore_vals);

    [zs_sorted, si] = sort(zscore_vals); n = numel(zs_sorted);

    figure('Position', [100 100 1100 500], 'Visible', 'off'); hold on;
    h = gobjects(0); leg = {};

    idx = find(is_other(si));
    if ~isempty(idx), h(end+1) = plot(idx, zs_sorted(idx), '.', 'Color', [0.75 0.75 0.75], 'MarkerSize', 8); leg{end+1} = sprintf('other (%d)', sum(is_other)); end

    idx = find(is_dual_smaller(si));
    if ~isempty(idx), h(end+1) = plot(idx, zs_sorted(idx), '.', 'Color', c_dual_smaller, 'MarkerSize', 12); leg{end+1} = sprintf('dual smaller (%d)', sum(is_dual_smaller)); end

    idx = find(is_dual_bigger(si));
    if ~isempty(idx), h(end+1) = plot(idx, zs_sorted(idx), '.', 'Color', c_dual_bigger, 'MarkerSize', 12); leg{end+1} = sprintf('dual bigger (%d)', sum(is_dual_bigger)); end

    idx = find(is_control(si));
    if ~isempty(idx), h(end+1) = plot(idx, zs_sorted(idx), 'k.', 'MarkerSize', 12); leg{end+1} = 'control'; end

    yline(0, 'k:', 'LineWidth', 0.5);
    ymax = max(abs(zs_sorted)); ylim([-ymax ymax]);
    legend(h, leg, 'Location', 'southeast', 'Interpreter', 'none', 'FontSize', 8);
    xlabel('lines (sorted)', 'Interpreter', 'none');
    ylabel('STDs of control', 'Interpreter', 'none');
    title(sprintf('%s', label), 'Interpreter', 'none');

    plotfilename = fullfile(plotdir4, sprintf('dual_%s', statfn));
    exportgraphics(gcf, [plotfilename '.png']);
    exportgraphics(gcf, [plotfilename '.pdf']);
    close(gcf);

    fprintf('Plotted %s: dual smaller=%d, dual bigger=%d\n', label, sum(is_dual_smaller), sum(is_dual_bigger));
end

fprintf('\nPlots saved to %s\n', plotdir4);

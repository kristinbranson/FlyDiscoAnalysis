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
    'qvalue_bigger_adj', 'qvalue_smaller_adj', 'idxcontrol', 'isdeltastat');

plotdir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Locomotion_analysis/plots_20260406_dualcriteria';
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

    % dual-criteria: significant in both LEDon-vs-control AND LEDon-vs-LEDoff
    is_dual_smaller = sig_smaller_on & (sig_smaller_delta | sig_bigger_delta);
    is_dual_bigger = sig_bigger_on & (sig_smaller_delta | sig_bigger_delta);

    % LEDon-vs-control only (significant vs control but NOT vs LEDoff)
    is_onlycontrol_smaller = sig_smaller_on & ~is_dual_smaller;
    is_onlycontrol_bigger = sig_bigger_on & ~is_dual_bigger;

    is_nonsig = ~sig_bigger_on & ~sig_smaller_on & ~is_control & ~isnan(zscore_vals);

    % sort by z-score
    [zs_sorted, si] = sort(zscore_vals);
    n = numel(zs_sorted);

    % counts
    n_dual_smaller = sum(is_dual_smaller);
    n_dual_bigger = sum(is_dual_bigger);
    n_onlycontrol_smaller = sum(is_onlycontrol_smaller);
    n_onlycontrol_bigger = sum(is_onlycontrol_bigger);

    % plot
    figure('Position', [100 100 1100 500], 'Visible', 'off');
    hold on;

    h = gobjects(0);
    leg = {};

    % gray: not significant
    idx = find(is_nonsig(si));
    if ~isempty(idx)
        h(end+1) = plot(idx, zs_sorted(idx), '.', 'Color', [0.8 0.8 0.8], 'MarkerSize', 8);
        leg{end+1} = 'not significant';
    end

    % light blue: LEDon-vs-control only (smaller)
    idx = find(is_onlycontrol_smaller(si));
    if ~isempty(idx)
        h(end+1) = plot(idx, zs_sorted(idx), '.', 'Color', [0.6 0.7 1.0], 'MarkerSize', 10);
        leg{end+1} = sprintf('vs control only smaller (%d)', n_onlycontrol_smaller);
    end

    % light red: LEDon-vs-control only (bigger)
    idx = find(is_onlycontrol_bigger(si));
    if ~isempty(idx)
        h(end+1) = plot(idx, zs_sorted(idx), '.', 'Color', [1.0 0.6 0.6], 'MarkerSize', 10);
        leg{end+1} = sprintf('vs control only bigger (%d)', n_onlycontrol_bigger);
    end

    % dark blue: dual-criteria (smaller)
    idx = find(is_dual_smaller(si));
    if ~isempty(idx)
        h(end+1) = plot(idx, zs_sorted(idx), '.', 'Color', [0 0 0.8], 'MarkerSize', 12);
        leg{end+1} = sprintf('dual-criteria smaller (%d)', n_dual_smaller);
    end

    % dark red: dual-criteria (bigger)
    idx = find(is_dual_bigger(si));
    if ~isempty(idx)
        h(end+1) = plot(idx, zs_sorted(idx), '.', 'Color', [0.8 0 0], 'MarkerSize', 12);
        leg{end+1} = sprintf('dual-criteria bigger (%d)', n_dual_bigger);
    end

    % black: control
    idx = find(is_control(si));
    if ~isempty(idx)
        h(end+1) = plot(idx, zs_sorted(idx), 'k.', 'MarkerSize', 12);
        leg{end+1} = 'control';
    end

    yline(0, 'k:', 'LineWidth', 0.5);
    ymax = max(abs(zs_sorted));
    ylim([-ymax ymax]);
    legend(h, leg, 'Location', 'southeast', 'Interpreter', 'none', 'FontSize', 8);
    xlabel('lines (sorted by z-score)', 'Interpreter', 'none');
    ylabel('z-score (normalized difference)', 'Interpreter', 'none');
    title(sprintf('%s', label), 'Interpreter', 'none');

    saveas(gcf, fullfile(plotdir, sprintf('dual_%s.png', statfn)));
    saveas(gcf, fullfile(plotdir, sprintf('dual_%s.pdf', statfn)));
    close(gcf);

    fprintf('Plotted %s: dual smaller=%d, dual bigger=%d, control-only smaller=%d, control-only bigger=%d\n', ...
        label, n_dual_smaller, n_dual_bigger, n_onlycontrol_smaller, n_onlycontrol_bigger);
end

fprintf('\nPlots saved to %s\n', plotdir);

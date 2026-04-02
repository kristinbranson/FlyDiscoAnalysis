%% script_plot_sorted_linemeans.m
% Plot sorted normalized line means colored by significance for
% representative metrics.

modpath;

%% load data
S = load('/groups/branson/bransonlab/flydisco_linelevel_VNC/CollectedVNC23PerFrameStats20260331.mat', ...
    'linestats', 'line_names', 'nlines', 'statfns', ...
    'qvalue_bigger_adj', 'qvalue_smaller_adj', 'idxcontrol');

plotdir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Locomotion_analysis/plots_20260401_sorteddata';
if ~isfolder(plotdir), mkdir(plotdir); end

fdr_alpha = 0.1;

%% metrics to plot
metrics = {
    'velmag_ctr__walk__LEDon__all'                          'body speed (walk)'
    'durations_time__swing__LEDon__all'                     'swing duration'
    'durations_time__stance__LEDon__all'                    'stance duration'
    'amplitude_BL__step__LEDon__all'                        'step amplitude'
    'instataeous_frequency_steps__step__LEDon__all'         'step frequency'
    'AEPy__step__LEDon__all'                                'foot placement (AEPy)'
    'mean_tips_speed_bodyref__swing__LEDon__all'            'leg tip speed (swing)'
    'CoM_stability__walk__LEDon__all'                       'CoM stability'
    'phasediff_hilbert__walk__LEDon__tripods_4'             'phase diff (tripods, signed)'
    'absphasediff_hilbert__walk__LEDon__absRM_LM'           'abs phase diff (RM-LM)'
    'absdtheta__walk__LEDon__all'                           'turning rate'
};

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
    q_bigger = S.qvalue_bigger_adj(:, stati);
    q_smaller = S.qvalue_smaller_adj(:, stati);

    % HACK: signed circular stats have swapped bigger/smaller in ComputePValueBySampling
    % due to circ_dist sign convention. Swap them here until the root cause is fixed.
    % Only affects signed phase (uses circ_dist), NOT absphase (uses linear stats).
    is_circstat = startsWith(statfn, 'phase');
    if is_circstat
        q_tmp = q_bigger;
        q_bigger = q_smaller;
        q_smaller = q_tmp;
    end

    % classify each line
    is_sig_bigger = q_bigger < fdr_alpha & ~isnan(q_bigger);
    is_sig_smaller = q_smaller < fdr_alpha & ~isnan(q_smaller);
    is_control = false(S.nlines, 1);
    is_control(S.idxcontrol) = true;
    is_nonsig = ~is_sig_bigger & ~is_sig_smaller & ~is_control & ~isnan(vals);

    % sort by value
    [vals_sorted, si] = sort(vals);
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
        h(end+1) = plot(idx, vals_sorted(idx), '.', 'Color', [0.7 0.7 0.7], 'MarkerSize', 8);
        leg{end+1} = 'not significant';
    end

    idx = find(is_sig_smaller_sorted);
    if ~isempty(idx)
        h(end+1) = plot(idx, vals_sorted(idx), 'b.', 'MarkerSize', 10);
        leg{end+1} = sprintf('sig smaller (%d, FDR<%.1f)', n_smaller, fdr_alpha);
    end

    idx = find(is_sig_bigger_sorted);
    if ~isempty(idx)
        h(end+1) = plot(idx, vals_sorted(idx), 'r.', 'MarkerSize', 10);
        leg{end+1} = sprintf('sig bigger (%d, FDR<%.1f)', n_bigger, fdr_alpha);
    end

    idx = find(is_control_sorted);
    if ~isempty(idx)
        h(end+1) = plot(idx, vals_sorted(idx), 'k.', 'MarkerSize', 12);
        leg{end+1} = 'control';
    end

    legend(h, leg, 'Location', 'southeast', 'Interpreter', 'none');

    xlabel('lines (sorted by value)', 'Interpreter', 'none');
    ylabel('normalized mean', 'Interpreter', 'none');
    title(sprintf('%s', label), 'Interpreter', 'none');

    saveas(gcf, fullfile(plotdir, sprintf('sorted_%s.png', statfn)));
    saveas(gcf, fullfile(plotdir, sprintf('sorted_%s.pdf', statfn)));
    close(gcf);

    fprintf('Plotted %s: %d smaller, %d bigger\n', label, n_smaller, n_bigger);
end

fprintf('\nPlots saved to %s\n', plotdir);

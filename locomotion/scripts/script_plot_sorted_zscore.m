%% script_plot_sorted_zscore.m
% Plot sorted z-scores (normalized difference from control) colored by
% significance for representative metrics.
% z-score = (line_normmean - controlmean) / controlstd

modpath;

%% load data
S = load('/groups/branson/bransonlab/flydisco_linelevel_VNC/CollectedVNC23PerFrameStats20260331.mat', ...
    'linestats', 'line_names', 'nlines', 'statfns', 'controlmean', 'controlstd', ...
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
    zscore_vals = (vals - S.controlmean(stati)) / S.controlstd(stati);

    q_bigger = S.qvalue_bigger_adj(:, stati);
    q_smaller = S.qvalue_smaller_adj(:, stati);

    % HACK: signed circular stats have swapped bigger/smaller in ComputePValueBySampling
    % due to circ_dist sign convention. Swap them here until root cause is fixed.
    % Only affects signed phase (uses circ_dist), NOT absphase (uses linear stats).
    if startsWith(statfn, 'phase')
        q_tmp = q_bigger;
        q_bigger = q_smaller;
        q_smaller = q_tmp;
    end

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

%% script_plot_pvalue_summary.m
% Plot sorted p-values and FDR-corrected q-values for representative metrics.
% Reproduces style from plots_20260114 with current onfloor-filtered data.

modpath;

%% load data
S = load('/groups/branson/bransonlab/flydisco_linelevel_VNC/CollectedVNC23PerFrameStats20260331.mat', ...
    'pvalue_smaller_expected', 'pvalue_bigger_expected', ...
    'qvalue_smaller_adj', 'qvalue_bigger_adj', ...
    'statfns', 'line_names', 'nlines');

plotdir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Locomotion_analysis/plots_20260401_pvalues';
if ~isfolder(plotdir), mkdir(plotdir); end

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

    for direction = {'bigger', 'smaller'}
        dir = direction{1};

        % HACK: signed circular stats have swapped bigger/smaller in ComputePValueBySampling
        % due to circ_dist sign convention. Swap them here until root cause is fixed.
        % Only affects signed phase (uses circ_dist), NOT absphase (uses linear stats).
        is_circstat = startsWith(statfn, 'phase');
        if is_circstat
            swap_dir = struct('bigger','smaller','smaller','bigger');
            use_dir = swap_dir.(dir);
        else
            use_dir = dir;
        end

        if strcmp(use_dir, 'bigger')
            pvals = S.pvalue_bigger_expected(:, stati);
            qvals = S.qvalue_bigger_adj(:, stati);
            dir_label = 'higher than control';
        else
            pvals = S.pvalue_smaller_expected(:, stati);
            qvals = S.qvalue_smaller_adj(:, stati);
            dir_label = 'lower than control';
        end

        % remove NaN (control line)
        valid = ~isnan(pvals);
        pv = pvals(valid);
        qv = qvals(valid);

        % sort by p-value
        [pv_sorted, si] = sort(pv);
        qv_sorted = qv(si);

        figure('Position', [100 100 800 500], 'Visible', 'off');
        semilogy(1:numel(pv_sorted), pv_sorted, 'k-', 'LineWidth', 1.5);
        hold on;
        semilogy(1:numel(qv_sorted), qv_sorted, 'r-', 'LineWidth', 1);
        yline(0.1, 'r--', 'LineWidth', 0.5);
        yline(0.05, 'b--', 'LineWidth', 0.5);

        xlabel('lines (sorted)', 'Interpreter', 'none');
        ylabel('p-values', 'Interpreter', 'none');
        title(sprintf('Lines with %s %s during LEDon', label, dir_label), 'Interpreter', 'none');
        legend({'raw p-values', 'FDR corrected', 'FDR=0.1', 'p=0.05'}, 'Location', 'east');
        ylim([1e-5 1]);

        % count significant
        n_sig_raw = sum(pv_sorted < 0.05);
        n_sig_fdr = sum(qv_sorted < 0.1);
        text(numel(pv_sorted)*0.5, 2e-5, sprintf('p<0.05: %d lines, FDR<0.1: %d lines', n_sig_raw, n_sig_fdr), ...
            'Interpreter', 'none', 'FontSize', 10);

        % save
        fname = sprintf('pvalue_%s_%s', dir, statfn);
        saveas(gcf, fullfile(plotdir, [fname '.png']));
        saveas(gcf, fullfile(plotdir, [fname '.pdf']));
        close(gcf);
    end

    fprintf('Plotted %s\n', statfn);
end

fprintf('\n%d plots saved to %s\n', 2*size(metrics,1)*2, plotdir);

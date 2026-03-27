%% script_plot_allVNC_onceiling_v2
% Plot per-experiment onceiling % for all VNC/VNC2/VNC3 experiments.
% Separate panels by screen type and control vs SS lines.
%
% Plot 1: over time, 3 rows (VNC/VNC2/VNC3) x 2 cols (YNA controls / SS lines)
% Plot 2: sorted by median onceiling, same layout. Shading = IQR.
%
% Excludes VGLUT experiments.

datafile = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/allVNC_onceiling_resnet_v2_summary.mat';
plotdir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Locomotion_analysis/claudeplots_20260319_onceiling_allVNC/';
if ~isfolder(plotdir), mkdir(plotdir); end

load(datafile, 'expnames', 'pct_oc', 'dates_dt', 'screen_type', ...
    'is_control', 'fly_oc_per_exp');

is_vglut = contains(expnames, 'VGLUT');
is_ss = ~is_control & ~is_vglut;

% Two datasets: VNC alone, VNC2+VNC3 together
datasets = {
    {'VNC'},        'VNC'
    {'VNC2','VNC3'}, 'VNC2+VNC3'
};

% Fixed limits
ylims = [0 100];

for di = 1:size(datasets, 1)
    st_set = datasets{di, 1};
    ds_label = datasets{di, 2};

    in_dataset = ismember(screen_type, st_set) & ~is_vglut;

    %% Over time: 2 rows (controls, SS lines) x 1 col
    figure('Position', [100 100 1200 600]);
    tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

    for col = 1:2
        nexttile;
        hold on;

        if col == 1
            idx = in_dataset & is_control;
            label = [ds_label ' controls (YNA)'];
            color = [0.8 0.2 0.2];
        else
            idx = in_dataset & is_ss;
            label = [ds_label ' SS lines'];
            color = [0.4 0.4 0.4];
        end

        scatter(dates_dt(idx), pct_oc(idx), 8, color, 'filled', 'MarkerFaceAlpha', 0.4);

        xlabel('Date');
        ylabel('% onceiling', 'Interpreter', 'none');
        title(sprintf('%s (n=%d)', label, sum(idx)), 'Interpreter', 'none');
        ax = gca;
        ax.XAxis.TickLabelFormat = 'MMM yyyy';
        xtickangle(45);
        ylim(ylims);
        grid on;
    end

    fname = lower(strrep(ds_label, '+', ''));
    exportgraphics(gcf, fullfile(plotdir, [fname '_onceiling_over_time.png']), 'Resolution', 150);
    exportgraphics(gcf, fullfile(plotdir, [fname '_onceiling_over_time.pdf']), 'ContentType', 'vector');
    fprintf('Time plot saved for %s.\n', ds_label);

    %% Sorted by median onceiling: 2 rows (controls, SS lines) x 1 col
    figure('Position', [100 100 1200 600]);
    tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

    for col = 1:2
        nexttile;
        hold on;

        if col == 1
            idx = find(in_dataset & is_control);
            label = [ds_label ' controls (YNA)'];
            color = [0.8 0.2 0.2];
        else
            idx = find(in_dataset & is_ss);
            label = [ds_label ' SS lines'];
            color = [0.4 0.4 0.4];
        end

        n_sub = numel(idx);
        if n_sub == 0, continue; end

        % Compute per-experiment quartiles
        med_sub = nan(n_sub, 1);
        q25_sub = nan(n_sub, 1);
        q75_sub = nan(n_sub, 1);
        for j = 1:n_sub
            fp = fly_oc_per_exp{idx(j)};
            med_sub(j) = median(fp);
            q25_sub(j) = prctile(fp, 25);
            q75_sub(j) = prctile(fp, 75);
        end

        [~, si] = sort(med_sub);
        med_sub = med_sub(si);
        q25_sub = q25_sub(si);
        q75_sub = q75_sub(si);

        % IQR shading
        x_fill = [1:n_sub, n_sub:-1:1];
        y_fill = [q25_sub(:)', fliplr(q75_sub(:)')];
        fill(x_fill, y_fill, [0.85 0.85 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.5);

        % Median dots
        scatter(1:n_sub, med_sub, 10, color, 'filled', 'MarkerFaceAlpha', 0.6);

        % Median line
        plot(1:n_sub, med_sub, 'k-', 'LineWidth', 0.5);

        xlabel('Experiments (sorted)');
        ylabel('% onceiling per fly', 'Interpreter', 'none');
        title(sprintf('%s (n=%d)', label, n_sub), 'Interpreter', 'none');
        legend({'IQR', 'median'}, 'Location', 'northwest');
        xlim([0 n_sub+1]);
        ylim(ylims);
        grid on;

        % Threshold line
        n_above = sum(med_sub > 50);
        if n_above > 0
            xline(n_sub - n_above + 0.5, 'r--', 'LineWidth', 1);
            text(n_sub - n_above, 90, sprintf('  >50%%: %d (%.0f%%)', ...
                n_above, 100*n_above/n_sub), 'Color', 'r', 'FontSize', 8);
        end
    end

    exportgraphics(gcf, fullfile(plotdir, [fname '_onceiling_sorted.png']), 'Resolution', 150);
    exportgraphics(gcf, fullfile(plotdir, [fname '_onceiling_sorted.pdf']), 'ContentType', 'image');
    fprintf('Sorted plot saved for %s.\n', ds_label);
end

%% Summary
fprintf('\nSummary (excluding VGLUT):\n');
for di = 1:size(datasets, 1)
    st_set = datasets{di, 1};
    ds_label = datasets{di, 2};
    in_ds = ismember(screen_type, st_set) & ~is_vglut;
    ctrl = in_ds & is_control;
    ss = in_ds & is_ss;
    fprintf('  %s controls: n=%d, median=%.1f%%, >50%%=%d (%.1f%%)\n', ...
        ds_label, sum(ctrl), median(pct_oc(ctrl)), sum(pct_oc(ctrl)>50), 100*mean(pct_oc(ctrl)>50));
    fprintf('  %s SS lines: n=%d, median=%.1f%%, >50%%=%d (%.1f%%)\n', ...
        ds_label, sum(ss), median(pct_oc(ss)), sum(pct_oc(ss)>50), 100*mean(pct_oc(ss)>50));
end

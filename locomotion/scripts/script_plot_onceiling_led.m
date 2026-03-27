%% script_plot_onceiling_led
% Compare onceiling % during LED on vs LED off periods.
% Each dot = one experiment. Diagonal = no LED effect.
%
% Data: 258 experiments across 5 days (2 high, 3 low ceiling).
% LED on = optogenetic inactivation active.
% pct_oc_led_on/off: % of fly-frames classified as onceiling (pooled across flies)
%   during LED on vs LED off periods.

datadir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20260317_oncelingHiLowDays/';
plotdir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Locomotion_analysis/claudeplots_20260317_onceiling_hilow/';
if ~isfolder(plotdir), mkdir(plotdir); end

load(fullfile(datadir, 'onceiling_led_summary.mat'), 'results');

grid_dates = {
    '20210406', ''
    '20220525', '20230706'
    '20240806', '20240724'
};
grid_labels = {
    '2021-04-06 (VNC, low)',   ''
    '2022-05-25 (VNC2, low)',  '2023-07-06 (VNC2, HIGH)'
    '2024-08-06 (VNC3, low)',  '2024-07-24 (VNC3, HIGH)'
};

% Fixed limits
axlims = [0 100];

%% Plot 1: LED on vs LED off scatter, colored by control/VGLUT/other
figure('Position', [100 100 1000 900]);
tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for row = 1:3
    for col = 1:2
        nexttile;
        hold on;
        d = grid_dates{row, col};
        if isempty(d)
            text(0.5, 0.5, 'No high-ceiling VNC day', ...
                'HorizontalAlignment', 'center', 'Units', 'normalized');
            xlim(axlims); ylim(axlims);
            grid on;
            continue;
        end

        idx = strcmp({results.date}, d);
        rd = results(idx);

        % Other lines
        other = ~[rd.is_control] & ~[rd.is_vglut];
        scatter([rd(other).pct_oc_led_off], [rd(other).pct_oc_led_on], 20, ...
            [0.6 0.6 0.6], 'filled', 'MarkerFaceAlpha', 0.5);

        % Controls
        ctrl = [rd.is_control];
        scatter([rd(ctrl).pct_oc_led_off], [rd(ctrl).pct_oc_led_on], 30, ...
            [0.8 0.2 0.2], 'filled', 'MarkerFaceAlpha', 0.7);

        % VGLUT
        vg = [rd.is_vglut];
        scatter([rd(vg).pct_oc_led_off], [rd(vg).pct_oc_led_on], 30, ...
            [0.2 0.6 0.2], 'filled', 'MarkerFaceAlpha', 0.7);

        plot(axlims, axlims, 'k--', 'LineWidth', 0.5);

        xlim(axlims); ylim(axlims);
        xlabel('% onceiling LED off');
        ylabel('% onceiling LED on');
        title(grid_labels{row, col}, 'Interpreter', 'none');
        legend({'other', 'control', 'VGLUT', 'identity'}, 'Location', 'best');
        axis square;
        grid on;
    end
end

exportgraphics(gcf, fullfile(plotdir, 'onceiling_led_on_vs_off.png'), 'Resolution', 150);
exportgraphics(gcf, fullfile(plotdir, 'onceiling_led_on_vs_off.pdf'), 'ContentType', 'vector');
fprintf('Plot 1 saved.\n');

%% Plot 2: difference (LED on - LED off) vs overall OC%
figure('Position', [100 100 1000 900]);
tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% Compute common y limits for difference plot
all_diff = [results.pct_oc_led_on] - [results.pct_oc_led_off];
diff_lim = max(abs(all_diff)) * 1.1;
diff_ylims = [-diff_lim diff_lim];

for row = 1:3
    for col = 1:2
        nexttile;
        hold on;
        d = grid_dates{row, col};
        if isempty(d)
            text(0.5, 0.5, 'No high-ceiling VNC day', ...
                'HorizontalAlignment', 'center', 'Units', 'normalized');
            xlim(axlims); ylim(diff_ylims);
            grid on;
            continue;
        end

        idx = strcmp({results.date}, d);
        rd = results(idx);
        diffs = [rd.pct_oc_led_on] - [rd.pct_oc_led_off];

        other = ~[rd.is_control] & ~[rd.is_vglut];
        scatter([rd(other).pct_oc_all], diffs(other), 20, ...
            [0.6 0.6 0.6], 'filled', 'MarkerFaceAlpha', 0.5);

        ctrl = [rd.is_control];
        scatter([rd(ctrl).pct_oc_all], diffs(ctrl), 30, ...
            [0.8 0.2 0.2], 'filled', 'MarkerFaceAlpha', 0.7);

        vg = [rd.is_vglut];
        scatter([rd(vg).pct_oc_all], diffs(vg), 30, ...
            [0.2 0.6 0.2], 'filled', 'MarkerFaceAlpha', 0.7);

        plot(axlims, [0 0], 'k--', 'LineWidth', 0.5);

        xlim(axlims); ylim(diff_ylims);
        xlabel('% onceiling (overall)');
        ylabel('LED on - LED off (pp)');
        title(grid_labels{row, col}, 'Interpreter', 'none');
        legend({'other', 'control', 'VGLUT', 'zero'}, 'Location', 'best');
        grid on;
    end
end

exportgraphics(gcf, fullfile(plotdir, 'onceiling_led_diff_vs_overall.png'), 'Resolution', 150);
exportgraphics(gcf, fullfile(plotdir, 'onceiling_led_diff_vs_overall.pdf'), 'ContentType', 'vector');
fprintf('Plot 2 saved.\n');

%% Print experiments with largest LED effect
fprintf('\nTop experiments by |LED on - LED off| difference:\n');
all_diffs = [results.pct_oc_led_on] - [results.pct_oc_led_off];
[~, si] = sort(abs(all_diffs), 'descend');
fprintf('%-55s %7s %7s %8s\n', 'Experiment', 'LEDoff', 'LEDon', 'Diff');
fprintf('%s\n', repmat('-', 1, 80));
for i = 1:15
    j = si(i);
    fprintf('%-55s %6.1f%% %6.1f%% %+7.1f%%\n', ...
        results(j).expname, results(j).pct_oc_led_off, results(j).pct_oc_led_on, all_diffs(j));
end

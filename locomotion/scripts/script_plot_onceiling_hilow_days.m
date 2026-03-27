%% script_plot_onceiling_hilow_days
% Per-fly onceiling % vs time of day for 5 selected days (2 high, 3 low ceiling).
% Each dot = one fly trajectory. Controls highlighted.
%
% Data: 258 experiments, 2542 flies from:
%   High ceiling: 2023-07-06 (VNC2, 95% ctrl median), 2024-07-24 (VNC3, 91%)
%   Low ceiling:  2021-04-06 (VNC, 0.6%), 2022-05-25 (VNC2, 0.2%), 2024-08-06 (VNC3, 0.9%)
%
% pct_oc per fly: % of that fly's frames with scores_onceiling_resnet_v2 > 0

datadir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20260317_oncelingHiLowDays/';
plotdir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Locomotion_analysis/claudeplots_20260317_onceiling_hilow/';
if ~isfolder(plotdir), mkdir(plotdir); end

load(fullfile(datadir, 'perfly_onceiling_summary.mat'), 'fly_data');

% Layout: 3 rows (VNC, VNC2, VNC3) x 2 columns (low, high)
% No high-ceiling VNC day exists, so row 1 col 2 is empty
grid_dates = {
    '20210406', ''          % VNC:  low, no high
    '20220525', '20230706'  % VNC2: low, high
    '20240806', '20240724'  % VNC3: low, high
};
grid_labels = {
    '2021-04-06 (VNC, low)',   ''
    '2022-05-25 (VNC2, low)',  '2023-07-06 (VNC2, HIGH)'
    '2024-08-06 (VNC3, low)',  '2024-07-24 (VNC3, HIGH)'
};
row_labels = {'VNC', 'VNC2', 'VNC3'};

% Fixed axis limits across all panels
ylims = [0 100];
xlims = [9.5 16];

%% Plot 1: controls vs non-controls
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
            xlim(xlims); ylim(ylims);
            if col == 1, ylabel(row_labels{row}, 'FontWeight', 'bold'); end
            grid on;
            continue;
        end

        idx = strcmp({fly_data.date}, d);
        fd = fly_data(idx);

        nc = ~[fd.is_control];
        scatter([fd(nc).time_of_day], [fd(nc).pct_oc], 15, [0.6 0.6 0.6], ...
            'filled', 'MarkerFaceAlpha', 0.4);
        ct = [fd.is_control];
        scatter([fd(ct).time_of_day], [fd(ct).pct_oc], 25, [0.8 0.2 0.2], ...
            'filled', 'MarkerFaceAlpha', 0.7);

        xlim(xlims); ylim(ylims);
        xlabel('Time of day (hours)');
        ylabel('% onceiling per fly');
        title(grid_labels{row, col}, 'Interpreter', 'none');
        legend({'non-control', 'control (YNA)'}, 'Location', 'best');
        grid on;
    end
end

exportgraphics(gcf, fullfile(plotdir, 'perfly_onceiling_vs_timeofday_bydate.png'), 'Resolution', 150);
exportgraphics(gcf, fullfile(plotdir, 'perfly_onceiling_vs_timeofday_bydate.pdf'), 'ContentType', 'vector');
fprintf('Plot 1 saved.\n');

%% Plot 2: colored by rig
rig_colors = struct('A', [0.2 0.4 0.8], 'B', [0.8 0.3 0.2], ...
    'C', [0.2 0.7 0.3], 'D', [0.7 0.4 0.8]);
rig_list = {'A', 'B', 'C', 'D'};

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
            xlim(xlims); ylim(ylims);
            grid on;
            continue;
        end

        idx = strcmp({fly_data.date}, d);
        fd = fly_data(idx);

        for ri = 1:4
            r = rig_list{ri};
            ridx = strcmp({fd.rig}, r);
            scatter([fd(ridx).time_of_day], [fd(ridx).pct_oc], 15, rig_colors.(r), ...
                'filled', 'MarkerFaceAlpha', 0.5);
        end

        xlim(xlims); ylim(ylims);
        xlabel('Time of day (hours)');
        ylabel('% onceiling per fly');
        title(grid_labels{row, col}, 'Interpreter', 'none');
        legend(rig_list, 'Location', 'best');
        grid on;
    end
end

exportgraphics(gcf, fullfile(plotdir, 'perfly_onceiling_vs_timeofday_byrig.png'), 'Resolution', 150);
exportgraphics(gcf, fullfile(plotdir, 'perfly_onceiling_vs_timeofday_byrig.pdf'), 'ContentType', 'vector');
fprintf('Plot 2 saved.\n');

%% Plot 3: colored by sex
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
            xlim(xlims); ylim(ylims);
            grid on;
            continue;
        end

        idx = strcmp({fly_data.date}, d);
        fd = fly_data(idx);

        males = strcmp({fd.sex}, 'M');
        females = strcmp({fd.sex}, 'F');

        scatter([fd(males).time_of_day], [fd(males).pct_oc], 15, [0.2 0.4 0.8], ...
            'filled', 'MarkerFaceAlpha', 0.4);
        scatter([fd(females).time_of_day], [fd(females).pct_oc], 15, [0.8 0.3 0.5], ...
            'filled', 'MarkerFaceAlpha', 0.4);

        xlim(xlims); ylim(ylims);
        xlabel('Time of day (hours)');
        ylabel('% onceiling per fly');
        title(grid_labels{row, col}, 'Interpreter', 'none');
        legend({'Male', 'Female'}, 'Location', 'best');
        grid on;
    end
end

exportgraphics(gcf, fullfile(plotdir, 'perfly_onceiling_vs_timeofday_bysex.png'), 'Resolution', 150);
exportgraphics(gcf, fullfile(plotdir, 'perfly_onceiling_vs_timeofday_bysex.pdf'), 'ContentType', 'vector');
fprintf('Plot 3 saved.\n');

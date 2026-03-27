%% script_plot_walk_onceiling_overlap
% Explore walk/onceiling overlap for VNC2+VNC3 experiments with 25-90% onceiling.
% Goal: identify experiments where floor walks can be salvaged vs fully discarded.
%
% Per fly:
%   pct_walk_floor = % frames walking AND NOT onceiling (usable data)
%   pct_walk_ceil  = % frames walking AND onceiling (discard)
%   pct_walk       = total walking %
%   pct_oc         = total onceiling %

datafile = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/walk_onceiling_overlap_vnc23_25to90.mat';
plotdir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Locomotion_analysis/claudeplots_20260319_onceiling_allVNC/';
if ~isfolder(plotdir), mkdir(plotdir); end

load(datafile, 'exp_names', 'exp_pct_oc', 'exp_screen_type', 'fly_data');

n_exp = numel(exp_names);

% Compute per-experiment summaries from fly_data
exp_med_oc = nan(n_exp, 1);
exp_med_walk = nan(n_exp, 1);
exp_med_walk_floor = nan(n_exp, 1);
exp_med_walk_ceil = nan(n_exp, 1);
exp_total_walk_floor = nan(n_exp, 1);  % mean across flies
exp_total_walk = nan(n_exp, 1);
exp_n_flies_floor_walk = nan(n_exp, 1); % flies with >5% floor walking

for i = 1:n_exp
    idx = [fly_data.exp_i] == i;
    fd = fly_data(idx);
    exp_med_oc(i) = median([fd.pct_oc]);
    exp_med_walk(i) = median([fd.pct_walk]);
    exp_med_walk_floor(i) = median([fd.pct_walk_floor]);
    exp_med_walk_ceil(i) = median([fd.pct_walk_ceil]);
    exp_total_walk_floor(i) = mean([fd.pct_walk_floor]);
    exp_total_walk(i) = mean([fd.pct_walk]);
    exp_n_flies_floor_walk(i) = sum([fd.pct_walk_floor] > 5);
end

%% Plot 1: per-experiment scatter — % onceiling vs % floor walking
figure('Position', [100 100 800 600]);
hold on;
scatter(exp_pct_oc, exp_total_walk_floor, 15, [0.4 0.4 0.4], 'filled', 'MarkerFaceAlpha', 0.5);
xlabel('% onceiling (experiment average)', 'Interpreter', 'none');
ylabel('% floor walking (experiment average)', 'Interpreter', 'none');
title('Onceiling vs salvageable floor walking (VNC2+VNC3, 25-90% OC)', 'Interpreter', 'none');
xlim([25 90]); ylim([0 max(exp_total_walk_floor)*1.1]);
grid on;

exportgraphics(gcf, fullfile(plotdir, 'overlap_oc_vs_floor_walk.png'), 'Resolution', 150);
fprintf('Plot 1 saved.\n');

%% Plot 2: stacked bar-like — sorted by onceiling, show walk breakdown
[~, si] = sort(exp_pct_oc);

figure('Position', [100 100 1400 500]);
hold on;

% Floor walking
bar_floor = exp_total_walk_floor(si);
bar_ceil = exp_med_walk_ceil(si);

bar(1:n_exp, bar_floor, 1, 'FaceColor', [0.2 0.6 0.2], 'EdgeColor', 'none', 'FaceAlpha', 0.7);
bar(1:n_exp, -bar_ceil, 1, 'FaceColor', [0.8 0.2 0.2], 'EdgeColor', 'none', 'FaceAlpha', 0.7);

% Overlay onceiling % as line (scaled to fit)
plot(1:n_exp, exp_pct_oc(si), 'k-', 'LineWidth', 1);

xlabel('Experiments (sorted by % onceiling)');
ylabel('% of frames', 'Interpreter', 'none');
title('Green = floor walking (up), Red = ceiling walking (down), Black = % onceiling', 'Interpreter', 'none');
legend({'floor walking', 'ceiling walking', '% onceiling'}, 'Location', 'northwest');
xlim([0 n_exp+1]);
grid on;

exportgraphics(gcf, fullfile(plotdir, 'overlap_walk_breakdown_sorted.png'), 'Resolution', 150);
fprintf('Plot 2 saved.\n');

%% Plot 3: per-fly scatter — onceiling vs floor walking, one dot per fly
figure('Position', [100 100 800 600]);
hold on;
scatter([fly_data.pct_oc], [fly_data.pct_walk_floor], 4, [0.4 0.4 0.4], ...
    'filled', 'MarkerFaceAlpha', 0.15);
xlabel('% onceiling per fly', 'Interpreter', 'none');
ylabel('% floor walking per fly', 'Interpreter', 'none');
title('Per-fly: onceiling vs salvageable floor walking (n=11952)', 'Interpreter', 'none');
xlim([0 100]); ylim([0 60]);
grid on;

exportgraphics(gcf, fullfile(plotdir, 'overlap_perfly_oc_vs_floor_walk.png'), 'Resolution', 150);
fprintf('Plot 3 saved.\n');

%% Plot 4: how many flies per experiment have usable floor walking?
figure('Position', [100 100 800 600]);
hold on;
scatter(exp_pct_oc, exp_n_flies_floor_walk, 15, [0.4 0.4 0.4], 'filled', 'MarkerFaceAlpha', 0.5);
xlabel('% onceiling (experiment average)', 'Interpreter', 'none');
ylabel('Number of flies with >5% floor walking', 'Interpreter', 'none');
title('Flies with usable floor walks vs experiment onceiling level', 'Interpreter', 'none');
xlim([25 90]);
grid on;

exportgraphics(gcf, fullfile(plotdir, 'overlap_nflies_usable_vs_oc.png'), 'Resolution', 150);
fprintf('Plot 4 saved.\n');

%% Plot 5: % floor walking vs % ceiling walking per fly
figure('Position', [100 100 800 600]);
hold on;
scatter([fly_data.pct_walk_ceil], [fly_data.pct_walk_floor], 4, [0.4 0.4 0.4], ...
    'filled', 'MarkerFaceAlpha', 0.1);

% Reference: baseline walking from low-ceiling VNC2 (~28.6%)
yline(28.6, 'b--', 'LineWidth', 1);
text(1, 29.5, '  VNC2 baseline walking (28.6%)', 'Color', 'b', 'FontSize', 9);

% Identity line
plot([0 60], [0 60], 'k--', 'LineWidth', 0.5);

xlabel('% ceiling walking per fly', 'Interpreter', 'none');
ylabel('% floor walking per fly', 'Interpreter', 'none');
title(sprintf('Per-fly: floor vs ceiling walking (n=%d)', numel(fly_data)), 'Interpreter', 'none');
xlim([0 60]); ylim([0 60]);
legend({'flies', 'VNC2 baseline', 'identity'}, 'Location', 'northeast');
grid on;

exportgraphics(gcf, fullfile(plotdir, 'overlap_perfly_floor_vs_ceil_walking.png'), 'Resolution', 150);
fprintf('Plot 5a saved.\n');

%% Plot 5b: per-fly floor walking vs experiment % onceiling (2D histogram)
fly_exp_oc = exp_pct_oc([fly_data.exp_i]);
fly_floor_walk = [fly_data.pct_walk_floor];

x_edges = 25:1:90;
y_edges = 0:1:60;
counts = histcounts2(fly_exp_oc(:), fly_floor_walk(:), x_edges, y_edges);

figure('Position', [100 100 800 600]);
imagesc(x_edges(1:end-1)+0.5, y_edges(1:end-1)+0.5, counts');
set(gca, 'YDir', 'normal');
colormap(flipud(gray));
cb = colorbar;
cb.Label.String = 'Number of flies';
clim([0 prctile(counts(:), 99)]);  % clip top 1% for visibility

hold on;
yline(28.6, 'b--', 'LineWidth', 1.5);
text(26, 29.5, '  VNC2 baseline', 'Color', 'b', 'FontSize', 9);

xlabel('Experiment % onceiling', 'Interpreter', 'none');
ylabel('% floor walking per fly', 'Interpreter', 'none');
title(sprintf('Per-fly floor walking vs experiment onceiling (n=%d)', numel(fly_data)), 'Interpreter', 'none');
xlim([25 90]); ylim([0 60]);

exportgraphics(gcf, fullfile(plotdir, 'overlap_perfly_floorwalk_vs_exp_oc.png'), 'Resolution', 150);
fprintf('Plot 5b saved.\n');

%% Plot 5c: % floor walking vs % ceiling walking per experiment
figure('Position', [100 100 800 600]);
hold on;
scatter(exp_med_walk_ceil, exp_total_walk_floor, 15, [0.4 0.4 0.4], 'filled', 'MarkerFaceAlpha', 0.5);

% Reference: baseline walking from low-ceiling VNC2 (~28.6%)
yline(28.6, 'b--', 'LineWidth', 1);
text(1, 29.5, '  VNC2 baseline walking (28.6%)', 'Color', 'b', 'FontSize', 9);

% Identity line
maxval = max([exp_med_walk_ceil; exp_total_walk_floor]) * 1.05;
plot([0 maxval], [0 maxval], 'k--', 'LineWidth', 0.5);

xlabel('% ceiling walking (experiment avg)', 'Interpreter', 'none');
ylabel('% floor walking (experiment avg)', 'Interpreter', 'none');
title('Floor walking vs ceiling walking (VNC2+VNC3, 25-90% OC)', 'Interpreter', 'none');
xlim([0 maxval]); ylim([0 35]);
legend({'experiments', 'VNC2 baseline', 'identity'}, 'Location', 'northeast');
grid on;

exportgraphics(gcf, fullfile(plotdir, 'overlap_floor_vs_ceil_walking.png'), 'Resolution', 150);
fprintf('Plot 5 saved.\n');

%% Summary
fprintf('\nSummary (1216 experiments, 25-90%% onceiling):\n');
edges = [25 40 50 60 75 90];
for i = 1:numel(edges)-1
    in_bin = exp_pct_oc >= edges(i) & exp_pct_oc < edges(i+1);
    fprintf('  %d-%d%% OC: n=%d, median floor walk=%.1f%%, median flies w/>5%% floor walk=%.0f\n', ...
        edges(i), edges(i+1), sum(in_bin), ...
        median(exp_total_walk_floor(in_bin)), median(exp_n_flies_floor_walk(in_bin)));
end

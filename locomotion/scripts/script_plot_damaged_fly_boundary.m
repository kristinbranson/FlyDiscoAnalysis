%% script_plot_damaged_fly_boundary
% Plot OC% vs walking% for positive cagemate deviants (>30pp above cagemates).
% Highlights VGLUT flies and shows both total walking and LED-off-only walking.
%
% Purpose: identify damaged/smooshed flies and determine detection threshold,
% accounting for VGLUT immobilization during LED activation.

datafile = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/damaged_fly_deviants.mat';
plotdir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Locomotion_analysis/claudeplots_20260319_onceiling_allVNC/';

load(datafile, 'all_devs');

oc_vals = [all_devs.pct_oc];
wk_vals = [all_devs.pct_walk];
wk_ledoff = [all_devs.pct_walk_ledoff];
is_vglut = [all_devs.is_vglut];

%% Plot 1: OC vs total walking, VGLUT colored
figure('Position', [100 100 800 600]);
hold on;

% Non-VGLUT
scatter(oc_vals(~is_vglut), wk_vals(~is_vglut), 15, [0.5 0.5 0.5], ...
    'filled', 'MarkerFaceAlpha', 0.4);

% VGLUT
scatter(oc_vals(is_vglut), wk_vals(is_vglut), 25, [0.2 0.6 0.2], ...
    'filled', 'MarkerFaceAlpha', 0.8);

% Previously reviewed damaged (>80% OC, <2% walk, non-VGLUT)
reviewed = oc_vals > 80 & wk_vals < 2 & ~is_vglut;
scatter(oc_vals(reviewed), wk_vals(reviewed), 25, [0.8 0.2 0.2], 'filled');

yline(5, 'r--', 'LineWidth', 1);
yline(10, 'b--', 'LineWidth', 1);
text(32, 5.5, 'walk = 5%', 'Color', 'r');
text(32, 10.5, 'walk = 10%', 'Color', 'b');

xlabel('% onceiling per fly', 'Interpreter', 'none');
ylabel('% walking per fly (all frames)', 'Interpreter', 'none');
title(sprintf('Positive deviants (>30pp, n=%d): total walking', numel(all_devs)), 'Interpreter', 'none');
legend({'other lines', 'VGLUT', 'confirmed damaged', '5% walk', '10% walk'}, 'Location', 'northeast');
xlim([30 100]); ylim([0 50]);
grid on;

exportgraphics(gcf, fullfile(plotdir, 'damaged_fly_boundary_vglut.png'), 'Resolution', 150);
fprintf('Plot 1 saved.\n');

%% Plot 2: OC vs LED-off-only walking, VGLUT colored
figure('Position', [100 100 800 600]);
hold on;

scatter(oc_vals(~is_vglut), wk_ledoff(~is_vglut), 15, [0.5 0.5 0.5], ...
    'filled', 'MarkerFaceAlpha', 0.4);

scatter(oc_vals(is_vglut), wk_ledoff(is_vglut), 25, [0.2 0.6 0.2], ...
    'filled', 'MarkerFaceAlpha', 0.8);

reviewed = oc_vals > 80 & wk_vals < 2 & ~is_vglut;
scatter(oc_vals(reviewed), wk_ledoff(reviewed), 25, [0.8 0.2 0.2], 'filled');

yline(5, 'r--', 'LineWidth', 1);
yline(10, 'b--', 'LineWidth', 1);
text(32, 5.5, 'walk = 5%', 'Color', 'r');
text(32, 10.5, 'walk = 10%', 'Color', 'b');

xlabel('% onceiling per fly', 'Interpreter', 'none');
ylabel('% walking per fly (LED off only)', 'Interpreter', 'none');
title(sprintf('Positive deviants (>30pp, n=%d): LED-off walking', numel(all_devs)), 'Interpreter', 'none');
legend({'other lines', 'VGLUT', 'confirmed damaged', '5% walk', '10% walk'}, 'Location', 'northeast');
xlim([30 100]); ylim([0 50]);
grid on;

exportgraphics(gcf, fullfile(plotdir, 'damaged_fly_boundary_vglut_ledoff.png'), 'Resolution', 150);
fprintf('Plot 2 saved.\n');

%% Summary
fprintf('\nVGLUT deviants:\n');
vg = all_devs(is_vglut);
fprintf('  n = %d\n', numel(vg));
fprintf('  Total walk: median=%.1f%%, mean=%.1f%%\n', median([vg.pct_walk]), mean([vg.pct_walk]));
fprintf('  LED-off walk: median=%.1f%%, mean=%.1f%%\n', median([vg.pct_walk_ledoff]), mean([vg.pct_walk_ledoff]));
fprintf('  With <3%% total walk: %d\n', sum([vg.pct_walk] < 3));
fprintf('  With <3%% LED-off walk: %d\n', sum([vg.pct_walk_ledoff] < 3));

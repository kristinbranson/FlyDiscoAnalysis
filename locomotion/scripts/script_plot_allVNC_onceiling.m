%% script_plot_allVNC_onceiling
% Plot per-experiment onceiling % for all 8663 VNC/VNC2/VNC3 experiments.
% Each dot = one experiment (% of all fly-frames with scores_onceiling_resnet_v2 > 0).
%
% Plot 1: onceiling % over time, colored by screen type
% Plot 2: experiments sorted by onceiling %, colored by screen type

datafile = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/allVNC_onceiling_resnet_v2_summary.mat';
plotdir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Locomotion_analysis/claudeplots_20260319_onceiling_allVNC/';
if ~isfolder(plotdir), mkdir(plotdir); end

load(datafile, 'expnames', 'pct_oc', 'dates_dt', 'screen_type', 'is_control', 'fly_oc_per_exp');

is_vglut = contains(expnames, 'VGLUT');

colors = struct('VNC', [0.2 0.6 0.2], 'VNC2', [0.2 0.4 0.8], 'VNC3', [0.8 0.3 0.2]);
st_list = {'VNC', 'VNC2', 'VNC3'};

%% Plot 1: onceiling % over time
figure('Position', [100 100 1400 500]);
hold on;

for si = 1:3
    st = st_list{si};
    idx = strcmp(screen_type, st);
    scatter(dates_dt(idx), pct_oc(idx), 8, colors.(st), 'filled', 'MarkerFaceAlpha', 0.3);
end

xlabel('Date');
ylabel('% frames onceiling (ResNet v2)', 'Interpreter', 'none');
title('All VNC/VNC2/VNC3 experiments: onceiling % over time', 'Interpreter', 'none');
legend(st_list, 'Location', 'best');
ax = gca;
ax.XAxis.TickLabelFormat = 'MMM yyyy';
xtickangle(45);
ylim([0 100]);
grid on;

exportgraphics(gcf, fullfile(plotdir, 'allVNC_onceiling_over_time.png'), 'Resolution', 150);
exportgraphics(gcf, fullfile(plotdir, 'allVNC_onceiling_over_time.pdf'), 'ContentType', 'vector');
fprintf('Plot 1 saved.\n');

%% Plot 2: sorted by median onceiling %, with per-fly quartiles
% Compute per-experiment quartiles
n_exp = numel(pct_oc);
med_oc = nan(n_exp, 1);
q25_oc = nan(n_exp, 1);
q75_oc = nan(n_exp, 1);
for i = 1:n_exp
    fp = fly_oc_per_exp{i};
    med_oc(i) = median(fp);
    q25_oc(i) = prctile(fp, 25);
    q75_oc(i) = prctile(fp, 75);
end

[~, si] = sort(med_oc);
sorted_med = med_oc(si);
sorted_q25 = q25_oc(si);
sorted_q75 = q75_oc(si);
sorted_ctrl = is_control(si);
sorted_vglut = is_vglut(si);
sorted_other = ~sorted_ctrl & ~sorted_vglut;

figure('Position', [100 100 1400 500]);
hold on;

% IQR as filled region
x_fill = [1:n_exp, n_exp:-1:1];
y_fill = [sorted_q25(:)', fliplr(sorted_q75(:)')];
fill(x_fill, y_fill, [0.85 0.85 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.5);

% Median line
plot(1:n_exp, sorted_med, 'k-', 'LineWidth', 0.5);

% Median dots: other lines in grey
scatter(find(sorted_other), sorted_med(sorted_other), 8, [0.5 0.5 0.5], ...
    'filled', 'MarkerFaceAlpha', 0.4);

% Control (YNA) in red
scatter(find(sorted_ctrl), sorted_med(sorted_ctrl), 15, [0.8 0.2 0.2], ...
    'filled', 'MarkerFaceAlpha', 0.8);

% VGLUT in green
scatter(find(sorted_vglut), sorted_med(sorted_vglut), 15, [0.2 0.6 0.2], ...
    'filled', 'MarkerFaceAlpha', 0.8);

xlabel('Experiments (sorted by median % onceiling)');
ylabel('% onceiling per fly', 'Interpreter', 'none');
title('All VNC/VNC2/VNC3 sorted by median onceiling. Shading = IQR.', 'Interpreter', 'none');
legend({'IQR (Q25-Q75)', 'median', 'other lines', 'control (YNA)', 'VGLUT'}, 'Location', 'northwest');
xlim([0 n_exp+1]);
ylim([0 100]);
grid on;

% Mark threshold
n_above_50 = sum(sorted_med > 50);
if n_above_50 > 0
    thresh_x = n_exp - n_above_50 + 0.5;
    xline(thresh_x, 'r--', 'LineWidth', 1);
    text(thresh_x, 95, sprintf('  >50%% median: %d (%.1f%%)', ...
        n_above_50, 100*n_above_50/n_exp), 'Color', 'r', 'FontSize', 9);
end

exportgraphics(gcf, fullfile(plotdir, 'allVNC_onceiling_sorted.png'), 'Resolution', 150);
exportgraphics(gcf, fullfile(plotdir, 'allVNC_onceiling_sorted.pdf'), 'ContentType', 'vector');
fprintf('Plot 2 saved.\n');

%% Summary stats
fprintf('\nSummary:\n');
for si = 1:3
    st = st_list{si};
    idx = strcmp(screen_type, st);
    fprintf('  %s (n=%d): median=%.1f%%, mean=%.1f%%, >50%%=%d (%.1f%%)\n', ...
        st, sum(idx), median(pct_oc(idx)), mean(pct_oc(idx)), ...
        sum(pct_oc(idx) > 50), 100*mean(pct_oc(idx) > 50));
end
fprintf('  Total: n=%d, >50%%=%d (%.1f%%)\n', numel(pct_oc), ...
    sum(pct_oc > 50), 100*mean(pct_oc > 50));

%% script_plot_onceiling_threshold
% Experiments sorted by median % onceiling, showing individual fly values.
% Each column = one experiment, each dot = one fly's % onceiling.
% Red columns: experiments where >=50% of flies have >50% walking-on-ceiling.
%
% Purpose: help select a threshold for discarding high-ceiling experiments.

datadir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20260317_oncelingHiLowDays/';
plotdir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Locomotion_analysis/claudeplots_20260317_onceiling_hilow/';
if ~isfolder(plotdir), mkdir(plotdir); end

load(fullfile(datadir, 'experiment_threshold_data.mat'), 'exp_data');

% Sort by median onceiling
[~, si] = sort([exp_data.median_oc]);
exp_data = exp_data(si);
n_exp = numel(exp_data);

figure('Position', [100 100 1600 600]);
hold on;

for i = 1:n_exp
    fly_pcts = exp_data(i).fly_oc_pcts;
    x = repmat(i, 1, numel(fly_pcts));

    % Threshold: >=50% of flies have >50% onceiling
    frac_above = mean(fly_pcts > 50);
    if frac_above >= 0.50
        color = [0.8 0.2 0.2];
        alpha = 0.6;
    else
        color = [0.4 0.4 0.4];
        alpha = 0.4;
    end

    scatter(x, fly_pcts, 10, color, 'filled', 'MarkerFaceAlpha', alpha);
end

% Overlay median line
medians = [exp_data.median_oc];
plot(1:n_exp, medians, 'k-', 'LineWidth', 1);

% Recompute threshold per experiment after sorting
exceeds = false(1, n_exp);
for i = 1:n_exp
    exceeds(i) = mean(exp_data(i).fly_oc_pcts > 50) >= 0.50;
end
thresh_idx = find(exceeds, 1, 'first');
if ~isempty(thresh_idx)
    xline(thresh_idx - 0.5, 'r--', 'LineWidth', 1.5);
    n_discard = sum(exceeds);
    text(thresh_idx, 95, sprintf('  %d discarded (%.0f%%)', ...
        n_discard, 100*n_discard/n_exp), 'Color', 'r', 'FontSize', 10);
end

xlabel('Experiments (sorted by median % onceiling)');
ylabel('% onceiling per fly');
title('Per-fly onceiling % sorted by experiment median. Red = >=50% of flies have >50% onceiling', ...
    'Interpreter', 'none');
ylim([0 100]);
xlim([0 n_exp+1]);
grid on;

% Add legend
scatter(nan, nan, 30, [0.4 0.4 0.4], 'filled');
scatter(nan, nan, 30, [0.8 0.2 0.2], 'filled');
plot(nan, nan, 'k-', 'LineWidth', 1);
legend({'below threshold', 'above threshold (>=50% flies >50% OC)', 'median'}, ...
    'Location', 'northwest');

exportgraphics(gcf, fullfile(plotdir, 'onceiling_threshold_sorted.png'), 'Resolution', 150);
exportgraphics(gcf, fullfile(plotdir, 'onceiling_threshold_sorted.pdf'), 'ContentType', 'vector');
fprintf('Plot saved.\n');

% Print summary
fprintf('\nThreshold: >=50%% of flies with >50%% onceiling\n');
fprintf('Experiments exceeding threshold: %d / %d (%.1f%%)\n', ...
    sum(exceeds), n_exp, 100*mean(exceeds));
if ~isempty(thresh_idx)
    fprintf('Median OC at threshold boundary: %.1f%%\n', exp_data(thresh_idx).median_oc);
end

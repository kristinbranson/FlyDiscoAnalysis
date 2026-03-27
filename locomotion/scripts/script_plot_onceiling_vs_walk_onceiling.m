%% script_plot_onceiling_vs_walk_onceiling
% Scatter plot of per-experiment % onceiling vs % walking-while-onceiling
% for 624 YNA_K_162984 controls. Each dot = one experiment.
% Tests whether walking-while-onceiling adds information beyond onceiling alone.
%
% Metrics:
%   pct_onceiling: % of all fly-frames with scores_onceiling_resnet_v2 > 0
%   pct_walk_onceiling: % of all fly-frames with BOTH scores_Walk2 > 0 AND
%                       scores_onceiling_resnet_v2 > 0
% Both are pooled across all flies in an experiment (weighted by fly duration).

datadir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20260225_VNC23controldata4onceilingwalkexplore/';
plotdir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Locomotion_analysis/claudeplots_20260313_onceiling/';

load(fullfile(datadir, 'onceiling_resnet_v2_summary.mat'), ...
    'pct_onceiling', 'pct_walk_onceiling', 'screen_type');

if ~isfolder(plotdir), mkdir(plotdir); end

colors = struct('VNC', [0.2 0.6 0.2], 'VNC2', [0.2 0.4 0.8], 'VNC3', [0.8 0.3 0.2]);
st_list = {'VNC', 'VNC2', 'VNC3'};

% Correlations
[r, p] = corr(pct_onceiling, pct_walk_onceiling);
fprintf('Overall: r=%.3f, p=%.1e, n=%d\n', r, p, numel(pct_onceiling));
for si = 1:3
    st = st_list{si};
    idx = strcmp(screen_type, st);
    [ri, pi] = corr(pct_onceiling(idx), pct_walk_onceiling(idx));
    fprintf('  %s: r=%.3f, p=%.1e, n=%d\n', st, ri, pi, sum(idx));
end

% Scatter plot
figure('Position', [100 100 600 500]);
hold on;
for si = 1:3
    st = st_list{si};
    idx = strcmp(screen_type, st);
    scatter(pct_onceiling(idx), pct_walk_onceiling(idx), 25, colors.(st), 'filled', 'MarkerFaceAlpha', 0.6);
end

% Identity line and linear fit
xl = [0 100];
plot(xl, xl, 'k--', 'LineWidth', 0.5);
coeffs = polyfit(pct_onceiling, pct_walk_onceiling, 1);
xfit = linspace(0, max(pct_onceiling), 100);
plot(xfit, polyval(coeffs, xfit), 'k-', 'LineWidth', 1.5);

xlabel('% frames onceiling', 'Interpreter', 'none');
ylabel('% frames walking AND onceiling', 'Interpreter', 'none');
title(sprintf('onceiling vs walking-while-onceiling (r=%.3f)', r), 'Interpreter', 'none');
legend([st_list, {'identity', sprintf('fit: y=%.2fx%+.2f', coeffs(1), coeffs(2))}], 'Location', 'northwest');
axis equal;
xlim([0 100]); ylim([0 55]);
grid on;

exportgraphics(gcf, fullfile(plotdir, 'ctrl_onceiling_vs_walk_onceiling.png'));
exportgraphics(gcf, fullfile(plotdir, 'ctrl_onceiling_vs_walk_onceiling.pdf'), 'ContentType', 'vector');
fprintf('Plot saved.\n');

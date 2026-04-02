%% script_explore_onfloor_Z_counts
% Explore Z count distributions from locostatsperexp_onfloor.mat
% Only considers experiments with pct_onfloor >= 10% (not already excluded).

modpath;

%% Load experiment list and pct_onfloor
M = load('/groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/expdirs_VNC23_20251011.mat');
expdirs = {M.metadata.file_system_path}';
n = numel(expdirs);

OF = load('/groups/branson/bransonlab/flydisco_linelevel_VNC/pct_onfloor_VNC23.mat');

% Only process experiments passing the onfloor filter
idx_keep = OF.pct_onfloor >= 10 & OF.issuccess;
fprintf('Experiments passing pct_onfloor >= 10%%: %d/%d\n', sum(idx_keep), n);

%% Load Z counts for key stats
% Walk: Z = number of valid frames (pooled across flies)
% Bout (swing/step): Z = number of bouts (pooled across flies and limbs for "all")
stats_to_check = {
    'velmag_ctr__walk__LEDoff__all',           'walk_LEDoff'
    'velmag_ctr__walk__LEDon__all',            'walk_LEDon'
    'durations_time__swing__LEDoff__all',       'swing_LEDoff'
    'durations_time__swing__LEDon__all',        'swing_LEDon'
    'durations_time__step__LEDoff__all',        'step_LEDoff'
    'durations_time__step__LEDon__all',         'step_LEDon'
    'durations_time__swing__LEDoff__limb1',     'swing_LEDoff_limb1'
    'durations_time__swing__LEDon__limb1',      'swing_LEDon_limb1'
};
nstats = size(stats_to_check, 1);

fprintf('Loading Z counts from %d experiments...\n', sum(idx_keep));
Z = nan(n, nstats);
load_success = true(n, 1);

parfor i = 1:n
    if ~idx_keep(i), continue; end
    try
        S = load(fullfile(expdirs{i}, 'locostatsperexp_onfloor.mat'), 'locostatsperexp');
        stats = S.locostatsperexp;
        z_row = nan(1, nstats);
        for si = 1:nstats
            fn = stats_to_check{si, 1};
            if isfield(stats, fn)
                z_row(si) = stats.(fn).Z;
            end
        end
        Z(i,:) = z_row;
    catch
        load_success(i) = false;
    end
    if mod(i, 1000) == 0, fprintf('  %d/%d\n', i, n); end
end
fprintf('Done. %d/%d loaded.\n', sum(idx_keep & load_success), n);

%% Load if not already in memory
% if ~exist('Z', 'var')
%     load(...);
% end

%% Summary tables
idx_valid = idx_keep & load_success;
idx_nonvglut = idx_valid & ~OF.is_poscontrol;
idx_vglut = idx_valid & OF.is_poscontrol;
idx_control = idx_valid & OF.is_control;

for grp = {'Non-VGLUT', 'VGLUT', 'Control'}
    switch grp{1}
        case 'Non-VGLUT', idx = idx_nonvglut;
        case 'VGLUT', idx = idx_vglut;
        case 'Control', idx = idx_control;
    end
    fprintf('\n=== %s (%d experiments, pct_onfloor >= 10%%) ===\n', grp{1}, sum(idx));
    fprintf('%-30s %8s %8s %8s', 'Stat', 'Median', 'Min', 'Max');
    for thresh = [0 20 50 100 200 500]
        fprintf(' %7s', sprintf('Z<=%d', thresh));
    end
    fprintf('\n');
    fprintf('%s\n', repmat('-', 1, 110));

    for si = 1:nstats
        z = Z(idx, si);
        fprintf('%-30s %8.0f %8.0f %8.0f', stats_to_check{si,2}, ...
            median(z, 'omitnan'), min(z), max(z));
        for thresh = [0 20 50 100 200 500]
            fprintf(' %3d/%3d', sum(z <= thresh), numel(z));
        end
        fprintf('\n');
    end
end

%% Plots
plotdir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Locomotion_analysis/claudeplots_20260327_onfloor_spotchecks';
if ~isfolder(plotdir), mkdir(plotdir); end

% Plot 1: Z count vs pct_onfloor for walk LEDon (non-VGLUT)
figure('Position', [100 100 900 400]);
scatter(OF.pct_onfloor(idx_nonvglut), Z(idx_nonvglut, 2), 10, 'filled', 'MarkerFaceAlpha', 0.3);
xlabel('pct_onfloor (%)', 'Interpreter', 'none');
ylabel('Z count (walk frames LEDon)', 'Interpreter', 'none');
title('Walk LEDon Z count vs pct_onfloor (non-VGLUT, pct_onfloor>=10%)', 'Interpreter', 'none');
hold on;
yline(20, 'r--', 'LineWidth', 1.5);
yline(100, 'g--', 'LineWidth', 1.5);
yline(200, 'm--', 'LineWidth', 1.5);
legend({'Experiments', 'Z=20', 'Z=100', 'Z=200'}, 'Location', 'northwest', 'Interpreter', 'none');
saveas(gcf, fullfile(plotdir, 'Z_walk_LEDon_vs_pct_onfloor.png'));

% Plot 2: Z count vs pct_onfloor for swing LEDon per limb (non-VGLUT)
figure('Position', [100 100 900 400]);
scatter(OF.pct_onfloor(idx_nonvglut), Z(idx_nonvglut, 8), 10, 'filled', 'MarkerFaceAlpha', 0.3);
xlabel('pct_onfloor (%)', 'Interpreter', 'none');
ylabel('Z count (swing bouts LEDon, limb1)', 'Interpreter', 'none');
title('Swing LEDon limb1 Z vs pct_onfloor (non-VGLUT, pct_onfloor>=10%)', 'Interpreter', 'none');
hold on;
yline(20, 'r--', 'LineWidth', 1.5);
yline(100, 'g--', 'LineWidth', 1.5);
legend({'Experiments', 'Z=20', 'Z=100'}, 'Location', 'northwest', 'Interpreter', 'none');
saveas(gcf, fullfile(plotdir, 'Z_swing_LEDon_limb1_vs_pct_onfloor.png'));

% Plot 3: Histogram of walk LEDon Z (non-VGLUT, zoomed)
figure('Position', [100 100 900 400]);
z_nv = Z(idx_nonvglut, 2);
z_low = z_nv(z_nv < 2000);
histogram(z_low, 0:50:2000, 'FaceColor', [0.3 0.5 0.8]);
hold on;
xline(20, 'r--', 'LineWidth', 2);
xline(100, 'g--', 'LineWidth', 2);
xline(200, 'm--', 'LineWidth', 2);
xlabel('Z count (walk frames LEDon)', 'Interpreter', 'none');
ylabel('Number of experiments');
title('Walk LEDon Z distribution (non-VGLUT, pct_onfloor>=10%)', 'Interpreter', 'none');
legend({'Histogram', 'Z=20', 'Z=100', 'Z=200'}, 'Location', 'northeast');
saveas(gcf, fullfile(plotdir, 'Z_walk_LEDon_histogram.png'));

fprintf('\nPlots saved to %s\n', plotdir);

%% Fly-wise phase comparison: High day (Dec 6) vs Low day (Dec 5)
% Plot 1: Each row = one fly, x-axis = phase difference
% Color by day, saturation by LED on/off

clear fly_data fly_counter

modpath

%% Output directory
outputdir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Locomotion_analysis/claudeplots_20260211';

%% Setup
settingsdir = '/groups/branson/home/robiea/Code_versioned/BransonFlyDiscoSettings/settings';
analysis_protocol = '20251009_flybubble_LED_VNC2';

% Low day: Dec 5, 2023 (~-2.5 rad)
expdirs_low = {
    '/groups/branson/bransonlab/flydisco_data/VNC2_YNA_K_162984_RigA_20231205T114519', ...
    '/groups/branson/bransonlab/flydisco_data/VNC2_YNA_K_162984_RigB_20231205T114631', ...
    '/groups/branson/bransonlab/flydisco_data/VNC2_YNA_K_162984_RigC_20231205T114709', ...
    '/groups/branson/bransonlab/flydisco_data/VNC2_YNA_K_162984_RigD_20231205T114743'
};

% High day: Dec 6, 2023 (~-2.1 rad)
expdirs_high = {
    '/groups/branson/bransonlab/flydisco_data/VNC2_YNA_K_162984_RigA_20231206T125420', ...
    '/groups/branson/bransonlab/flydisco_data/VNC2_YNA_K_162984_RigB_20231206T125502', ...
    '/groups/branson/bransonlab/flydisco_data/VNC2_YNA_K_162984_RigC_20231206T125605', ...
    '/groups/branson/bransonlab/flydisco_data/VNC2_YNA_K_162984_RigD_20231206T125700'
};

all_expdirs = [expdirs_low, expdirs_high];
exp_labels = {'Dec5_RigA', 'Dec5_RigB', 'Dec5_RigC', 'Dec5_RigD', ...
              'Dec6_RigA', 'Dec6_RigB', 'Dec6_RigC', 'Dec6_RigD'};
day_labels = [repmat({'Dec5'}, 1, 4), repmat({'Dec6'}, 1, 4)];
rig_labels = {'A', 'B', 'C', 'D', 'A', 'B', 'C', 'D'};

% Read stage params
trx_temp = FBATrx('analysis_protocol', analysis_protocol, 'settingsdir', settingsdir, ...
    'datalocparamsfilestr', 'dataloc_params.txt');
stageparamsfile = fullfile(trx_temp.settingsdir, trx_temp.analysis_protocol, trx_temp.dataloc_params.locomotionmetricsparamsfilestr);
stage_params = ReadParams(stageparamsfile);
legtip_landmarknums = stage_params.legtip_landmarknums;
gc_threshold_low = stage_params.gc_threshold_low;
gc_threshold_high = stage_params.gc_threshold_high;
pairs = stage_params.pairs;
minimum_bout = stage_params.minimum_bout_groundcontact;

%% Collect per-fly phase data from all experiments (LED on and off)
fly_data = [];
fly_counter = 0;

for e = 1:numel(all_expdirs)
    expdir = all_expdirs{e};
    exp_label = exp_labels{e};
    day = day_labels{e};
    rig = rig_labels{e};

    fprintf('Processing %s...\n', exp_label);

    % Initialize trx
    trx = FBATrx('analysis_protocol', analysis_protocol, 'settingsdir', settingsdir, ...
        'datalocparamsfilestr', 'dataloc_params.txt');
    trx.AddExpDir(expdir, 'dooverwrite', false, 'openmovie', false);

    % Load data
    aptfile = trx.dataloc_params.apttrkfilestr;
    aptdata = TrkFile.load(fullfile(expdir, aptfile));

    load(fullfile(expdir, 'tips_velmag.mat'), 'tips_velmag');
    load(fullfile(expdir, 'tips_pos_body.mat'), 'tips_pos_body');

    if exist(fullfile(expdir, 'groundcontact.mat'), 'file')
        load(fullfile(expdir, 'groundcontact.mat'), 'groundcontact');
    else
        [groundcontact] = compute_groundcontact(tips_velmag, 'pairs', pairs, ...
            'gc_threshold_low', gc_threshold_low, 'gc_threshold_high', gc_threshold_high, ...
            'minimum_bout', minimum_bout);
    end

    [~, walking_scores] = LoadScoresFromFile(trx, 'scores_Walk2', 1);
    indicatordata = trx.getIndicatorLED(1);
    digitalindicator = indicatordata.indicatordigital;

    % Run LimbBoutAnalyzer
    loco_analyzer = LimbBoutAnalyzer(trx, aptdata, tips_pos_body, legtip_landmarknums, ...
        groundcontact, digitalindicator, walking_scores, ...
        'phase_methods', {'phasediff_hilbert'});

    loco_analyzer.analyzeWalkAndStimConditions();

    % Get number of flies
    nflies = numel(walking_scores);

    % Extract for both LED conditions
    led_conditions = {'led_off_traj', 'led_on_traj'};
    led_names = {'LEDoff', 'LEDon'};

    for f = 1:nflies
        fly_counter = fly_counter + 1;

        for led = 1:2
            led_cond = led_conditions{led};
            led_name = led_names{led};

            if ~loco_analyzer.walkMetrics.isKey(led_cond)
                fly_data(fly_counter).(led_name).phase_walk_mean = NaN;
                fly_data(fly_counter).(led_name).phase_frm_mean = NaN;
                fly_data(fly_counter).(led_name).phase_tripods_4 = NaN;
                fly_data(fly_counter).(led_name).n_walks = 0;
                fly_data(fly_counter).(led_name).n_frames = 0;
                continue
            end

            walk_metrics = loco_analyzer.walkMetrics(led_cond);
            perfly = walk_metrics.perfly;

            if f > numel(perfly) || ~isfield(perfly(f), 'phasediff_hilbert')
                fly_data(fly_counter).(led_name).phase_walk_mean = NaN;
                fly_data(fly_counter).(led_name).phase_frm_mean = NaN;
                fly_data(fly_counter).(led_name).phase_tripods_4 = NaN;
                fly_data(fly_counter).(led_name).n_walks = 0;
                fly_data(fly_counter).(led_name).n_frames = 0;
                continue
            end

            % ipsi_post_2 - walk means
            if isfield(perfly(f).phasediff_hilbert, 'ipsi_post_2') && ...
               isfield(perfly(f).phasediff_hilbert.ipsi_post_2, 'walk_mean_fly')
                fly_data(fly_counter).(led_name).phase_walk_mean = perfly(f).phasediff_hilbert.ipsi_post_2.walk_mean_fly;
                fly_data(fly_counter).(led_name).n_walks = perfly(f).phasediff_hilbert.ipsi_post_2.walk_n_fly;
            else
                fly_data(fly_counter).(led_name).phase_walk_mean = NaN;
                fly_data(fly_counter).(led_name).n_walks = 0;
            end

            % ipsi_post_2 - frame means
            if isfield(perfly(f).phasediff_hilbert, 'ipsi_post_2') && ...
               isfield(perfly(f).phasediff_hilbert.ipsi_post_2, 'frm_mean_fly')
                fly_data(fly_counter).(led_name).phase_frm_mean = perfly(f).phasediff_hilbert.ipsi_post_2.frm_mean_fly;
                fly_data(fly_counter).(led_name).n_frames = perfly(f).phasediff_hilbert.ipsi_post_2.frm_n_fly;
            else
                fly_data(fly_counter).(led_name).phase_frm_mean = NaN;
                fly_data(fly_counter).(led_name).n_frames = 0;
            end

            % tripods_4
            if isfield(perfly(f).phasediff_hilbert, 'tripods_4') && ...
               isfield(perfly(f).phasediff_hilbert.tripods_4, 'walk_mean_fly')
                fly_data(fly_counter).(led_name).phase_tripods_4 = perfly(f).phasediff_hilbert.tripods_4.walk_mean_fly;
            else
                fly_data(fly_counter).(led_name).phase_tripods_4 = NaN;
            end
        end

        fly_data(fly_counter).exp_idx = e;
        fly_data(fly_counter).exp_label = exp_label;
        fly_data(fly_counter).day = day;
        fly_data(fly_counter).rig = rig;
        fly_data(fly_counter).fly_in_exp = f;
    end

    fprintf('  %d flies extracted\n', nflies);
end

fprintf('\nTotal flies: %d\n', fly_counter);

%% PLOT 1a: Each row = one fly, x-axis = phase difference (WALK MEANS)
% Color by day, saturation by LED

figure('Position', [100, 100, 900, 900]);

% Colors: Dec5 = blue, Dec6 = red
% LED off = saturated, LED on = lighter
color_dec5_off = [0.1, 0.3, 0.8];   % dark blue
color_dec5_on = [0.6, 0.7, 0.95];   % light blue
color_dec6_off = [0.8, 0.1, 0.1];   % dark red
color_dec6_on = [0.95, 0.6, 0.6];   % light red

hold on;

% Reorder flies: Dec5 on top, Dec6 on bottom
dec5_idx = find(strcmp({fly_data.day}, 'Dec5'));
dec6_idx = find(strcmp({fly_data.day}, 'Dec6'));
fly_order = [dec5_idx, dec6_idx];  % Dec5 first (will be on top)

y_labels = {};
n_labels_off_x = -1.1;  % x position for LED off n labels
n_labels_on_x = -0.8;   % x position for LED on n labels

for i = 1:fly_counter
    fi = fly_order(i);
    y = fly_counter - i + 1;  % Reverse so first in list is at top

    day = fly_data(fi).day;

    % Determine colors
    if strcmp(day, 'Dec5')
        color_off = color_dec5_off;
        color_on = color_dec5_on;
    else
        color_off = color_dec6_off;
        color_on = color_dec6_on;
    end

    % Plot LED off (walk mean)
    phase_off = fly_data(fi).LEDoff.phase_walk_mean;
    n_off = fly_data(fi).LEDoff.n_walks;
    if ~isnan(phase_off)
        scatter(phase_off, y, 50, color_off, 'filled', 'MarkerFaceAlpha', 0.8);
        % n label at fixed x position
        text(n_labels_off_x, y, sprintf('%d', n_off), 'FontSize', 5, 'Color', color_off, 'HorizontalAlignment', 'right');
    end

    % Plot LED on (slightly offset in y for visibility)
    phase_on = fly_data(fi).LEDon.phase_walk_mean;
    n_on = fly_data(fi).LEDon.n_walks;
    if ~isnan(phase_on)
        scatter(phase_on, y + 0.3, 30, color_on, 'filled', 'MarkerFaceAlpha', 0.8);
        % n label at fixed x position
        text(n_labels_on_x, y + 0.3, sprintf('%d', n_on), 'FontSize', 5, 'Color', color_on, 'HorizontalAlignment', 'left');
    end

    % Connect LED off and on with a line
    if ~isnan(phase_off) && ~isnan(phase_on)
        plot([phase_off, phase_on], [y, y + 0.3], '-', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.5);
    end

    y_labels{i} = sprintf('%s_f%d', fly_data(fi).exp_label, fly_data(fi).fly_in_exp);
end

% Reverse y_labels to match reversed y positions
y_labels = flip(y_labels);

% Add separator between Dec5 and Dec6
n_dec5_flies = numel(dec5_idx);
yline(fly_counter - n_dec5_flies + 0.5, 'k--', 'LineWidth', 2);

% Add vertical lines for day means
dec5_phases_off = arrayfun(@(x) x.LEDoff.phase_walk_mean, fly_data(dec5_idx));
dec6_phases_off = arrayfun(@(x) x.LEDoff.phase_walk_mean, fly_data(dec6_idx));
mean_dec5 = circ_mean(dec5_phases_off(~isnan(dec5_phases_off))');
mean_dec6 = circ_mean(dec6_phases_off(~isnan(dec6_phases_off))');

xline(mean_dec5, '-', 'Color', color_dec5_off, 'LineWidth', 2, 'Label', sprintf('Dec5 mean: %.2f', mean_dec5), 'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'top', 'Interpreter', 'none');
xline(mean_dec6, '-', 'Color', color_dec6_off, 'LineWidth', 2, 'Label', sprintf('Dec6 mean: %.2f', mean_dec6), 'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'bottom', 'Interpreter', 'none');

xlabel('Phase difference ipsi_post_2 (rad)', 'Interpreter', 'none');
ylabel('Fly (experiment_fly#)', 'Interpreter', 'none');
xlim([-pi, 0]);
ylim([0, fly_counter + 1]);

set(gca, 'YTick', 1:fly_counter, 'YTickLabel', y_labels, 'FontSize', 6, 'TickLabelInterpreter', 'none');
t = title('Per-fly WALK MEAN phase: Dec5 (blue) vs Dec6 (red), n=walks', 'Interpreter', 'none');
t.Position(2) = t.Position(2) + 3;  % Move title up
grid on;

% Add column headers for n values
text(n_labels_off_x, fly_counter + 1.5, 'n(off)', 'FontSize', 7, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
text(n_labels_on_x, fly_counter + 1.5, 'n(on)', 'FontSize', 7, 'FontWeight', 'bold', 'HorizontalAlignment', 'left');

% Add legend with only 4 entries
h1 = scatter(NaN, NaN, 50, color_dec5_off, 'filled');
h2 = scatter(NaN, NaN, 30, color_dec5_on, 'filled');
h3 = scatter(NaN, NaN, 50, color_dec6_off, 'filled');
h4 = scatter(NaN, NaN, 30, color_dec6_on, 'filled');
legend([h1, h2, h3, h4], {'Dec5 LED off', 'Dec5 LED on', 'Dec6 LED off', 'Dec6 LED on'}, 'Location', 'southeast');

exportgraphics(gcf, fullfile(outputdir, 'plot1a_flywise_walk_mean.png'), 'Resolution', 150);

%% PLOT 1b: Each row = one fly, x-axis = phase difference (FRAME MEANS)
% Color by day, saturation by LED

figure('Position', [100, 100, 900, 900]);

hold on;

for i = 1:fly_counter
    fi = fly_order(i);
    y = fly_counter - i + 1;  % Reverse so first in list is at top

    day = fly_data(fi).day;

    % Determine colors
    if strcmp(day, 'Dec5')
        color_off = color_dec5_off;
        color_on = color_dec5_on;
    else
        color_off = color_dec6_off;
        color_on = color_dec6_on;
    end

    % Plot LED off (frame mean)
    phase_off = fly_data(fi).LEDoff.phase_frm_mean;
    n_off = fly_data(fi).LEDoff.n_frames;
    if ~isnan(phase_off)
        scatter(phase_off, y, 50, color_off, 'filled', 'MarkerFaceAlpha', 0.8);
        % n label at fixed x position
        text(n_labels_off_x, y, sprintf('%d', n_off), 'FontSize', 5, 'Color', color_off, 'HorizontalAlignment', 'right');
    end

    % Plot LED on (slightly offset in y for visibility)
    phase_on = fly_data(fi).LEDon.phase_frm_mean;
    n_on = fly_data(fi).LEDon.n_frames;
    if ~isnan(phase_on)
        scatter(phase_on, y + 0.3, 30, color_on, 'filled', 'MarkerFaceAlpha', 0.8);
        % n label at fixed x position
        text(n_labels_on_x, y + 0.3, sprintf('%d', n_on), 'FontSize', 5, 'Color', color_on, 'HorizontalAlignment', 'left');
    end

    % Connect LED off and on with a line
    if ~isnan(phase_off) && ~isnan(phase_on)
        plot([phase_off, phase_on], [y, y + 0.3], '-', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.5);
    end
end

% Add separator between Dec5 and Dec6
yline(fly_counter - n_dec5_flies + 0.5, 'k--', 'LineWidth', 2);

% Add vertical lines for day means (frame-based)
dec5_phases_frm = arrayfun(@(x) x.LEDoff.phase_frm_mean, fly_data(dec5_idx));
dec6_phases_frm = arrayfun(@(x) x.LEDoff.phase_frm_mean, fly_data(dec6_idx));
mean_dec5_frm = circ_mean(dec5_phases_frm(~isnan(dec5_phases_frm))');
mean_dec6_frm = circ_mean(dec6_phases_frm(~isnan(dec6_phases_frm))');

xline(mean_dec5_frm, '-', 'Color', color_dec5_off, 'LineWidth', 2, 'Label', sprintf('Dec5 mean: %.2f', mean_dec5_frm), 'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'top', 'Interpreter', 'none');
xline(mean_dec6_frm, '-', 'Color', color_dec6_off, 'LineWidth', 2, 'Label', sprintf('Dec6 mean: %.2f', mean_dec6_frm), 'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'bottom', 'Interpreter', 'none');

xlabel('Phase difference ipsi_post_2 (rad)', 'Interpreter', 'none');
ylabel('Fly (experiment_fly#)', 'Interpreter', 'none');
xlim([-pi, 0]);
ylim([0, fly_counter + 1]);

set(gca, 'YTick', 1:fly_counter, 'YTickLabel', y_labels, 'FontSize', 6, 'TickLabelInterpreter', 'none');
t = title('Per-fly FRAME MEAN phase: Dec5 (blue) vs Dec6 (red), n=frames', 'Interpreter', 'none');
t.Position(2) = t.Position(2) + 3;  % Move title up
grid on;

% Add column headers for n values
text(n_labels_off_x, fly_counter + 1.5, 'n(off)', 'FontSize', 7, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
text(n_labels_on_x, fly_counter + 1.5, 'n(on)', 'FontSize', 7, 'FontWeight', 'bold', 'HorizontalAlignment', 'left');

% Add legend with only 4 entries
h1 = scatter(NaN, NaN, 50, color_dec5_off, 'filled');
h2 = scatter(NaN, NaN, 30, color_dec5_on, 'filled');
h3 = scatter(NaN, NaN, 50, color_dec6_off, 'filled');
h4 = scatter(NaN, NaN, 30, color_dec6_on, 'filled');
legend([h1, h2, h3, h4], {'Dec5 LED off', 'Dec5 LED on', 'Dec6 LED off', 'Dec6 LED on'}, 'Location', 'southeast');

exportgraphics(gcf, fullfile(outputdir, 'plot1b_flywise_frame_mean.png'), 'Resolution', 150);

%% Summary stats
fprintf('\n=== Summary by Day (LED off, walk means) ===\n');
dec5_phases = [fly_data(strcmp({fly_data.day}, 'Dec5')).LEDoff];
dec5_phases = [dec5_phases.phase_walk_mean];
dec6_phases = [fly_data(strcmp({fly_data.day}, 'Dec6')).LEDoff];
dec6_phases = [dec6_phases.phase_walk_mean];

fprintf('Dec5: n=%d flies, mean=%.3f, std=%.3f\n', ...
    sum(~isnan(dec5_phases)), circ_mean(dec5_phases(~isnan(dec5_phases))'), circ_std(dec5_phases(~isnan(dec5_phases))'));
fprintf('Dec6: n=%d flies, mean=%.3f, std=%.3f\n', ...
    sum(~isnan(dec6_phases)), circ_mean(dec6_phases(~isnan(dec6_phases))'), circ_std(dec6_phases(~isnan(dec6_phases))'));
fprintf('Difference: %.3f rad\n', circ_dist(circ_mean(dec5_phases(~isnan(dec5_phases))'), circ_mean(dec6_phases(~isnan(dec6_phases))')));

fprintf('\n=== Summary by Day (LED off, frame means) ===\n');
dec5_phases_frm = [fly_data(strcmp({fly_data.day}, 'Dec5')).LEDoff];
dec5_phases_frm = [dec5_phases_frm.phase_frm_mean];
dec6_phases_frm = [fly_data(strcmp({fly_data.day}, 'Dec6')).LEDoff];
dec6_phases_frm = [dec6_phases_frm.phase_frm_mean];

fprintf('Dec5: n=%d flies, mean=%.3f, std=%.3f\n', ...
    sum(~isnan(dec5_phases_frm)), circ_mean(dec5_phases_frm(~isnan(dec5_phases_frm))'), circ_std(dec5_phases_frm(~isnan(dec5_phases_frm))'));
fprintf('Dec6: n=%d flies, mean=%.3f, std=%.3f\n', ...
    sum(~isnan(dec6_phases_frm)), circ_mean(dec6_phases_frm(~isnan(dec6_phases_frm))'), circ_std(dec6_phases_frm(~isnan(dec6_phases_frm))'));
fprintf('Difference: %.3f rad\n', circ_dist(circ_mean(dec5_phases_frm(~isnan(dec5_phases_frm))'), circ_mean(dec6_phases_frm(~isnan(dec6_phases_frm))')));

%% PLOT 2: Walks over time - each walk as a point, colored by phase
% Each row = one fly, x-axis = walk start time, color = phase difference

fprintf('\n=== Collecting per-walk data for Plot 2 ===\n');

% Collect per-walk data
walk_data = [];
walk_counter = 0;
global_fly_counter = 0;

for e = 1:numel(all_expdirs)
    expdir = all_expdirs{e};
    exp_label = exp_labels{e};
    day = day_labels{e};

    fprintf('Processing %s for walk data...\n', exp_label);

    % Initialize trx
    trx = FBATrx('analysis_protocol', analysis_protocol, 'settingsdir', settingsdir, ...
        'datalocparamsfilestr', 'dataloc_params.txt');
    trx.AddExpDir(expdir, 'dooverwrite', false, 'openmovie', false);

    % Load data
    aptfile = trx.dataloc_params.apttrkfilestr;
    aptdata = TrkFile.load(fullfile(expdir, aptfile));

    load(fullfile(expdir, 'tips_velmag.mat'), 'tips_velmag');
    load(fullfile(expdir, 'tips_pos_body.mat'), 'tips_pos_body');

    if exist(fullfile(expdir, 'groundcontact.mat'), 'file')
        load(fullfile(expdir, 'groundcontact.mat'), 'groundcontact');
    else
        [groundcontact] = compute_groundcontact(tips_velmag, 'pairs', pairs, ...
            'gc_threshold_low', gc_threshold_low, 'gc_threshold_high', gc_threshold_high, ...
            'minimum_bout', minimum_bout);
    end

    [~, walking_scores] = LoadScoresFromFile(trx, 'scores_Walk2', 1);
    indicatordata = trx.getIndicatorLED(1);
    digitalindicator = indicatordata.indicatordigital;

    % Run LimbBoutAnalyzer
    loco_analyzer = LimbBoutAnalyzer(trx, aptdata, tips_pos_body, legtip_landmarknums, ...
        groundcontact, digitalindicator, walking_scores, ...
        'phase_methods', {'phasediff_hilbert'});

    loco_analyzer.analyzeWalkAndStimConditions();

    nflies = numel(walking_scores);

    % Only LED off for simplicity
    if ~loco_analyzer.walkMetrics.isKey('led_off_traj')
        continue
    end

    walk_metrics = loco_analyzer.walkMetrics('led_off_traj');
    perwalk = walk_metrics.perwalk;

    for w = 1:numel(perwalk)
        if ~isfield(perwalk(w), 'phasediff_hilbert'), continue; end
        if ~isfield(perwalk(w).phasediff_hilbert, 'ipsi_post_2'), continue; end

        walk_counter = walk_counter + 1;

        walk_data(walk_counter).exp_idx = e;
        walk_data(walk_counter).exp_label = exp_label;
        walk_data(walk_counter).day = day;
        walk_data(walk_counter).fly = perwalk(w).fly;
        walk_data(walk_counter).walk_t0 = perwalk(w).walk_t0;
        walk_data(walk_counter).walk_t1 = perwalk(w).walk_t1;
        walk_data(walk_counter).phase_mean = perwalk(w).phasediff_hilbert.ipsi_post_2.mean;

        % Create global fly ID (unique across experiments)
        walk_data(walk_counter).global_fly_id = global_fly_counter + perwalk(w).fly;
    end

    global_fly_counter = global_fly_counter + nflies;
    fprintf('  %d walks extracted\n', numel(perwalk));
end

fprintf('Total walks: %d\n', walk_counter);

%% Create Plot 2
figure('Position', [100, 100, 1200, 900]);

% Colormap settings - blue (more negative) to red (less negative)
% Create custom blue-white-red colormap
n_colors = 256;
half = n_colors / 2;
blue_to_white = [linspace(0, 1, half)', linspace(0, 1, half)', ones(half, 1)];
white_to_red = [ones(half, 1), linspace(1, 0, half)', linspace(1, 0, half)'];
cmap = [blue_to_white; white_to_red];
clim_range = [-2.6, -2.1];

% Create unique fly list in same order as Plot 1
unique_flies = [];
fly_labels_plot2 = {};
fly_days = {};

for e = 1:numel(all_expdirs)
    exp_label = exp_labels{e};
    day = day_labels{e};

    % Find walks from this experiment
    exp_walks = walk_data(strcmp({walk_data.exp_label}, exp_label));
    if isempty(exp_walks), continue; end

    flies_in_exp = unique([exp_walks.fly]);
    for f = flies_in_exp
        unique_flies(end+1).exp_label = exp_label;
        unique_flies(end).fly = f;
        unique_flies(end).day = day;
        fly_labels_plot2{end+1} = sprintf('%s_f%d', exp_label, f);
        fly_days{end+1} = day;
    end
end

% Reorder: Dec5 on top
dec5_fly_idx = find(strcmp(fly_days, 'Dec5'));
dec6_fly_idx = find(strcmp(fly_days, 'Dec6'));
fly_order_plot2 = [dec5_fly_idx, dec6_fly_idx];

nflies_plot2 = numel(unique_flies);

hold on;

% Plot each walk
for i = 1:nflies_plot2
    fi = fly_order_plot2(i);
    y = nflies_plot2 - i + 1;  % Reverse so first is at top

    exp_label = unique_flies(fi).exp_label;
    fly_id = unique_flies(fi).fly;

    % Find walks for this fly
    fly_walks_idx = find(strcmp({walk_data.exp_label}, exp_label) & [walk_data.fly] == fly_id);

    for wi = fly_walks_idx
        t0 = walk_data(wi).walk_t0;
        phase = walk_data(wi).phase_mean;

        if isnan(phase), continue; end

        % Map phase to color
        phase_norm = (phase - clim_range(1)) / (clim_range(2) - clim_range(1));
        phase_norm = max(0, min(1, phase_norm));  % Clamp to [0,1]
        color_idx = round(phase_norm * 255) + 1;
        color_idx = max(1, min(256, color_idx));

        scatter(t0, y, 20, cmap(color_idx, :), 'filled', 'MarkerFaceAlpha', 0.7);
    end
end

% Reorder labels to match plot
fly_labels_ordered = fly_labels_plot2(fly_order_plot2);
fly_labels_ordered = flip(fly_labels_ordered);

% Add separator between Dec5 and Dec6
n_dec5_flies_plot2 = numel(dec5_fly_idx);
yline(nflies_plot2 - n_dec5_flies_plot2 + 0.5, 'k--', 'LineWidth', 2);

% Add day labels
text(-500, nflies_plot2 - n_dec5_flies_plot2/2, 'Dec5', 'FontSize', 12, 'FontWeight', 'bold', 'Color', color_dec5_off, 'HorizontalAlignment', 'right');
text(-500, (nflies_plot2 - n_dec5_flies_plot2)/2, 'Dec6', 'FontSize', 12, 'FontWeight', 'bold', 'Color', color_dec6_off, 'HorizontalAlignment', 'right');

xlabel('Walk start time (frames)', 'Interpreter', 'none');
ylabel('Fly', 'Interpreter', 'none');
ylim([0, nflies_plot2 + 1]);

set(gca, 'YTick', 1:nflies_plot2, 'YTickLabel', fly_labels_ordered, 'FontSize', 5, 'TickLabelInterpreter', 'none');
title('Per-walk phase over time: color = phase difference', 'Interpreter', 'none');

% Colorbar
colormap(cmap);
caxis(clim_range);
cb = colorbar;
cb.Label.String = 'Phase diff ipsi_post_2 (rad)';
cb.Label.Interpreter = 'none';

exportgraphics(gcf, fullfile(outputdir, 'plot2_walks_over_time.png'), 'Resolution', 150);

fprintf('\nPlot 2 saved.\n');

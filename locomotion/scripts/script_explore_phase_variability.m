%% Explore phase difference variability across walks and flies
% Investigate sources of bimodality in phase difference data
% Compare Dec 5 vs Dec 6 experiments, LED ON vs LED OFF

modpath

%% Output directory for saving figures
outputdir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Locomotion_analysis/claudeplots_20260211';
if ~exist(outputdir, 'dir')
    mkdir(outputdir);
end

%% Setup - experiments to compare
settingsdir = '/groups/branson/home/robiea/Code_versioned/BransonFlyDiscoSettings/settings';
analysis_protocol = '20251009_flybubble_LED_VNC2';

% Experiments from Dec 5 and Dec 6
expdirs = {
    '/groups/branson/bransonlab/flydisco_data/VNC2_YNA_K_162984_RigA_20231205T114519', ...
    '/groups/branson/bransonlab/flydisco_data/VNC2_YNA_K_162984_RigC_20231206T125605'
};
exp_labels = {'Dec5_RigA', 'Dec6_RigC'};

% LED conditions to analyze
led_conditions = {'led_off_traj', 'led_on_traj'};
led_labels = {'LEDoff', 'LEDon'};

% Read stage params (same for all experiments)
trx_temp = FBATrx('analysis_protocol', analysis_protocol, 'settingsdir', settingsdir, ...
    'datalocparamsfilestr', 'dataloc_params.txt');
stageparamsfile = fullfile(trx_temp.settingsdir, trx_temp.analysis_protocol, trx_temp.dataloc_params.locomotionmetricsparamsfilestr);
stage_params = ReadParams(stageparamsfile);
legtip_landmarknums = stage_params.legtip_landmarknums;
gc_threshold_low = stage_params.gc_threshold_low;
gc_threshold_high = stage_params.gc_threshold_high;
pairs = stage_params.pairs;
minimum_bout = stage_params.minimum_bout_groundcontact;

%% Loop over experiments and conditions to collect data
all_walk_data = struct();
condition_idx = 0;

for e = 1:numel(expdirs)
    expdir = expdirs{e};
    exp_label = exp_labels{e};

    fprintf('\n=== Processing %s ===\n', exp_label);

    % Initialize trx
    trx = FBATrx('analysis_protocol', analysis_protocol, 'settingsdir', settingsdir, ...
        'datalocparamsfilestr', 'dataloc_params.txt');
    trx.AddExpDir(expdir, 'dooverwrite', false, 'openmovie', false);

    % Load APT data
    aptfile = trx.dataloc_params.apttrkfilestr;
    aptdata = TrkFile.load(fullfile(expdir, aptfile));

    % Load cached data
    if exist(fullfile(expdir, 'tips_velmag.mat'), 'file')
        load(fullfile(expdir, 'tips_velmag.mat'), 'tips_velmag');
    else
        error('tips_velmag.mat not found in %s', expdir);
    end

    if exist(fullfile(expdir, 'tips_pos_body.mat'), 'file')
        load(fullfile(expdir, 'tips_pos_body.mat'), 'tips_pos_body');
    else
        error('tips_pos_body.mat not found in %s', expdir);
    end

    % Load or compute groundcontact
    if exist(fullfile(expdir, 'groundcontact.mat'), 'file')
        load(fullfile(expdir, 'groundcontact.mat'), 'groundcontact');
    else
        [groundcontact] = compute_groundcontact(tips_velmag, 'pairs', pairs, ...
            'gc_threshold_low', gc_threshold_low, 'gc_threshold_high', gc_threshold_high, ...
            'minimum_bout', minimum_bout);
    end

    % Walking scores and LED indicator
    [~, walking_scores] = LoadScoresFromFile(trx, 'scores_Walk2', 1);
    indicatordata = trx.getIndicatorLED(1);
    digitalindicator = indicatordata.indicatordigital;

    % Initialize LimbBoutAnalyzer
    loco_analyzer = LimbBoutAnalyzer(trx, aptdata, tips_pos_body, legtip_landmarknums, ...
        groundcontact, digitalindicator, walking_scores, ...
        'phase_methods', {'phasediff_hilbert'});

    % Run analysis
    loco_analyzer.analyzeWalkAndStimConditions();

    % Extract data for each LED condition
    for led = 1:numel(led_conditions)
        led_cond = led_conditions{led};
        led_label = led_labels{led};

        condition_idx = condition_idx + 1;
        cond_label = sprintf('%s_%s', exp_label, led_label);

        fprintf('  Extracting %s...\n', cond_label);

        % Get walk metrics
        if ~loco_analyzer.walkMetrics.isKey(led_cond)
            fprintf('    Warning: %s not found, skipping\n', led_cond);
            continue
        end
        walk_metrics = loco_analyzer.walkMetrics(led_cond);
        perwalk = walk_metrics.perwalk;
        nwalks = numel(perwalk);

        if nwalks == 0
            fprintf('    No walks found for %s\n', cond_label);
            continue
        end

        % Build walk data struct
        walk_data = struct();
        walk_data.condition = repmat({cond_label}, nwalks, 1);
        walk_data.exp_label = repmat({exp_label}, nwalks, 1);
        walk_data.led_label = repmat({led_label}, nwalks, 1);
        walk_data.fly = [perwalk.fly]';
        walk_data.walk_t0 = [perwalk.walk_t0]';
        walk_data.walk_t1 = [perwalk.walk_t1]';
        walk_data.walk_duration = walk_data.walk_t1 - walk_data.walk_t0;

        % Extract speed
        walk_data.speed_mean = nan(nwalks, 1);
        for w = 1:nwalks
            if isfield(perwalk(w), 'velmag_ctr') && isfield(perwalk(w).velmag_ctr, 'mean')
                walk_data.speed_mean(w) = perwalk(w).velmag_ctr.mean;
            end
        end

        % Phase pairs and groups
        phase_pairs = {'LM_LH', 'RM_RH', 'LF_LM', 'RF_RM'};
        phase_groups = {'ipsi_post_2', 'ipsi_ant_2', 'ipsi_P2A_4', 'tripods_4'};

        % Initialize phase arrays
        for pp = 1:numel(phase_pairs)
            pairname = phase_pairs{pp};
            walk_data.(['phase_' pairname '_mean']) = nan(nwalks, 1);
        end
        for pg = 1:numel(phase_groups)
            groupname = phase_groups{pg};
            walk_data.(['phase_' groupname '_mean']) = nan(nwalks, 1);
        end

        % Extract phase data
        for w = 1:nwalks
            if ~isfield(perwalk(w), 'phasediff_hilbert')
                continue
            end
            ph = perwalk(w).phasediff_hilbert;

            for pp = 1:numel(phase_pairs)
                pairname = phase_pairs{pp};
                if isfield(ph, pairname) && isfield(ph.(pairname), 'mean')
                    walk_data.(['phase_' pairname '_mean'])(w) = ph.(pairname).mean;
                end
            end

            for pg = 1:numel(phase_groups)
                groupname = phase_groups{pg};
                if isfield(ph, groupname) && isfield(ph.(groupname), 'mean')
                    walk_data.(['phase_' groupname '_mean'])(w) = ph.(groupname).mean;
                end
            end
        end

        % Store in all_walk_data
        all_walk_data(condition_idx).label = cond_label;
        all_walk_data(condition_idx).exp_label = exp_label;
        all_walk_data(condition_idx).led_label = led_label;
        all_walk_data(condition_idx).data = walk_data;
        all_walk_data(condition_idx).nwalks = nwalks;

        fprintf('    %d walks extracted\n', nwalks);
    end
end

%% Combine all data into one table
all_tables = {};
for c = 1:numel(all_walk_data)
    if ~isempty(all_walk_data(c).data)
        all_tables{end+1} = struct2table(all_walk_data(c).data);
    end
end
T_all = vertcat(all_tables{:});
fprintf('\nTotal walks across all conditions: %d\n', height(T_all));

%% Define colors for conditions
cmap = lines(4);
cond_colors = containers.Map();
cond_colors('Dec5_RigA_LEDoff') = cmap(1,:);
cond_colors('Dec5_RigA_LEDon') = cmap(2,:);
cond_colors('Dec6_RigC_LEDoff') = cmap(3,:);
cond_colors('Dec6_RigC_LEDon') = cmap(4,:);

%% PLOT 1: Distribution comparison - stair histograms for clarity
figure('Position', [100, 100, 1400, 600]);
tiledlayout(2, 3, 'TileSpacing', 'compact');

% Fixed axis limits
xlims = [-pi, pi];
nbins = 30;
bin_edges = linspace(-pi, pi, nbins+1);

% 1a. ipsi_post_2 - all 4 conditions overlaid
nexttile;
hold on;
for c = 1:numel(all_walk_data)
    if isempty(all_walk_data(c).data), continue; end
    data = all_walk_data(c).data.phase_ipsi_post_2_mean;
    histogram(data, bin_edges, 'DisplayStyle', 'stairs', 'LineWidth', 2, ...
        'EdgeColor', cond_colors(all_walk_data(c).label), 'DisplayName', all_walk_data(c).label);
end
xlim(xlims);
xlabel('Phase difference (rad)');
ylabel('Count');
title('ipsi\_post\_2 (all conditions)');
legend('Location', 'best');

% 1b. LM_LH (left posterior) - all conditions
nexttile;
hold on;
for c = 1:numel(all_walk_data)
    if isempty(all_walk_data(c).data), continue; end
    data = all_walk_data(c).data.phase_LM_LH_mean;
    histogram(data, bin_edges, 'DisplayStyle', 'stairs', 'LineWidth', 2, ...
        'EdgeColor', cond_colors(all_walk_data(c).label), 'DisplayName', all_walk_data(c).label);
end
xlim(xlims);
xlabel('Phase difference (rad)');
ylabel('Count');
title('LM\_LH (left posterior)');
legend('Location', 'best');

% 1c. RM_RH (right posterior) - all conditions
nexttile;
hold on;
for c = 1:numel(all_walk_data)
    if isempty(all_walk_data(c).data), continue; end
    data = all_walk_data(c).data.phase_RM_RH_mean;
    histogram(data, bin_edges, 'DisplayStyle', 'stairs', 'LineWidth', 2, ...
        'EdgeColor', cond_colors(all_walk_data(c).label), 'DisplayName', all_walk_data(c).label);
end
xlim(xlims);
xlabel('Phase difference (rad)');
ylabel('Count');
title('RM\_RH (right posterior)');
legend('Location', 'best');

% 1d. tripods_4 - all conditions
nexttile;
hold on;
for c = 1:numel(all_walk_data)
    if isempty(all_walk_data(c).data), continue; end
    data = all_walk_data(c).data.phase_tripods_4_mean;
    histogram(data, bin_edges, 'DisplayStyle', 'stairs', 'LineWidth', 2, ...
        'EdgeColor', cond_colors(all_walk_data(c).label), 'DisplayName', all_walk_data(c).label);
end
xlim(xlims);
xlabel('Phase difference (rad)');
ylabel('Count');
title('tripods\_4');
legend('Location', 'best');

% 1e. ipsi_ant_2 (anterior pairs) - all conditions
nexttile;
hold on;
for c = 1:numel(all_walk_data)
    if isempty(all_walk_data(c).data), continue; end
    data = all_walk_data(c).data.phase_ipsi_ant_2_mean;
    histogram(data, bin_edges, 'DisplayStyle', 'stairs', 'LineWidth', 2, ...
        'EdgeColor', cond_colors(all_walk_data(c).label), 'DisplayName', all_walk_data(c).label);
end
xlim(xlims);
xlabel('Phase difference (rad)');
ylabel('Count');
title('ipsi\_ant\_2 (LF\_LM + RF\_RM)');
legend('Location', 'best');

% 1f. LF_LM and RF_RM individually
nexttile;
hold on;
for c = 1:numel(all_walk_data)
    if isempty(all_walk_data(c).data), continue; end
    data = all_walk_data(c).data.phase_LF_LM_mean;
    histogram(data, bin_edges, 'DisplayStyle', 'stairs', 'LineWidth', 2, ...
        'EdgeColor', cond_colors(all_walk_data(c).label), 'DisplayName', [all_walk_data(c).label ' LF\_LM']);
end
xlim(xlims);
xlabel('Phase difference (rad)');
ylabel('Count');
title('LF\_LM (left anterior)');
legend('Location', 'best');

sgtitle('Phase Difference Distributions - Dec 5 vs Dec 6, LED OFF vs ON');
exportgraphics(gcf, fullfile(outputdir, 'plot1_distributions.png'), 'Resolution', 150);

%% PLOT 2: Summary bar plot with means and stds
figure('Position', [100, 100, 1200, 500]);
tiledlayout(1, 3);

phase_features = {'ipsi_post_2', 'ipsi_ant_2', 'tripods_4'};
feature_titles = {'ipsi\_post\_2', 'ipsi\_ant\_2', 'tripods\_4'};

for pf = 1:numel(phase_features)
    nexttile;
    hold on;

    means = nan(numel(all_walk_data), 1);
    stds = nan(numel(all_walk_data), 1);
    labels = {};

    for c = 1:numel(all_walk_data)
        if isempty(all_walk_data(c).data), continue; end
        data = all_walk_data(c).data.(['phase_' phase_features{pf} '_mean']);
        valid = data(~isnan(data));
        means(c) = circ_mean(valid);
        stds(c) = circ_std(valid);
        labels{c} = all_walk_data(c).label;
    end

    bar_colors = zeros(numel(all_walk_data), 3);
    for c = 1:numel(all_walk_data)
        bar_colors(c,:) = cond_colors(all_walk_data(c).label);
    end

    for c = 1:numel(all_walk_data)
        bar(c, means(c), 'FaceColor', bar_colors(c,:));
    end
    errorbar(1:numel(all_walk_data), means, stds, 'k.', 'LineWidth', 1.5);

    ylim(xlims);
    ylabel('Phase difference (rad)');
    set(gca, 'XTick', 1:numel(labels), 'XTickLabel', labels, 'XTickLabelRotation', 45);
    title(feature_titles{pf});
end

sgtitle('Mean Phase Differences by Condition (circ\_mean +/- circ\_std)');
exportgraphics(gcf, fullfile(outputdir, 'plot2_means_bars.png'), 'Resolution', 150);

%% PLOT 3: Phase vs Speed by condition
figure('Position', [100, 100, 1200, 400]);
tiledlayout(1, 4);

for c = 1:numel(all_walk_data)
    if isempty(all_walk_data(c).data), continue; end

    nexttile;
    scatter(all_walk_data(c).data.speed_mean, all_walk_data(c).data.phase_ipsi_post_2_mean, ...
        20, cond_colors(all_walk_data(c).label), 'filled', 'MarkerFaceAlpha', 0.6);

    xlabel('Speed (mm/s)');
    ylabel('Phase diff (rad)');
    ylim(xlims);
    title(all_walk_data(c).label, 'Interpreter', 'none');

    % Add correlation
    valid = ~isnan(all_walk_data(c).data.speed_mean) & ~isnan(all_walk_data(c).data.phase_ipsi_post_2_mean);
    if sum(valid) > 2
        [r, p] = corr(all_walk_data(c).data.speed_mean(valid), ...
            all_walk_data(c).data.phase_ipsi_post_2_mean(valid), 'type', 'Spearman');
        text(0.05, 0.95, sprintf('r=%.2f\np=%.2g', r, p), 'Units', 'normalized', ...
            'VerticalAlignment', 'top', 'FontSize', 9);
    end
end

sgtitle('ipsi\_post\_2 Phase vs Speed');
exportgraphics(gcf, fullfile(outputdir, 'plot3_phase_vs_speed.png'), 'Resolution', 150);

%% PLOT 4: Left vs Right comparison by condition
figure('Position', [100, 100, 1200, 400]);
tiledlayout(1, 4);

for c = 1:numel(all_walk_data)
    if isempty(all_walk_data(c).data), continue; end

    nexttile;
    hold on;

    L_data = all_walk_data(c).data.phase_LM_LH_mean;
    R_data = all_walk_data(c).data.phase_RM_RH_mean;

    scatter(L_data, R_data, 20, all_walk_data(c).data.fly, 'filled', 'MarkerFaceAlpha', 0.6);
    plot(xlims, xlims, 'k--', 'LineWidth', 1);

    xlabel('LM\_LH (left)');
    ylabel('RM\_RH (right)');
    xlim(xlims);
    ylim(xlims);
    axis square;
    title(all_walk_data(c).label, 'Interpreter', 'none');
    colorbar;

    % Compute L-R difference
    L_mean = circ_mean(L_data(~isnan(L_data)));
    R_mean = circ_mean(R_data(~isnan(R_data)));
    LR_diff = circ_dist(L_mean, R_mean);
    text(0.05, 0.95, sprintf('L-R=%.3f', LR_diff), 'Units', 'normalized', ...
        'VerticalAlignment', 'top', 'FontSize', 9);
end

sgtitle('Left vs Right Posterior Phase (colored by fly)');
exportgraphics(gcf, fullfile(outputdir, 'plot4_left_vs_right.png'), 'Resolution', 150);

%% PLOT 5: Raw phase LM vs RM (legs 5 and 2)
% Leg order: {'RF','RM','RH','LH','LM','LF'} = legs 1-6
% LM = leg 5, RM = leg 2

figure('Position', [100, 100, 1400, 800]);
tiledlayout(2, 4, 'TileSpacing', 'compact');

% Extract raw phase data from perwalk for each condition
for c = 1:numel(all_walk_data)
    if isempty(all_walk_data(c).data), continue; end

    % Get perwalk from walk_metrics (need to re-access)
    % We stored walk_metrics earlier, but let's extract phasedata
end

% Re-extract phasedata from the stored perwalk data
% Need to go back to the loco_analyzer to get raw phases
% For now, let's extract from the perwalk structure we built

% Actually, we need to re-run to get phasedata. Let me extract it during the main loop.
% For this plot, let's compute the phase values directly from all_walk_data's perwalk

% First, collect all raw phase data per condition
raw_phase_data = struct();

for c = 1:numel(all_walk_data)
    raw_phase_data(c).LM_phases = [];
    raw_phase_data(c).RM_phases = [];
    raw_phase_data(c).label = all_walk_data(c).label;
end

% We need to re-run the analysis to get phasedata, or modify the extraction loop
% For now, let's add phasedata extraction. Since we can't go back, let's re-run just the LED OFF conditions

fprintf('\nExtracting raw phase data for LM vs RM plot...\n');

for e = 1:numel(expdirs)
    expdir = expdirs{e};
    exp_label = exp_labels{e};

    % Re-initialize for this experiment
    trx = FBATrx('analysis_protocol', analysis_protocol, 'settingsdir', settingsdir, ...
        'datalocparamsfilestr', 'dataloc_params.txt');
    trx.AddExpDir(expdir, 'dooverwrite', false, 'openmovie', false);

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

    loco_analyzer = LimbBoutAnalyzer(trx, aptdata, tips_pos_body, legtip_landmarknums, ...
        groundcontact, digitalindicator, walking_scores, ...
        'phase_methods', {'phasediff_hilbert'});

    loco_analyzer.analyzeWalkAndStimConditions();

    % Extract for LED OFF only (for speed)
    for led = 1:numel(led_conditions)
        led_cond = led_conditions{led};
        led_label = led_labels{led};
        cond_label = sprintf('%s_%s', exp_label, led_label);

        % Find matching condition index
        cond_idx = find(strcmp({all_walk_data.label}, cond_label));
        if isempty(cond_idx), continue; end

        if ~loco_analyzer.walkMetrics.isKey(led_cond), continue; end
        walk_metrics = loco_analyzer.walkMetrics(led_cond);
        perwalk = walk_metrics.perwalk;

        LM_all = [];
        RM_all = [];

        for w = 1:numel(perwalk)
            if ~isfield(perwalk(w), 'phasediff_hilbert'), continue; end
            if ~isfield(perwalk(w).phasediff_hilbert, 'phasedata'), continue; end

            phasedata = perwalk(w).phasediff_hilbert.phasedata;
            % phasedata is 6 x nframes, legs: RF, RM, RH, LH, LM, LF
            LM_phase = phasedata(5, :);  % LM = leg 5
            RM_phase = phasedata(2, :);  % RM = leg 2

            % Only keep valid (non-NaN) paired data
            valid = ~isnan(LM_phase) & ~isnan(RM_phase);
            LM_all = [LM_all, LM_phase(valid)];
            RM_all = [RM_all, RM_phase(valid)];
        end

        raw_phase_data(cond_idx).LM_phases = LM_all;
        raw_phase_data(cond_idx).RM_phases = RM_all;
        fprintf('  %s: %d phase frames\n', cond_label, numel(LM_all));
    end
end

% Plot raw phases: LM vs RM scatter
for c = 1:numel(raw_phase_data)
    if isempty(raw_phase_data(c).LM_phases), continue; end

    nexttile;
    % Subsample for plotting if too many points
    n = numel(raw_phase_data(c).LM_phases);
    if n > 5000
        idx = randperm(n, 5000);
    else
        idx = 1:n;
    end

    scatter(raw_phase_data(c).LM_phases(idx), raw_phase_data(c).RM_phases(idx), ...
        3, cond_colors(raw_phase_data(c).label), 'filled', 'MarkerFaceAlpha', 0.3);
    hold on;
    plot([-pi, pi], [-pi, pi], 'k--', 'LineWidth', 1);

    xlabel('LM phase (rad)');
    ylabel('RM phase (rad)');
    xlim([-pi, pi]);
    ylim([-pi, pi]);
    axis square;
    title(raw_phase_data(c).label, 'Interpreter', 'none');

    % Show sample size
    text(0.05, 0.95, sprintf('n=%d', n), 'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 9);
end

% Plot raw phase distributions for LM and RM
for c = 1:numel(raw_phase_data)
    if isempty(raw_phase_data(c).LM_phases), continue; end

    nexttile;
    hold on;
    histogram(raw_phase_data(c).LM_phases, bin_edges, 'DisplayStyle', 'stairs', 'LineWidth', 2, ...
        'EdgeColor', 'b', 'DisplayName', 'LM (leg 5)');
    histogram(raw_phase_data(c).RM_phases, bin_edges, 'DisplayStyle', 'stairs', 'LineWidth', 2, ...
        'EdgeColor', 'r', 'DisplayName', 'RM (leg 2)');
    xlim([-pi, pi]);
    xlabel('Phase (rad)');
    ylabel('Count');
    title([raw_phase_data(c).label ' raw phases'], 'Interpreter', 'none');
    legend('Location', 'best');
end

sgtitle('Raw Phase: LM (leg 5) vs RM (leg 2)');
exportgraphics(gcf, fullfile(outputdir, 'plot5_raw_phase_LM_vs_RM.png'), 'Resolution', 150);

%% PLOT 6: Individual walk bout phase means - check for within-condition bimodality
% Shows each walk bout as a point to see if bimodality exists within conditions

figure('Position', [100, 100, 1200, 600]);
tiledlayout(2, 2, 'TileSpacing', 'compact');

% For ipsi_post_2 (combined LM_LH and RM_RH)
nexttile([1, 2]);
hold on;
x_offset = 0;
xtick_positions = [];
xtick_labels = {};

for c = 1:numel(all_walk_data)
    if isempty(all_walk_data(c).data), continue; end

    data = all_walk_data(c).data.phase_ipsi_post_2_mean;
    data = data(~isnan(data));
    n = numel(data);

    % Jitter x positions
    x = x_offset + 0.3 * (rand(n, 1) - 0.5);

    scatter(x, data, 15, cond_colors(all_walk_data(c).label), 'filled', 'MarkerFaceAlpha', 0.5);

    % Add mean line
    m = circ_mean(data);
    plot([x_offset-0.3, x_offset+0.3], [m, m], 'k-', 'LineWidth', 2);

    xtick_positions(end+1) = x_offset;
    xtick_labels{end+1} = all_walk_data(c).label;
    x_offset = x_offset + 1;
end

ylabel('Phase diff ipsi_post_2 (rad)', 'Interpreter', 'none');
ylim([-pi, pi]);
set(gca, 'XTick', xtick_positions, 'XTickLabel', xtick_labels, 'XTickLabelRotation', 45, 'TickLabelInterpreter', 'none');
title('Individual walk bout means - ipsi_post_2 (each dot = one walk)', 'Interpreter', 'none');
grid on;

% Separate histograms for each condition with more bins
for c = 1:numel(all_walk_data)
    if isempty(all_walk_data(c).data), continue; end
    if c > 2, break; end  % Only show first 2 (Dec5 LEDoff and LEDon or Dec5 and Dec6 LEDoff)

    nexttile;
    data = all_walk_data(c).data.phase_ipsi_post_2_mean;
    data = data(~isnan(data));

    histogram(data, 40, 'FaceColor', cond_colors(all_walk_data(c).label), 'FaceAlpha', 0.7);
    xlim([-pi, pi]);
    xlabel('Phase diff (rad)');
    ylabel('Count');
    title(sprintf('%s (n=%d walks)', all_walk_data(c).label, numel(data)), 'Interpreter', 'none');

    % Add vertical line at mean
    hold on;
    xline(circ_mean(data), 'k--', 'LineWidth', 2);

    % Check for bimodality hint: show std
    text(0.95, 0.95, sprintf('mean=%.2f\nstd=%.2f', circ_mean(data), circ_std(data)), ...
        'Units', 'normalized', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
end

sgtitle('Within-condition distribution of walk bout phase means', 'Interpreter', 'none');
exportgraphics(gcf, fullfile(outputdir, 'plot6_walkbout_means_bimodality_check.png'), 'Resolution', 150);

%% Summary statistics table
fprintf('\n=== Summary Statistics by Condition ===\n');
fprintf('%-25s %8s %8s %8s %8s %8s\n', 'Condition', 'N_walks', 'ipsi_p2', 'std', 'speed', 'L-R_diff');
fprintf('%s\n', repmat('-', 1, 70));

for c = 1:numel(all_walk_data)
    if isempty(all_walk_data(c).data), continue; end

    d = all_walk_data(c).data;
    nw = all_walk_data(c).nwalks;

    phase_valid = d.phase_ipsi_post_2_mean(~isnan(d.phase_ipsi_post_2_mean));
    phase_mean = circ_mean(phase_valid);
    phase_std = circ_std(phase_valid);

    speed_mean = mean(d.speed_mean, 'omitnan');

    L_mean = circ_mean(d.phase_LM_LH_mean(~isnan(d.phase_LM_LH_mean)));
    R_mean = circ_mean(d.phase_RM_RH_mean(~isnan(d.phase_RM_RH_mean)));
    LR_diff = circ_dist(L_mean, R_mean);

    fprintf('%-25s %8d %8.4f %8.4f %8.2f %8.4f\n', ...
        all_walk_data(c).label, nw, phase_mean, phase_std, speed_mean, LR_diff);
end

%% Direct Dec5 vs Dec6 comparison (LED OFF only)
fprintf('\n=== Dec 5 vs Dec 6 Comparison (LED OFF) ===\n');
dec5_idx = find(strcmp({all_walk_data.label}, 'Dec5_RigA_LEDoff'));
dec6_idx = find(strcmp({all_walk_data.label}, 'Dec6_RigC_LEDoff'));

if ~isempty(dec5_idx) && ~isempty(dec6_idx)
    d5 = all_walk_data(dec5_idx).data;
    d6 = all_walk_data(dec6_idx).data;

    phase5 = d5.phase_ipsi_post_2_mean(~isnan(d5.phase_ipsi_post_2_mean));
    phase6 = d6.phase_ipsi_post_2_mean(~isnan(d6.phase_ipsi_post_2_mean));

    fprintf('Dec5 ipsi_post_2: circ_mean=%.4f, circ_std=%.4f, n=%d\n', ...
        circ_mean(phase5), circ_std(phase5), numel(phase5));
    fprintf('Dec6 ipsi_post_2: circ_mean=%.4f, circ_std=%.4f, n=%d\n', ...
        circ_mean(phase6), circ_std(phase6), numel(phase6));
    fprintf('Difference in means: %.4f rad\n', circ_dist(circ_mean(phase5), circ_mean(phase6)));

    % Speed comparison
    speed5 = mean(d5.speed_mean, 'omitnan');
    speed6 = mean(d6.speed_mean, 'omitnan');
    fprintf('\nDec5 speed: %.2f mm/s\n', speed5);
    fprintf('Dec6 speed: %.2f mm/s\n', speed6);
    fprintf('Speed difference: %.2f mm/s\n', speed6 - speed5);
end

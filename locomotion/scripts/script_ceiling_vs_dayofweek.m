%% script_ceiling_vs_dayofweek.m
% Analyze ceiling walking patterns across 400 control experiments.
% Investigates: day-in-batch, position in session, experimenter,
% and residual floor-only AEPy variability after removing ceiling walks.

%% Setup
rootdir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20260225_VNC23controldata4onceilingwalkexplore';
savedir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Locomotion_analysis/claudeplots_20260225_ceiling_dayofweek';
if ~exist(savedir, 'dir'), mkdir(savedir); end

dd = dir(fullfile(rootdir, 'VNC*'));
dd = dd([dd.isdir]);
expnames = {dd.name};
nexp = numel(expnames);

fprintf('Found %d experiments\n', nexp);

%% Collect per-experiment ceiling stats
exp_date = NaT(nexp, 1);
exp_dayofweek = zeros(nexp, 1);  % 1=Sun, 2=Mon, ..., 7=Sat
exp_dayname = cell(nexp, 1);
exp_nwalks = zeros(nexp, 1);
exp_nwalks_ceiling = zeros(nexp, 1);
exp_nwalks_floor = zeros(nexp, 1);
exp_frac_ceiling = nan(nexp, 1);
exp_name = cell(nexp, 1);

ceiling_threshold = 0.5;  % fraction of walk frames classified as ceiling

for i = 1:nexp
    if mod(i, 50) == 0
        fprintf('Processing %d / %d\n', i, nexp);
    end

    expdir = fullfile(rootdir, expnames{i});
    exp_name{i} = expnames{i};

    % Extract date from directory name
    tokens = regexp(expnames{i}, '(\d{8})T', 'tokens');
    datestr_val = tokens{1}{1};
    exp_date(i) = datetime(datestr_val, 'InputFormat', 'yyyyMMdd');
    exp_dayofweek(i) = weekday(exp_date(i));  % 1=Sun, ..., 7=Sat
    exp_dayname{i} = datestr(exp_date(i), 'ddd');  % Mon, Tue, etc.

    % Load ceiling scores
    S = load(fullfile(expdir, 'scores_onceiling.mat'), 'allScores');
    allScores = S.allScores;

    % Load walk struct
    W = load(fullfile(expdir, 'locomotion_walkstruct.mat'), 'walk_struct_OFF');
    ws = W.walk_struct_OFF;

    nwalks = numel(ws.fly);
    exp_nwalks(i) = nwalks;

    if nwalks == 0
        exp_frac_ceiling(i) = NaN;
        continue;
    end

    % Classify each walk as ceiling or floor
    is_ceiling = false(1, nwalks);
    for w = 1:nwalks
        fly_id = ws.fly(w);
        t0 = ws.walk_t0(w);
        t1 = ws.walk_t1(w);

        % Use postprocessed binary labels
        pp = allScores.postprocessed{fly_id};
        walk_frames = t0:t1;

        % Handle frame index bounds
        valid = walk_frames >= 1 & walk_frames <= numel(pp);
        if sum(valid) == 0
            continue;
        end

        frac_on_ceiling = mean(pp(walk_frames(valid)));
        is_ceiling(w) = frac_on_ceiling > ceiling_threshold;
    end

    exp_nwalks_ceiling(i) = sum(is_ceiling);
    exp_nwalks_floor(i) = sum(~is_ceiling);
    exp_frac_ceiling(i) = mean(is_ceiling);
end

fprintf('Done collecting. Total walks: %d, ceiling: %d, floor: %d\n', ...
    sum(exp_nwalks), sum(exp_nwalks_ceiling), sum(exp_nwalks_floor));

%% Compute day-in-batch (1st, 2nd, 3rd day of experiments that week)
daynames_order = {'Mon','Tue','Wed','Thu','Fri','Sat','Sun'};
daynum_order = [2 3 4 5 6 7 1];

dow_plot = zeros(nexp, 1);
for k = 1:7
    dow_plot(exp_dayofweek == daynum_order(k)) = k;
end

exp_weeknum = week(exp_date, 'weekofyear');
exp_year = year(exp_date);
exp_weekid = exp_year * 100 + exp_weeknum;

unique_weeks = unique(exp_weekid);
exp_day_in_batch = zeros(nexp, 1);

for i = 1:numel(unique_weeks)
    wk = unique_weeks(i);
    idx = find(exp_weekid == wk);
    dates_in_week = unique(exp_date(idx));
    dates_in_week = sort(dates_in_week);
    for j = 1:numel(dates_in_week)
        day_idx = idx(exp_date(idx) == dates_in_week(j));
        exp_day_in_batch(day_idx) = j;
    end
end

%% Compute position in daily session using full experiment list
S_full = load('/groups/branson/bransonlab/flydisco_linelevel_VNC/CollectedVNC23PerFrameStats20260205.mat', 'expdirs');
all_expdirs = S_full.expdirs;

exp_datetime = NaT(nexp, 1);
for i = 1:nexp
    tokens = regexp(expnames{i}, '(\d{8}T\d{6})', 'tokens');
    exp_datetime(i) = datetime(tokens{1}{1}, 'InputFormat', 'yyyyMMdd''T''HHmmss');
end

all_datetime = NaT(numel(all_expdirs), 1);
for i = 1:numel(all_expdirs)
    [~, ename] = fileparts(all_expdirs{i});
    tokens = regexp(ename, '(\d{8}T\d{6})', 'tokens');
    if ~isempty(tokens)
        all_datetime(i) = datetime(tokens{1}{1}, 'InputFormat', 'yyyyMMdd''T''HHmmss');
    end
end

all_dateonly = dateshift(all_datetime, 'start', 'day');
ctrl_dateonly = dateshift(exp_datetime, 'start', 'day');

exp_rank_in_day = nan(nexp, 1);
exp_nexp_that_day = zeros(nexp, 1);
exp_frac_position = nan(nexp, 1);

unique_ctrl_dates = unique(ctrl_dateonly);
for d = 1:numel(unique_ctrl_dates)
    this_date = unique_ctrl_dates(d);
    all_idx = find(all_dateonly == this_date);
    all_times_today = all_datetime(all_idx);
    [~, sort_order] = sort(all_times_today);
    n_total_today = numel(all_idx);

    rank_map = containers.Map('KeyType', 'char', 'ValueType', 'double');
    for k = 1:n_total_today
        [~, ename] = fileparts(all_expdirs{all_idx(sort_order(k))});
        rank_map(ename) = k;
    end

    ctrl_idx = find(ctrl_dateonly == this_date);
    for j = 1:numel(ctrl_idx)
        ci = ctrl_idx(j);
        exp_nexp_that_day(ci) = n_total_today;
        if rank_map.isKey(expnames{ci})
            exp_rank_in_day(ci) = rank_map(expnames{ci});
            exp_frac_position(ci) = (rank_map(expnames{ci}) - 1) / max(n_total_today - 1, 1);
        end
    end
end

%% Extract metadata from Metadata.xml
exp_experimenter = cell(nexp, 1);
exp_rig = cell(nexp, 1);
exp_temperature = nan(nexp, 1);
exp_humidity = nan(nexp, 1);

for i = 1:nexp
    xmlfile = fullfile(rootdir, expnames{i}, 'Metadata.xml');
    if ~exist(xmlfile, 'file'), continue; end
    txt = fileread(xmlfile);

    tok = regexp(txt, 'experimenter="([^"]*)"', 'tokens');
    if ~isempty(tok), exp_experimenter{i} = tok{1}{1}; end

    tok = regexp(txt, 'rig="([^"]*)"', 'tokens');
    if ~isempty(tok), exp_rig{i} = tok{1}{1}; end

    tok = regexp(txt, 'temperature="([^"]*)"', 'tokens');
    if ~isempty(tok), exp_temperature(i) = str2double(tok{1}{1}); end

    tok = regexp(txt, 'humidity="([^"]*)"', 'tokens');
    if ~isempty(tok), exp_humidity(i) = str2double(tok{1}{1}); end
end

%% Plot 1: Ceiling fraction vs day-in-batch (boxplot)
present_days_batch = unique(exp_day_in_batch);
present_days_batch = present_days_batch(present_days_batch > 0);
batch_labels = arrayfun(@(x) sprintf('Day %d', x), present_days_batch, 'uni', false);

figure('Position', [100 100 900 500]);
boxplot(exp_frac_ceiling, exp_day_in_batch, 'Labels', batch_labels, 'Widths', 0.5);
hold on;
rng(42);
cmap = [0.2 0.4 0.8; 0.8 0.4 0.2; 0.4 0.7 0.3];
for d = present_days_batch(:)'
    idx = exp_day_in_batch == d;
    n = sum(idx);
    jitter = 0.2 * (rand(n, 1) - 0.5);
    plot(find(present_days_batch == d) + jitter, exp_frac_ceiling(idx), '.', ...
        'MarkerSize', 10, 'Color', cmap(d,:));
end
hold off;
ylabel('Fraction of walks on ceiling', 'Interpreter', 'none');
xlabel('Day in weekly batch', 'Interpreter', 'none');
title('Ceiling walk fraction vs day-in-batch (400 control experiments)', 'Interpreter', 'none');
set(gca, 'TickLabelInterpreter', 'none');
grid on;
for d = present_days_batch(:)'
    n = sum(exp_day_in_batch == d);
    text(find(present_days_batch == d), -0.08, sprintf('n=%d', n), ...
        'HorizontalAlignment', 'center', 'FontSize', 9);
end
saveas(gcf, fullfile(savedir, 'ceiling_frac_vs_day_in_batch_boxplot.png'));

%% Plot 2: Ceiling fraction over time, colored by day-in-batch
figure('Position', [100 100 1200 500]);
hold on;
for d = 1:2
    idx = exp_day_in_batch == d;
    scatter(exp_date(idx), exp_frac_ceiling(idx), 40, cmap(d,:), 'filled', 'MarkerFaceAlpha', 0.6);
end
hold off;
legend({'Day 1 of batch', 'Day 2 of batch'}, 'Location', 'best');
ylabel('Fraction of walks on ceiling', 'Interpreter', 'none');
xlabel('Date', 'Interpreter', 'none');
title('Ceiling walk fraction over time, colored by day-in-batch', 'Interpreter', 'none');
grid on; ylim([-0.05 1.05]);
saveas(gcf, fullfile(savedir, 'ceiling_frac_over_time_by_batch.png'));

%% Plot 3: Ceiling fraction vs position in daily session
figure('Position', [100 100 1000 500]);
hold on;
for d = 1:2
    idx = exp_day_in_batch == d;
    scatter(exp_frac_position(idx), exp_frac_ceiling(idx), 40, cmap(d,:), 'filled', 'MarkerFaceAlpha', 0.5);
end
hold off;
legend({'Day 1 of batch', 'Day 2 of batch'}, 'Location', 'best');
xlabel('Fractional position in day''s session (0=first, 1=last)', 'Interpreter', 'none');
ylabel('Fraction of walks on ceiling', 'Interpreter', 'none');
title('Ceiling walk fraction vs position in daily session', 'Interpreter', 'none');
grid on; ylim([-0.05 1.05]);
saveas(gcf, fullfile(savedir, 'ceiling_frac_vs_position_in_day.png'));

%% Plot 4: Ceiling fraction over time, colored by experimenter
figure('Position', [100 100 1200 500]);
unique_exp = unique(exp_experimenter(~cellfun(@isempty, exp_experimenter)));
cmap_exp = [0.8 0.2 0.2; 0.2 0.5 0.8; 0.5 0.7 0.2];
hold on;
for j = 1:numel(unique_exp)
    idx = strcmp(exp_experimenter, unique_exp{j});
    scatter(exp_date(idx), exp_frac_ceiling(idx), 40, cmap_exp(j,:), 'filled', 'MarkerFaceAlpha', 0.6);
end
hold off;
legend(unique_exp, 'Location', 'best', 'Interpreter', 'none');
ylabel('Fraction of walks on ceiling', 'Interpreter', 'none');
xlabel('Date', 'Interpreter', 'none');
title('Ceiling walk fraction over time, colored by experimenter', 'Interpreter', 'none');
grid on; ylim([-0.05 1.05]);
saveas(gcf, fullfile(savedir, 'ceiling_frac_over_time_by_experimenter.png'));

%% Collect per-walk AEPy and ceiling labels across all experiments
max_walks = 500000;
walk_AEPy_RM = nan(max_walks, 1);
walk_AEPy_LM = nan(max_walks, 1);
walk_AEPy_mid = nan(max_walks, 1);
walk_is_ceiling = false(max_walks, 1);
walk_expidx = zeros(max_walks, 1);
walk_velmag = nan(max_walks, 1);

total_walks = 0;
for i = 1:nexp
    if mod(i, 50) == 0
        fprintf('Collecting walks %d / %d\n', i, nexp);
    end
    expdir = fullfile(rootdir, expnames{i});
    W = load(fullfile(expdir, 'locomotion_walkstruct.mat'), 'walk_struct_OFF');
    ws = W.walk_struct_OFF;
    Sc = load(fullfile(expdir, 'scores_onceiling.mat'), 'allScores');
    allScores = Sc.allScores;

    nw = numel(ws.fly);
    if nw == 0, continue; end
    idx = total_walks + (1:nw);

    walk_AEPy_RM(idx) = ws.step_AEPy_RM;
    walk_AEPy_LM(idx) = ws.step_AEPy_LM;
    walk_AEPy_mid(idx) = (ws.step_AEPy_RM + ws.step_AEPy_LM) / 2;
    walk_velmag(idx) = ws.velmag_ctr;
    walk_expidx(idx) = i;

    for w = 1:nw
        fly_id = ws.fly(w);
        t0 = ws.walk_t0(w);
        t1 = ws.walk_t1(w);
        pp = allScores.postprocessed{fly_id};
        frames = t0:min(t1, numel(pp));
        if isempty(frames), continue; end
        walk_is_ceiling(total_walks + w) = mean(pp(frames)) > ceiling_threshold;
    end
    total_walks = total_walks + nw;
end

walk_AEPy_RM = walk_AEPy_RM(1:total_walks);
walk_AEPy_LM = walk_AEPy_LM(1:total_walks);
walk_AEPy_mid = walk_AEPy_mid(1:total_walks);
walk_is_ceiling = walk_is_ceiling(1:total_walks);
walk_expidx = walk_expidx(1:total_walks);
walk_velmag = walk_velmag(1:total_walks);

%% Compute per-experiment floor-only stats
exp_AEPy_floor_mean = nan(nexp, 1);
exp_AEPy_floor_std = nan(nexp, 1);
exp_AEPy_floor_n = zeros(nexp, 1);
exp_AEPy_all_mean = nan(nexp, 1);
exp_velmag_floor_mean = nan(nexp, 1);

for i = 1:nexp
    idx_all = walk_expidx == i;
    idx_floor = idx_all & ~walk_is_ceiling;
    if sum(idx_floor) > 0
        exp_AEPy_floor_mean(i) = mean(walk_AEPy_mid(idx_floor), 'omitnan');
        exp_AEPy_floor_std(i) = std(walk_AEPy_mid(idx_floor), 'omitnan');
        exp_AEPy_floor_n(i) = sum(idx_floor & ~isnan(walk_AEPy_mid));
        exp_velmag_floor_mean(i) = mean(walk_velmag(idx_floor), 'omitnan');
    end
    if sum(idx_all) > 0
        exp_AEPy_all_mean(i) = mean(walk_AEPy_mid(idx_all), 'omitnan');
    end
end

%% Plot 5: AEPy over time — all walks vs floor-only
figure('Position', [100 100 1400 800]);
ylims = [-6 4];

subplot(2,1,1);
has_data = ~isnan(exp_AEPy_all_mean);
errorbar(exp_date(has_data), exp_AEPy_all_mean(has_data), exp_AEPy_floor_std(has_data), ...
    'o', 'MarkerSize', 5, 'Color', [0.5 0.5 0.5], 'MarkerFaceColor', [0.5 0.5 0.5]);
ylabel('AEPy (mid legs)', 'Interpreter', 'none');
title('ALL walks - mean AEPy per experiment over time', 'Interpreter', 'none');
ylim(ylims); grid on;

subplot(2,1,2);
has_data = ~isnan(exp_AEPy_floor_mean);
errorbar(exp_date(has_data), exp_AEPy_floor_mean(has_data), exp_AEPy_floor_std(has_data), ...
    'o', 'MarkerSize', 5, 'Color', [0.2 0.5 0.8], 'MarkerFaceColor', [0.2 0.5 0.8]);
ylabel('AEPy (mid legs)', 'Interpreter', 'none');
xlabel('Date', 'Interpreter', 'none');
title('FLOOR-ONLY walks - mean AEPy per experiment over time', 'Interpreter', 'none');
ylim(ylims); grid on;
saveas(gcf, fullfile(savedir, 'AEPy_over_time_all_vs_floor.png'));

%% Plot 6: AEPy vs speed analysis
figure('Position', [100 100 1200 900]);
has_floor = ~isnan(exp_AEPy_floor_mean);
unique_expers = unique(exp_experimenter(~cellfun(@isempty, exp_experimenter)));

subplot(2,2,1);
hold on;
for j = 1:numel(unique_expers)
    idx = strcmp(exp_experimenter, unique_expers{j}) & has_floor;
    scatter(exp_velmag_floor_mean(idx), exp_AEPy_floor_mean(idx), 30, cmap_exp(j,:), 'filled', 'MarkerFaceAlpha', 0.6);
end
hold off;
legend(unique_expers, 'Location', 'best', 'Interpreter', 'none');
xlabel('Mean speed (floor walks)', 'Interpreter', 'none');
ylabel('Mean AEPy (floor walks)', 'Interpreter', 'none');
title('Floor-only AEPy vs speed, by experimenter', 'Interpreter', 'none');
grid on;

subplot(2,2,2);
hold on;
for j = 1:numel(unique_expers)
    idx = strcmp(exp_experimenter, unique_expers{j}) & has_floor;
    scatter(exp_date(idx), exp_velmag_floor_mean(idx), 30, cmap_exp(j,:), 'filled', 'MarkerFaceAlpha', 0.6);
end
hold off;
legend(unique_expers, 'Location', 'best', 'Interpreter', 'none');
xlabel('Date', 'Interpreter', 'none');
ylabel('Mean speed (floor walks)', 'Interpreter', 'none');
title('Floor-only speed over time', 'Interpreter', 'none');
grid on;

subplot(2,2,3);
valid = has_floor & ~isnan(exp_velmag_floor_mean);
p_fit = polyfit(exp_velmag_floor_mean(valid), exp_AEPy_floor_mean(valid), 1);
AEPy_residual = nan(nexp, 1);
AEPy_residual(valid) = exp_AEPy_floor_mean(valid) - polyval(p_fit, exp_velmag_floor_mean(valid));
hold on;
for j = 1:numel(unique_expers)
    idx = strcmp(exp_experimenter, unique_expers{j}) & valid;
    scatter(exp_date(idx), AEPy_residual(idx), 30, cmap_exp(j,:), 'filled', 'MarkerFaceAlpha', 0.6);
end
hold off;
legend(unique_expers, 'Location', 'best', 'Interpreter', 'none');
xlabel('Date', 'Interpreter', 'none');
ylabel('AEPy residual (speed-corrected)', 'Interpreter', 'none');
title('Floor-only AEPy after removing speed effect', 'Interpreter', 'none');
grid on;

subplot(2,2,4);
text(0.1, 0.9, sprintf('Original AEPy std: %.3f', std(exp_AEPy_floor_mean(valid))), 'FontSize', 12, 'Units', 'normalized');
text(0.1, 0.75, sprintf('Residual AEPy std: %.3f', std(AEPy_residual(valid))), 'FontSize', 12, 'Units', 'normalized');
text(0.1, 0.6, sprintf('Variance explained by speed: %.1f%%', ...
    100*(1 - var(AEPy_residual(valid))/var(exp_AEPy_floor_mean(valid)))), 'FontSize', 12, 'Units', 'normalized');
text(0.1, 0.4, sprintf('Speed slope: %.3f AEPy/mm/s', p_fit(1)), 'FontSize', 12, 'Units', 'normalized');
axis off;
title('Summary stats', 'Interpreter', 'none');
saveas(gcf, fullfile(savedir, 'AEPy_floor_vs_speed_analysis.png'));

%% Print summary tables
fprintf('\n--- Summary by day-in-batch ---\n');
fprintf('%-8s  %5s  %8s  %8s  %8s  %10s\n', 'Batch', 'Nexp', 'Nwalks', 'Ceiling', 'Floor', 'FracCeil');
for d = present_days_batch(:)'
    idx = exp_day_in_batch == d;
    fprintf('%-8s  %5d  %8d  %8d  %8d  %10.3f\n', ...
        sprintf('Day %d', d), sum(idx), sum(exp_nwalks(idx)), ...
        sum(exp_nwalks_ceiling(idx)), sum(exp_nwalks_floor(idx)), ...
        sum(exp_nwalks_ceiling(idx)) / sum(exp_nwalks(idx)));
end

fprintf('\n--- Ceiling fraction by experimenter ---\n');
fprintf('%-20s  %5s  %10s  %10s\n', 'Experimenter', 'Nexp', 'MedFrac', 'Frac>0.5');
for j = 1:numel(unique_expers)
    idx = strcmp(exp_experimenter, unique_expers{j});
    fprintf('%-20s  %5d  %10.3f  %10.3f\n', unique_expers{j}, sum(idx), ...
        median(exp_frac_ceiling(idx)), mean(exp_frac_ceiling(idx) > 0.5));
end

fprintf('\n--- Correlations with floor-only AEPy ---\n');
valid = has_floor & ~isnan(exp_velmag_floor_mean);
[r,p] = corr(exp_velmag_floor_mean(valid), exp_AEPy_floor_mean(valid));
fprintf('Speed:       r=%.3f, p=%.2e\n', r, p);
[r,p] = corr(exp_frac_position(valid), exp_AEPy_floor_mean(valid));
fprintf('Position:    r=%.3f, p=%.2e\n', r, p);
[r,p] = corr(exp_temperature(valid), exp_AEPy_floor_mean(valid));
fprintf('Temperature: r=%.3f, p=%.2e\n', r, p);

fprintf('\nAEPy variability: all walks std=%.2f, floor-only std=%.2f\n', ...
    std(exp_AEPy_all_mean, 'omitnan'), std(exp_AEPy_floor_mean, 'omitnan'));

%% Save workspace
save(fullfile(savedir, 'ceiling_dayofweek_data.mat'), ...
    'exp_date', 'exp_dayofweek', 'exp_dayname', 'exp_nwalks', ...
    'exp_nwalks_ceiling', 'exp_nwalks_floor', 'exp_frac_ceiling', ...
    'exp_name', 'dow_plot', 'daynames_order', 'exp_day_in_batch', ...
    'exp_frac_position', 'exp_rank_in_day', 'exp_experimenter', 'exp_rig', ...
    'exp_temperature', 'exp_humidity', 'exp_AEPy_floor_mean', ...
    'exp_AEPy_all_mean', 'exp_velmag_floor_mean');
fprintf('Saved workspace to %s\n', savedir);

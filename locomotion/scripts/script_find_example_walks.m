%% Find example walks for visual inspection
% Find walks that exemplify the phase-movement relationship for movie review
%
% Key trends from analysis:
%   - Low phase (more negative) + high speed + high turning
%   - High phase (less negative) + low speed + low turning
%
% Uses composite score: (zscore(velmag_ctr) + zscore(absdtheta)) / 2

modpath

%% Load per-walk data
rootdir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20260211_exploringphaseissue';
outputdir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Locomotion_analysis/claudeplots_20260211';
cache_file = fullfile(rootdir, 'all_walk_data_cache.mat');

if ~exist(cache_file, 'file')
    error('Run script_phase_correlation_perwalk.m first to generate cache');
end

load(cache_file, 'all_walk_data', 'walk_counter');
fprintf('Loaded %d walks\n', walk_counter);

%% Extract key variables
phase = [all_walk_data.phase_ipsi_post_2]';
speed = [all_walk_data.velmag_ctr]';
turning = [all_walk_data.absdtheta]';
walk_dur = [all_walk_data.walk_duration]';

% Filter: valid data and minimum walk duration (at least 30 frames = 1 sec at 30fps)
min_duration = 30;
valid = ~isnan(phase) & ~isnan(speed) & ~isnan(turning) & walk_dur >= min_duration;

fprintf('Valid walks with duration >= %d frames: %d\n', min_duration, sum(valid));

%% Compute composite score (zscore of speed + turning)
% Only compute on valid data
speed_z = nan(size(speed));
turning_z = nan(size(turning));
speed_z(valid) = zscore(speed(valid));
turning_z(valid) = zscore(turning(valid));
composite = (speed_z + turning_z) / 2;

fprintf('\nComposite score = (zscore(velmag_ctr) + zscore(absdtheta)) / 2\n');
fprintf('  Speed (velmag_ctr): mean=%.2f, std=%.2f\n', mean(speed(valid)), std(speed(valid)));
fprintf('  Turning (absdtheta): mean=%.2f, std=%.2f\n', mean(turning(valid)), std(turning(valid)));

%% Compute percentiles for selection
phase_valid = phase(valid);
composite_valid = composite(valid);

phase_10 = prctile(phase_valid, 10);  % low phase (more negative)
phase_90 = prctile(phase_valid, 90);  % high phase (less negative)
comp_25 = prctile(composite_valid, 25);  % low movement intensity
comp_75 = prctile(composite_valid, 75);  % high movement intensity

fprintf('\nPhase percentiles: 10th=%.3f, 90th=%.3f\n', phase_10, phase_90);
fprintf('Composite percentiles: 25th=%.2f, 75th=%.2f\n', comp_25, comp_75);

%% Find example walks

% Category 1: LOW phase + HIGH composite (expected by correlation)
cat1_idx = find(valid & phase < phase_10 & composite > comp_75);
fprintf('\n=== LOW phase + HIGH movement: %d walks ===\n', numel(cat1_idx));

% Category 2: HIGH phase + LOW composite (expected by correlation)
cat2_idx = find(valid & phase > phase_90 & composite < comp_25);
fprintf('=== HIGH phase + LOW movement: %d walks ===\n', numel(cat2_idx));

% Category 3: LOW phase + LOW composite (counter-trend)
cat3_idx = find(valid & phase < phase_10 & composite < comp_25);
fprintf('=== LOW phase + LOW movement (counter-trend): %d walks ===\n', numel(cat3_idx));

% Category 4: HIGH phase + HIGH composite (counter-trend)
cat4_idx = find(valid & phase > phase_90 & composite > comp_75);
fprintf('=== HIGH phase + HIGH movement (counter-trend): %d walks ===\n', numel(cat4_idx));

%% Print top examples for each category
n_examples = 5;

fprintf('\n\n======================================================\n');
fprintf('CATEGORY 1: LOW phase + HIGH movement (typical trend)\n');
fprintf('======================================================\n');
fprintf('Phase < %.2f rad, Composite > %.2f\n', phase_10, comp_75);
fprintf('These walks have more negative phase with fast walking and high turning.\n');
fprintf('Look for: vigorous, coordinated tripod gait\n\n');
if ~isempty(cat1_idx)
    % Sort by walk duration (prefer longer for easier viewing)
    [~, sort_idx] = sort(walk_dur(cat1_idx), 'descend');
    print_walk_details(all_walk_data, phase, speed, turning, composite, walk_dur, cat1_idx(sort_idx), n_examples);
end

fprintf('\n\n======================================================\n');
fprintf('CATEGORY 2: HIGH phase + LOW movement (typical trend)\n');
fprintf('======================================================\n');
fprintf('Phase > %.2f rad, Composite < %.2f\n', phase_90, comp_25);
fprintf('These walks have less negative phase with slow walking and low turning.\n');
fprintf('Look for: sluggish, possibly altered coordination\n\n');
if ~isempty(cat2_idx)
    [~, sort_idx] = sort(walk_dur(cat2_idx), 'descend');
    print_walk_details(all_walk_data, phase, speed, turning, composite, walk_dur, cat2_idx(sort_idx), n_examples);
end

fprintf('\n\n======================================================\n');
fprintf('CATEGORY 3: LOW phase + LOW movement (counter-trend)\n');
fprintf('======================================================\n');
fprintf('Phase < %.2f rad, Composite < %.2f\n', phase_10, comp_25);
fprintf('These walks have good phase coordination but slow movement.\n');
fprintf('Interesting exceptions - why slow despite good coordination?\n\n');
if ~isempty(cat3_idx)
    [~, sort_idx] = sort(walk_dur(cat3_idx), 'descend');
    print_walk_details(all_walk_data, phase, speed, turning, composite, walk_dur, cat3_idx(sort_idx), n_examples);
end

fprintf('\n\n======================================================\n');
fprintf('CATEGORY 4: HIGH phase + HIGH movement (counter-trend)\n');
fprintf('======================================================\n');
fprintf('Phase > %.2f rad, Composite > %.2f\n', phase_90, comp_75);
fprintf('These walks have poor phase coordination but fast movement.\n');
fprintf('Interesting exceptions - fast despite poor coordination?\n\n');
if ~isempty(cat4_idx)
    [~, sort_idx] = sort(walk_dur(cat4_idx), 'descend');
    print_walk_details(all_walk_data, phase, speed, turning, composite, walk_dur, cat4_idx(sort_idx), n_examples);
end

%% Summary table for easy reference
fprintf('\n\n======================================================\n');
fprintf('SUMMARY: Quick reference for movie review\n');
fprintf('======================================================\n\n');

fprintf('Experiment directories:\n');
fprintf('  %s\n\n', rootdir);

categories = {'LOW_phase_HIGH_movement', 'HIGH_phase_LOW_movement', ...
              'LOW_phase_LOW_movement', 'HIGH_phase_HIGH_movement'};
cat_indices = {cat1_idx, cat2_idx, cat3_idx, cat4_idx};
cat_descriptions = {'Typical: good coord + vigorous', 'Typical: poor coord + sluggish', ...
                    'Exception: good coord + sluggish', 'Exception: poor coord + vigorous'};

for c = 1:4
    fprintf('--- %s (%s) ---\n', categories{c}, cat_descriptions{c});
    idx = cat_indices{c};
    if isempty(idx)
        fprintf('  No examples found\n\n');
        continue;
    end

    % Sort by walk duration (prefer longer)
    [~, sort_idx] = sort(walk_dur(idx), 'descend');
    idx = idx(sort_idx);

    for i = 1:min(3, numel(idx))
        w = idx(i);
        fprintf('  %s\n', all_walk_data(w).expname);
        fprintf('    Fly %d, frames %d-%d (duration: %d frames)\n', ...
            all_walk_data(w).fly_in_exp, ...
            all_walk_data(w).walk_t0, ...
            all_walk_data(w).walk_t1, ...
            walk_dur(w));
        fprintf('    Phase=%.2f, Speed=%.1f, Turning=%.2f\n\n', ...
            phase(w), speed(w), turning(w));
    end
end

%% Save examples to file
examples = struct();
examples.cat1_low_phase_high_movement = cat1_idx;
examples.cat2_high_phase_low_movement = cat2_idx;
examples.cat3_low_phase_low_movement = cat3_idx;
examples.cat4_high_phase_high_movement = cat4_idx;
examples.phase = phase;
examples.speed = speed;
examples.turning = turning;
examples.composite = composite;
examples.walk_dur = walk_dur;
examples.phase_10 = phase_10;
examples.phase_90 = phase_90;
examples.comp_25 = comp_25;
examples.comp_75 = comp_75;

save(fullfile(outputdir, 'example_walks_for_review.mat'), 'examples', 'all_walk_data');
fprintf('\nExamples saved to: %s\n', fullfile(outputdir, 'example_walks_for_review.mat'));

%% Scatter plot highlighting examples
figure('Position', [100, 100, 900, 700]);
hold on;

% All walks (gray background)
scatter(phase(valid), composite(valid), 8, [0.7, 0.7, 0.7], 'filled', 'MarkerFaceAlpha', 0.3);

% Highlight categories with larger markers
if ~isempty(cat1_idx)
    scatter(phase(cat1_idx), composite(cat1_idx), 40, 'b', 'filled', 'MarkerEdgeColor', 'k');
end
if ~isempty(cat2_idx)
    scatter(phase(cat2_idx), composite(cat2_idx), 40, 'r', 'filled', 'MarkerEdgeColor', 'k');
end
if ~isempty(cat3_idx)
    scatter(phase(cat3_idx), composite(cat3_idx), 40, 'c', 'filled', 'MarkerEdgeColor', 'k');
end
if ~isempty(cat4_idx)
    scatter(phase(cat4_idx), composite(cat4_idx), 40, 'm', 'filled', 'MarkerEdgeColor', 'k');
end

% Reference lines for thresholds
xline(phase_10, 'b--', 'LineWidth', 1.5, 'Label', '10th pctl', 'LabelHorizontalAlignment', 'left');
xline(phase_90, 'r--', 'LineWidth', 1.5, 'Label', '90th pctl', 'LabelHorizontalAlignment', 'right');
yline(comp_25, 'k:', 'LineWidth', 1.5, 'Label', '25th pctl');
yline(comp_75, 'k:', 'LineWidth', 1.5, 'Label', '75th pctl');

xlabel('Phase difference ipsi\_post\_2 (rad)', 'Interpreter', 'none');
ylabel('Composite score (zscore speed + zscore turning)', 'Interpreter', 'none');
title('Example walks for visual inspection', 'Interpreter', 'none');

% Zoom to data range
xlim([prctile(phase_valid, 1), prctile(phase_valid, 99)]);
ylim([prctile(composite_valid, 1), prctile(composite_valid, 99)]);

legend({'All walks', ...
        'Cat1: Low phase + High movement', ...
        'Cat2: High phase + Low movement', ...
        'Cat3: Low phase + Low movement', ...
        'Cat4: High phase + High movement'}, ...
        'Location', 'northeast');

exportgraphics(gcf, fullfile(outputdir, 'example_walks_scatter.png'), 'Resolution', 150);
fprintf('Scatter plot saved.\n');

%% Additional plot: Show the two components
figure('Position', [100, 100, 1200, 500]);

subplot(1,2,1);
hold on;
scatter(phase(valid), speed(valid), 8, [0.7, 0.7, 0.7], 'filled', 'MarkerFaceAlpha', 0.3);
if ~isempty(cat1_idx), scatter(phase(cat1_idx), speed(cat1_idx), 30, 'b', 'filled'); end
if ~isempty(cat2_idx), scatter(phase(cat2_idx), speed(cat2_idx), 30, 'r', 'filled'); end
if ~isempty(cat3_idx), scatter(phase(cat3_idx), speed(cat3_idx), 30, 'c', 'filled'); end
if ~isempty(cat4_idx), scatter(phase(cat4_idx), speed(cat4_idx), 30, 'm', 'filled'); end
xlabel('Phase (rad)', 'Interpreter', 'none');
ylabel('Speed (velmag\_ctr)', 'Interpreter', 'none');
title(sprintf('Phase vs Speed (r=%.2f)', corr(phase(valid), speed(valid))), 'Interpreter', 'none');
xlim([prctile(phase_valid, 1), prctile(phase_valid, 99)]);

subplot(1,2,2);
hold on;
scatter(phase(valid), turning(valid), 8, [0.7, 0.7, 0.7], 'filled', 'MarkerFaceAlpha', 0.3);
if ~isempty(cat1_idx), scatter(phase(cat1_idx), turning(cat1_idx), 30, 'b', 'filled'); end
if ~isempty(cat2_idx), scatter(phase(cat2_idx), turning(cat2_idx), 30, 'r', 'filled'); end
if ~isempty(cat3_idx), scatter(phase(cat3_idx), turning(cat3_idx), 30, 'c', 'filled'); end
if ~isempty(cat4_idx), scatter(phase(cat4_idx), turning(cat4_idx), 30, 'm', 'filled'); end
xlabel('Phase (rad)', 'Interpreter', 'none');
ylabel('Turning (absdtheta)', 'Interpreter', 'none');
title(sprintf('Phase vs Turning (r=%.2f)', corr(phase(valid), turning(valid))), 'Interpreter', 'none');
xlim([prctile(phase_valid, 1), prctile(phase_valid, 99)]);

sgtitle('Phase vs individual components', 'Interpreter', 'none');
exportgraphics(gcf, fullfile(outputdir, 'example_walks_components.png'), 'Resolution', 150);
fprintf('Components plot saved.\n');

fprintf('\nDone.\n');

%% Local function to print walk details
function print_walk_details(all_walk_data, phase, speed, turning, composite, walk_dur, indices, n_print)
    fprintf('%-45s %5s %8s %8s %8s %8s %8s %8s %8s\n', ...
        'Experiment', 'Fly', 'Start', 'End', 'Dur', 'Phase', 'Speed', 'Turn', 'Comp');
    fprintf('%s\n', repmat('-', 1, 115));

    for i = 1:min(n_print, numel(indices))
        w = indices(i);
        fprintf('%-45s %5d %8d %8d %8d %8.3f %8.2f %8.2f %8.2f\n', ...
            all_walk_data(w).expname, ...
            all_walk_data(w).fly_in_exp, ...
            all_walk_data(w).walk_t0, ...
            all_walk_data(w).walk_t1, ...
            walk_dur(w), ...
            phase(w), ...
            speed(w), ...
            turning(w), ...
            composite(w));
    end
end

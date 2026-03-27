%% script_make_fly_gifs_3windows
% Generate GIFs with 3 clips per fly at 25%, 50%, 75% of trajectory.
% Each GIF shows the fly at early, middle, and late timepoints with
% a brief pause (black frame) between segments.
%
% White border = onceiling classified. Each segment ~66 frames.
%
% Builds review set from damaged_fly_deviants.mat:
%   - All flies with <10% walking (damaged candidates + grey zone)
%   - 5% random normal flies (>10% walk) mixed in as controls
%   - Shuffled so reviewer is blinded

modpath;

outdir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Locomotion_analysis/claudeplots_20260319_onceiling_allVNC/review_gifs_3win_full/';
if ~isfolder(outdir), mkdir(outdir); end

%% Build clip list from data
load('/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/damaged_fly_deviants.mat', 'all_devs');

wk = [all_devs.pct_walk];
oc = [all_devs.pct_oc];

% All candidates with <10% walking
candidate_idx = find(wk < 10);

% Random 5% normal flies (>10% walk) as controls
normal_idx = find(wk >= 10);
rng(42);  % reproducible
n_random = max(1, round(0.05 * numel(candidate_idx)));
rand_idx = normal_idx(randperm(numel(normal_idx), n_random));

all_idx = [candidate_idx, rand_idx];

% Build labels: encode sex, oc%, walk% but NOT whether it's a candidate or control
% Reviewer sees the same info for all clips
clips = cell(numel(all_idx), 5);
for i = 1:numel(all_idx)
    d = all_devs(all_idx(i));
    clips{i, 1} = d.expdir;
    clips{i, 2} = d.fly;
    clips{i, 3} = [];
    clips{i, 4} = [];
    clips{i, 5} = sprintf('%s_oc%d_wk%d', d.sex, round(d.pct_oc), round(d.pct_walk));
end

% Shuffle order
shuffle_order = randperm(size(clips, 1));
clips = clips(shuffle_order, :);

fprintf('Review set: %d candidates (<10%% walk) + %d random normals = %d total\n', ...
    numel(candidate_idx), n_random, size(clips, 1));

%% GIF generation parameters
crop_size = 128;
frames_per_segment = 66;  % ~66 frames per window, 3 windows = ~200 frames total
frame_step = 3;           % sample every Nth frame
delay_time = 0.05;        % seconds between GIF frames
pause_frames = 5;         % black frames between segments

n_clips = size(clips, 1);
for ci = 1:n_clips
    expdir = clips{ci, 1};
    fly = clips{ci, 2};
    label = clips{ci, 5};

    [~, ename] = fileparts(expdir);
    fprintf('[%d/%d] %s fly %d...', ci, n_clips, ename, fly);

    % Check movie exists
    moviefile = fullfile(expdir, 'movie.ufmf');
    if ~exist(moviefile, 'file')
        fprintf(' SKIPPED (no movie)\n');
        continue;
    end

    % Load trx
    trxdata = load(fullfile(expdir, 'registered_trx.mat'));
    trx = trxdata.trx;
    if fly > numel(trx)
        fprintf(' SKIPPED (fly %d > %d)\n', fly, numel(trx));
        continue;
    end
    ff = trx(fly).firstframe;
    ef = trx(fly).endframe;
    nf_total = ef - ff + 1;

    % Load onceiling scores
    scorefile = fullfile(expdir, 'scores_onceiling_resnet_v2.mat');
    if exist(scorefile, 'file')
        soc = load(scorefile);
        has_oc = true;
    else
        has_oc = false;
    end

    % Open movie
    [readfcn, ~, fid] = get_readframe_fcn(moviefile);

    % Compute 3 window centers at 25%, 50%, 75% of trajectory
    win_length = frames_per_segment * frame_step;
    centers = round([0.25 0.50 0.75] * nf_total) + ff;
    windows = zeros(3, 2);
    for wi = 1:3
        s = max(ff, centers(wi) - win_length/2);
        e = min(ef, s + win_length);
        windows(wi, :) = [s e];
    end

    % Build filename with first and last frame range
    gif_file = fullfile(outdir, sprintf('%s_fly%02d_fr%dto%d_%dto%d_%dto%d_%s.gif', ...
        ename, fly, windows(1,1), windows(1,2), windows(2,1), windows(2,2), ...
        windows(3,1), windows(3,2), label));

    % Skip if already generated
    if exist(gif_file, 'file')
        fprintf(' exists, skipping\n');
        continue;
    end

    % Pre-render trajectory map: full arena with fly path
    all_x = trx(fly).x;
    all_y = trx(fly).y;
    map_size = crop_size * 2;  % same height as crop panel

    % Get arena bounds from first frame
    frame1 = readfcn(ff);
    [arena_h, arena_w] = size(frame1);
    scale = map_size / max(arena_h, arena_w);

    frame_count = 0;
    for wi = 1:3
        sample_frames = windows(wi,1):frame_step:windows(wi,2);
        sample_frames = sample_frames(1:min(frames_per_segment, numel(sample_frames)));

        for fi = 1:numel(sample_frames)
            fr = sample_frames(fi);
            frame = readfcn(fr);
            [h, w] = size(frame);

            tidx = fr - ff + 1;
            if tidx < 1 || tidx > numel(trx(fly).x), continue; end
            x = round(trx(fly).x(tidx));
            y = round(trx(fly).y(tidx));

            % --- Left panel: fly crop ---
            y0 = max(1, y - crop_size/2);
            x0 = max(1, x - crop_size/2);
            y1 = min(h, y0 + crop_size - 1);
            x1 = min(w, x0 + crop_size - 1);
            crop = frame(y0:y1, x0:x1);
            if size(crop,1) < crop_size || size(crop,2) < crop_size
                crop = imresize(crop, [crop_size crop_size]);
            end
            crop = imresize(crop, 2, 'nearest');

            % White border if onceiling
            if has_oc
                oc_tidx = fr - ff + 1;
                if oc_tidx >= 1 && oc_tidx <= numel(soc.allScores.scores{fly})
                    is_oc = soc.allScores.scores{fly}(oc_tidx) > 0;
                else
                    is_oc = false;
                end
            else
                is_oc = false;
            end
            if is_oc
                crop(1:3, :) = 255; crop(end-2:end, :) = 255;
                crop(:, 1:3) = 255; crop(:, end-2:end) = 255;
            end

            % Segment indicator bar
            bar_vals = [150 200 250];
            crop(4:6, :) = bar_vals(wi);

            % Add frame number text using insertText
            crop_rgb = repmat(crop, [1 1 3]);
            crop_rgb = insertText(crop_rgb, [5 map_size-18], sprintf('fr %d', fr), ...
                'FontSize', 12, 'TextColor', 'black', 'BoxColor', 'white', 'BoxOpacity', 0.8);
            crop_gray = rgb2gray(crop_rgb);

            % --- Right panel: trajectory map (black on white) ---
            tmap = ones(map_size, map_size, 'uint8') * 255;  % white background

            % Draw full trajectory in black
            px = round(all_x * scale);
            py = round(all_y * scale);
            px = max(1, min(map_size, px));
            py = max(1, min(map_size, py));
            for ti = 1:numel(px)
                tmap(py(ti), px(ti)) = 0;
            end

            % Current position as gray dot
            cx = px(tidx); cy = py(tidx);
            r = 3;
            for dy = -r:r
                for dx = -r:r
                    if dx^2 + dy^2 <= r^2
                        yy = max(1, min(map_size, cy+dy));
                        xx = max(1, min(map_size, cx+dx));
                        tmap(yy, xx) = 120;
                    end
                end
            end

            % Add total frame count label
            tmap_rgb = repmat(tmap, [1 1 3]);
            tmap_rgb = insertText(tmap_rgb, [5 map_size-18], sprintf('%d frames', nf_total), ...
                'FontSize', 12, 'TextColor', 'black', 'BoxColor', 'white', 'BoxOpacity', 0.8);
            tmap = rgb2gray(tmap_rgb);

            % Combine panels side by side
            combined = [crop_gray, tmap];

            % Write to GIF
            frame_count = frame_count + 1;
            if frame_count == 1
                imwrite(combined, gif_file, 'gif', 'Loopcount', inf, 'DelayTime', delay_time);
            else
                imwrite(combined, gif_file, 'gif', 'WriteMode', 'append', 'DelayTime', delay_time);
            end
        end

        % Pause frames between segments
        if wi < 3
            black = zeros(map_size, map_size * 2, 'uint8');
            for pi = 1:pause_frames
                imwrite(black, gif_file, 'gif', 'WriteMode', 'append', 'DelayTime', 0.15);
                frame_count = frame_count + 1;
            end
        end
    end

    if exist('fid', 'var') && fid > 0
        try fclose(fid); catch, end
    end

    fprintf(' %d frames\n', frame_count);
end

%% List GIF filenames for review.html
gif_listing = dir(fullfile(outdir, '*.gif'));
fprintf('\nGenerated %d GIFs in %s\n', numel(gif_listing), outdir);
fprintf('\nGIF filenames:\n');
for gi = 1:numel(gif_listing)
    fprintf('  "%s",\n', gif_listing(gi).name);
end

fprintf('Done.\n');

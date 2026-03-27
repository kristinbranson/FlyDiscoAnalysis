%% script_make_walk_gifs.m
% Generate GIFs showing individual walks with color-coded border indicating
% when the union of onceiling | nottracking classifiers is positive.
%
% Green border = on floor (good frames)
% Red border = not on floor (onceiling or nottracking positive)
%
% Selects ~20 walks spread across frac_notonfloor bins (0, ~0.1, ~0.5, ~1.0)
% from experiments with 25-90% median on-ceiling.

modpath;

%% Config
rootdir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20260225_VNC23controldata4onceilingwalkexplore';
outdir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Locomotion_analysis/claudeplots_20260324_walk_classifier_filtering/walk_gifs_allexp';
if ~isfolder(outdir), mkdir(outdir); end

crop_size = 128;       % half-width of crop around fly
frame_step = 2;        % sample every Nth frame within walk
delay_time = 0.05;     % seconds between GIF frames
border_width = 4;      % pixels for color border

% Target bins by frac_notonfloor — focus on the ambiguous middle range
frac_bin_edges = [0 0.001; 0.01 0.10; 0.10 0.25; 0.25 0.50; 0.50 0.75; 0.75 0.99; 0.999 1.0];
frac_bin_labels = {'onfloor100%', 'onfloor90-99%', 'onfloor75-90%', 'onfloor50-75%', 'onfloor25-50%', 'onfloor1-25%', 'onfloor0%'};
walks_per_bin = [3 8 8 8 8 8 3];

%% Load cached walk data
cache_file = fullfile(rootdir, 'walk_classifier_overlap_cache_allexp.mat');
fprintf('Loading walk cache...\n');
C = load(cache_file);
sel_expnames = C.sel_expnames;

%% Build sex lookup from perfly_classifier_overlap.mat
fprintf('Loading sex info from perfly_classifier_overlap.mat...\n');
F = load(fullfile(rootdir, 'perfly_classifier_overlap.mat'), 'fly_data');
% Build map: expname_flyN -> sex
sex_map = containers.Map();
for fi = 1:numel(F.fly_data)
    key = sprintf('%s_fly%d', F.fly_data(fi).expname, F.fly_data(fi).fly);
    sex_map(key) = F.fly_data(fi).sex;
end

% Look up sex for each walk in the cache
all_sex = cell(1, numel(C.all_fly));
for wi = 1:numel(C.all_fly)
    key = sprintf('%s_fly%d', sel_expnames{C.all_exp_idx(wi)}, C.all_fly(wi));
    if sex_map.isKey(key)
        all_sex{wi} = sex_map(key);
    else
        all_sex{wi} = '?';
    end
end

%% Select walks from each frac_notonfloor bin
% Enforce diversity: max 1 walk per (experiment, fly) pair, mix M and F
rng(42);
selected = struct('exp_idx', {}, 'walk_idx_in_exp', {}, 'frac', {}, 'bin_label', {}, ...
    'global_idx', {}, 'fly', {}, 'velmag', {}, 'sex', {});

for bi = 1:size(frac_bin_edges, 1)
    lo = frac_bin_edges(bi, 1);
    hi = frac_bin_edges(bi, 2);
    in_bin = find(C.all_frac_notonfloor >= lo & C.all_frac_notonfloor <= hi & ~isnan(C.all_frac_notonfloor));

    if isempty(in_bin)
        fprintf('Bin %s: no walks found\n', frac_bin_labels{bi});
        continue;
    end

    % Shuffle candidates
    in_bin = in_bin(randperm(numel(in_bin)));

    % Pick walks enforcing: max 1 per (exp, fly), aim for ~half M half F
    target_n = walks_per_bin(bi);
    used_exp_fly = containers.Map();  % track used (exp,fly) pairs
    picked_M = {};
    picked_F = {};

    for ci = 1:numel(in_bin)
        idx = in_bin(ci);
        exp_fly_key = sprintf('%d_%d', C.all_exp_idx(idx), C.all_fly(idx));
        if used_exp_fly.isKey(exp_fly_key), continue; end

        s = struct();
        s.exp_idx = C.all_exp_idx(idx);
        s.frac = C.all_frac_notonfloor(idx);
        s.bin_label = frac_bin_labels{bi};
        same_exp = find(C.all_exp_idx == s.exp_idx);
        s.walk_idx_in_exp = find(same_exp == idx);
        s.global_idx = idx;
        s.fly = C.all_fly(idx);
        s.velmag = C.all_velmag(idx);
        s.sex = all_sex{idx};

        used_exp_fly(exp_fly_key) = true;

        if strcmp(s.sex, 'M')
            picked_M{end+1} = s; %#ok<SAGROW>
        else
            picked_F{end+1} = s; %#ok<SAGROW>
        end

        if numel(picked_M) + numel(picked_F) >= target_n * 3
            break;  % have enough candidates to select from
        end
    end

    % Balance M and F: take up to half from each, fill remainder from other
    half = ceil(target_n / 2);
    n_M = min(half, numel(picked_M));
    n_F = min(target_n - n_M, numel(picked_F));
    if n_M + n_F < target_n
        n_M = min(target_n - n_F, numel(picked_M));
    end

    bin_selected = [picked_M(1:n_M), picked_F(1:n_F)];
    for pi = 1:numel(bin_selected)
        selected(end+1) = bin_selected{pi}; %#ok<SAGROW>
    end
    fprintf('Bin %s [%.3f-%.3f]: selected %d walks (%d M, %d F)\n', ...
        frac_bin_labels{bi}, lo, hi, numel(bin_selected), n_M, n_F);
end

fprintf('Total walks to render: %d\n', numel(selected));

%% Generate GIFs
for si = 1:numel(selected)
    sel = selected(si);
    expname = sel_expnames{sel.exp_idx};
    expdir = fullfile(rootdir, expname);

    frac_onfloor = 1 - sel.frac;
    fprintf('[%d/%d] %s fly %d %s onfloor=%.0f%% (%s)...', ...
        si, numel(selected), expname, sel.fly, sel.sex, frac_onfloor*100, sel.bin_label);

    % Load walk struct to get walk boundaries
    W = load(fullfile(expdir, 'locomotion_walkstruct.mat'), 'walk_struct_OFF');
    ws = W.walk_struct_OFF;

    % Find the walk for this fly matching global index
    % Reconstruct: find all walks from this experiment in the flat array
    same_exp = find(C.all_exp_idx == sel.exp_idx);
    local_idx = find(same_exp == sel.global_idx);
    if isempty(local_idx)
        fprintf(' SKIP (index mismatch)\n');
        continue;
    end

    walk_t0 = ws.walk_t0(local_idx);
    walk_t1 = ws.walk_t1(local_idx);
    walk_fly = ws.fly(local_idx);

    % Load movie
    moviefile = fullfile(expdir, 'movie.ufmf');
    if ~exist(moviefile, 'file')
        fprintf(' SKIP (no movie)\n');
        continue;
    end

    % Load trx
    trxdata = load(fullfile(expdir, 'registered_trx.mat'));
    trx = trxdata.trx;
    if walk_fly > numel(trx)
        fprintf(' SKIP (fly %d > %d)\n', walk_fly, numel(trx));
        continue;
    end
    ff = trx(walk_fly).firstframe;

    % Load classifier scores
    OC = load(fullfile(expdir, 'scores_onceiling_resnet_v2.mat'), 'allScores');
    NT = load(fullfile(expdir, 'scores_nottracking_apt.mat'), 'allScores');
    pp_oc = OC.allScores.postprocessed{walk_fly};
    pp_nt = NT.allScores.postprocessed{walk_fly};

    % Open movie
    [readfcn, ~, fid] = get_readframe_fcn(moviefile);

    % Convert walk frame indices (trajectory-relative) to movie frames
    movie_t0 = walk_t0 + ff - 1;
    movie_t1 = walk_t1 + ff - 1;
    sample_frames = movie_t0:frame_step:movie_t1;

    frac_onfloor = round((1 - sel.frac) * 100);
    gif_file = fullfile(outdir, sprintf('%s_fly%02d_%s_walk_fr%dto%d_onfloor%dpct.gif', ...
        expname, walk_fly, sel.sex, walk_t0, walk_t1, frac_onfloor));

    if exist(gif_file, 'file')
        fprintf(' exists, skipping\n');
        if exist('fid', 'var') && fid > 0, try fclose(fid); catch, end, end
        continue;
    end

    frame_count = 0;
    for fi = 1:numel(sample_frames)
        fr = sample_frames(fi);

        frame = readfcn(fr);
        [h, w] = size(frame);

        tidx = fr - ff + 1;
        if tidx < 1 || tidx > numel(trx(walk_fly).x), continue; end
        x = round(trx(walk_fly).x(tidx));
        y = round(trx(walk_fly).y(tidx));

        % Crop around fly
        y0 = max(1, y - crop_size/2);
        x0 = max(1, x - crop_size/2);
        y1 = min(h, y0 + crop_size - 1);
        x1 = min(w, x0 + crop_size - 1);
        crop = frame(y0:y1, x0:x1);
        if size(crop,1) < crop_size || size(crop,2) < crop_size
            crop = imresize(crop, [crop_size crop_size]);
        end
        crop = imresize(crop, 2, 'nearest');
        cs = size(crop, 1);  % crop size after resize

        % Check classifier at this frame (trajectory index)
        is_onfloor = true;
        if tidx >= 1 && tidx <= min(numel(pp_oc), numel(pp_nt))
            oc_val = pp_oc(tidx);
            nt_val = pp_nt(tidx);
            if isnan(oc_val), oc_val = 0; end
            if isnan(nt_val), nt_val = 0; end
            is_onfloor = ~((oc_val > 0) || (nt_val > 0));
        end

        % Create RGB crop with colored border based on per-frame classifier
        crop_rgb = repmat(crop, [1 1 3]);
        bw = border_width;
        if is_onfloor
            % Green border = on floor
            crop_rgb(1:bw, :, 1) = 0; crop_rgb(1:bw, :, 2) = 200; crop_rgb(1:bw, :, 3) = 0;
            crop_rgb(end-bw+1:end, :, 1) = 0; crop_rgb(end-bw+1:end, :, 2) = 200; crop_rgb(end-bw+1:end, :, 3) = 0;
            crop_rgb(:, 1:bw, 1) = 0; crop_rgb(:, 1:bw, 2) = 200; crop_rgb(:, 1:bw, 3) = 0;
            crop_rgb(:, end-bw+1:end, 1) = 0; crop_rgb(:, end-bw+1:end, 2) = 200; crop_rgb(:, end-bw+1:end, 3) = 0;
        else
            % Red border = not on floor (onceiling or nottracking)
            crop_rgb(1:bw, :, 1) = 255; crop_rgb(1:bw, :, 2) = 0; crop_rgb(1:bw, :, 3) = 0;
            crop_rgb(end-bw+1:end, :, 1) = 255; crop_rgb(end-bw+1:end, :, 2) = 0; crop_rgb(end-bw+1:end, :, 3) = 0;
            crop_rgb(:, 1:bw, 1) = 255; crop_rgb(:, 1:bw, 2) = 0; crop_rgb(:, 1:bw, 3) = 0;
            crop_rgb(:, end-bw+1:end, 1) = 255; crop_rgb(:, end-bw+1:end, 2) = 0; crop_rgb(:, end-bw+1:end, 3) = 0;
        end

        % Add text: frame number, sex, and frac_onfloor for this walk
        frac_onfloor_pct = round((1 - sel.frac) * 100);
        crop_rgb = insertText(crop_rgb, [5 cs-20], sprintf('fr %d  %s  onfloor: %d%%', fr, sel.sex, frac_onfloor_pct), ...
            'FontSize', 12, 'TextColor', 'white', 'BoxColor', 'black', 'BoxOpacity', 0.7);

        % Write to GIF
        frame_count = frame_count + 1;
        [ind, cmap] = rgb2ind(crop_rgb, 256);
        if frame_count == 1
            imwrite(ind, cmap, gif_file, 'gif', 'Loopcount', inf, 'DelayTime', delay_time);
        else
            imwrite(ind, cmap, gif_file, 'gif', 'WriteMode', 'append', 'DelayTime', delay_time);
        end
    end

    if exist('fid', 'var') && fid > 0
        try fclose(fid); catch, end
    end

    fprintf(' %d frames\n', frame_count);
end

fprintf('\nGenerated GIFs in %s\n', outdir);
fprintf('Done.\n');

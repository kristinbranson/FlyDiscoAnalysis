%% script_make_fly_gifs
% Generate short GIFs of cropped fly videos for reviewing classifier outliers.
% Reads movie.ufmf and registered_trx.mat, crops around the fly, saves as GIF.
%
% Edit the 'clips' cell array below to specify which fly/frame ranges to review.

modpath;

outdir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Locomotion_analysis/claudeplots_20260319_onceiling_allVNC/review_gifs_damaged2/';
if ~isfolder(outdir), mkdir(outdir); end

% Each row: {expdir, fly, start_frame, end_frame, label}
% start_frame/end_frame are movie frames (1-indexed)
% Set start_frame=[] to auto-pick a window from the middle of the trajectory
% Candidate damaged: <3% walk, >30pp deviant, not previously reviewed
clips = {
    '/groups/branson/bransonlab/flydisco_data/VNC2_EXT_VGLUT-GAL4_RigA_20230912T114817', 1, [], [], 'dmg2_m_oc53_wk0'
    '/groups/branson/bransonlab/flydisco_data/VNC2_EXT_VGLUT-GAL4_RigB_20230912T114850', 1, [], [], 'dmg2_m_oc63_wk0'
    '/groups/branson/bransonlab/flydisco_data/VNC2_EXT_VGLUT-GAL4_RigB_20231114T124058', 2, [], [], 'dmg2_m_oc42_wk1'
    '/groups/branson/bransonlab/flydisco_data/VNC2_EXT_VGLUT-GAL4_RigC_20220921T104115', 1, [], [], 'dmg2_m_oc54_wk0'
    '/groups/branson/bransonlab/flydisco_data/VNC2_EXT_VGLUT-GAL4_RigC_20230614T101434', 1, [], [], 'dmg2_m_oc38_wk1'
    '/groups/branson/bransonlab/flydisco_data/VNC2_JRC_SS86693_RigC_20231114T101936', 10, [], [], 'dmg2_f_oc42_wk3'
    '/groups/branson/bransonlab/flydisco_data/VNC2_JRC_SS89969_RigC_20220906T121622', 1, [], [], 'dmg2_m_oc58_wk0'
    '/groups/branson/bransonlab/flydisco_data/VNC2_JRC_SS92103_RigC_20220906T120105', 3, [], [], 'dmg2_m_oc65_wk0'
    '/groups/branson/bransonlab/flydisco_data/VNC2_JRC_SS94914_RigD_20220913T105042', 6, [], [], 'dmg2_f_oc44_wk0'
    '/groups/branson/bransonlab/flydisco_data/VNC2_JRC_SS94914_RigD_20220913T105042', 9, [], [], 'dmg2_f_oc79_wk0'
    '/groups/branson/bransonlab/flydisco_data/VNC3_EXT_VGLUT-GAL4_RigB_20240625T114612', 1, [], [], 'dmg2_m_oc79_wk0'
    '/groups/branson/bransonlab/flydisco_data/VNC3_EXT_VGLUT-GAL4_RigB_20240903T120900', 10, [], [], 'dmg2_f_oc38_wk0'
    '/groups/branson/bransonlab/flydisco_data/VNC3_EXT_VGLUT-GAL4_RigB_20240919T123524', 1, [], [], 'dmg2_m_oc71_wk0'
    '/groups/branson/bransonlab/flydisco_data/VNC3_JRC_SS104403_RigC_20240820T104412', 2, [], [], 'dmg2_m_oc69_wk0'
    '/groups/branson/bransonlab/flydisco_data/VNC3_JRC_SS41570_RigA_20240702T121826', 9, [], [], 'dmg2_f_oc58_wk0'
    '/groups/branson/bransonlab/flydisco_data/VNC3_YNA_K_162984_RigC_20240702T123229', 10, [], [], 'dmg2_f_oc44_wk0'
    '/groups/branson/bransonlab/flydisco_data/VNC_EXT_VGLUT-GAL4_RigB_20210405T145151', 1, [], [], 'dmg2_m_oc44_wk0'
    '/groups/branson/bransonlab/flydisco_data/VNC_EXT_VGLUT-GAL4_RigB_20211201T160558', 1, [], [], 'dmg2_m_oc79_wk0'
    '/groups/branson/bransonlab/flydisco_data/VNC_EXT_VGLUT-GAL4_RigC_20210405T143235', 1, [], [], 'dmg2_m_oc41_wk0'
    '/groups/branson/bransonlab/flydisco_data/VNC_JRC_SS32312_RigB_20210628T135738', 1, [], [], 'dmg2_m_oc56_wk0'
    '/groups/branson/bransonlab/flydisco_data/VNC_JRC_SS32371_RigC_20210601T140734', 7, [], [], 'dmg2_f_oc47_wk0'
    '/groups/branson/bransonlab/flydisco_data/VNC_JRC_SS34537_RigD_20210426T130239', 6, [], [], 'dmg2_f_oc64_wk0'
    '/groups/branson/bransonlab/flydisco_data/VNC_JRC_SS34687_RigC_20210929T130756', 3, [], [], 'dmg2_m_oc61_wk1'
    '/groups/branson/bransonlab/flydisco_data/VNC_JRC_SS36194_RigD_20210428T133032', 1, [], [], 'dmg2_m_oc46_wk3'
    '/groups/branson/bransonlab/flydisco_data/VNC_JRC_SS37633_RigA_20210517T140104', 1, [], [], 'dmg2_m_oc73_wk0'
    '/groups/branson/bransonlab/flydisco_data/VNC_JRC_SS45635_RigB_20210413T135253', 1, [], [], 'dmg2_m_oc48_wk0'
    '/groups/branson/bransonlab/flydisco_data/VNC_JRC_SS46269_RigC_20210413T132456', 4, [], [], 'dmg2_m_oc33_wk0'
    '/groups/branson/bransonlab/flydisco_data/VNC_JRC_SS46730_RigD_20210426T150635', 1, [], [], 'dmg2_m_oc69_wk0'
    '/groups/branson/bransonlab/flydisco_data/VNC_JRC_SS49220_RigB_20210419T144428', 1, [], [], 'dmg2_m_oc67_wk0'
    '/groups/branson/bransonlab/flydisco_data/VNC_JRC_SS53029_RigB_20210427T141025', 10, [], [], 'dmg2_f_oc70_wk0'
    '/groups/branson/bransonlab/flydisco_data/VNC_JRC_SS57969_RigC_20210510T142434', 6, [], [], 'dmg2_f_oc37_wk3'
    '/groups/branson/bransonlab/flydisco_data/VNC_JRC_SS57980_RigB_20211028T133745', 1, [], [], 'dmg2_m_oc58_wk0'
    '/groups/branson/bransonlab/flydisco_data/VNC_JRC_SS61045_RigC_20210928T151244', 1, [], [], 'dmg2_m_oc33_wk0'
    '/groups/branson/bransonlab/flydisco_data/VNC_JRC_SS65710_RigA_20210511T151754', 1, [], [], 'dmg2_m_oc80_wk0'
    '/groups/branson/bransonlab/flydisco_data/VNC_JRC_SS94679_RigD_20211202T154538', 8, [], [], 'dmg2_f_oc67_wk0'
};

crop_size = 128;
n_frames_gif = 200;   % number of frames in GIF
frame_step = 3;       % sample every Nth frame
delay_time = 0.05;    % seconds between GIF frames

for ci = 1:size(clips, 1)
    expdir = clips{ci, 1};
    fly = clips{ci, 2};
    start_fr = clips{ci, 3};
    end_fr = clips{ci, 4};
    label = clips{ci, 5};

    [~, ename] = fileparts(expdir);
    fprintf('Processing %s fly %d...\n', ename, fly);

    % Load trx
    trx = load(fullfile(expdir, 'registered_trx.mat'));
    trx = trx.trx;
    ff = trx(fly).firstframe;
    ef = trx(fly).endframe;

    % Auto-pick window from middle of trajectory
    if isempty(start_fr)
        win = n_frames_gif * frame_step;
        mid = round((ff + ef) / 2);
        start_fr = max(ff, mid - win/2);
        end_fr = min(ef, start_fr + win);
    end

    % Load onceiling scores for border indicator
    soc = load(fullfile(expdir, 'scores_onceiling_resnet_v2.mat'));

    % Open movie
    [readfcn, nframes, fid] = get_readframe_fcn(fullfile(expdir, 'movie.ufmf'));

    % Sample frames
    sample_frames = start_fr:frame_step:min(end_fr, ef);
    sample_frames = sample_frames(1:min(n_frames_gif, numel(sample_frames)));

    gif_file = fullfile(outdir, sprintf('%s_fly%02d_fr%dto%d_%s.gif', ename, fly, sample_frames(1), sample_frames(end), label));

    for fi = 1:numel(sample_frames)
        fr = sample_frames(fi);
        frame = readfcn(fr);
        [h, w] = size(frame);

        % Get fly position
        tidx = fr - ff + 1;
        if tidx < 1 || tidx > numel(trx(fly).x)
            continue;
        end
        x = round(trx(fly).x(tidx));
        y = round(trx(fly).y(tidx));

        % Simple centered crop (no rotation)
        y0 = max(1, y - crop_size/2);
        x0 = max(1, x - crop_size/2);
        y1 = min(h, y0 + crop_size - 1);
        x1 = min(w, x0 + crop_size - 1);
        crop = frame(y0:y1, x0:x1);

        % Pad if needed
        if size(crop,1) < crop_size || size(crop,2) < crop_size
            crop = imresize(crop, [crop_size crop_size]);
        end

        % Scale up 2x for easier viewing
        crop = imresize(crop, 2, 'nearest');

        % Check onceiling score for this frame
        oc_tidx = fr - ff + 1;
        if exist('soc', 'var') && oc_tidx >= 1 && oc_tidx <= numel(soc.allScores.scores{fly})
            is_oc = soc.allScores.scores{fly}(oc_tidx) > 0;
        else
            is_oc = false;
        end

        % Add red border if onceiling
        if is_oc
            crop(1:2, :) = 255;
            crop(end-1:end, :) = 255;
            crop(:, 1:2) = 255;
            crop(:, end-1:end) = 255;
        end

        % Write to GIF
        if fi == 1
            imwrite(crop, gif_file, 'gif', 'Loopcount', inf, 'DelayTime', delay_time);
        else
            imwrite(crop, gif_file, 'gif', 'WriteMode', 'append', 'DelayTime', delay_time);
        end
    end

    if exist('fid', 'var') && fid > 0
        try fclose(fid); catch, end
    end

    fprintf('  Saved %s (%d frames)\n', gif_file, numel(sample_frames));
    clear soc;
end
%% Print GIF filenames for updating review.html GIFS array
fprintf('\nGIF filenames for review.html:\n');
gif_listing = dir(fullfile(outdir, '*.gif'));
for gi = 1:numel(gif_listing)
    fprintf('  "%s",\n', gif_listing(gi).name);
end

fprintf('Done.\n');

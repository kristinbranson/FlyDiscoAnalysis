%% script_compile_allVNC_onceiling
% Compile per-experiment and per-fly onceiling % for all VNC/VNC2/VNC3 experiments.
% Reads scores_onceiling_resnet_v2.mat from each production experiment directory.
%
% Input:
%   expdirlist_allVNC_production_8663.txt — list of production experiment paths
%
% Output:
%   allVNC_onceiling_resnet_v2_summary.mat containing:
%     expnames      — cell array of experiment directory names
%     expdirs       — cell array of full paths
%     pct_oc        — per-experiment % onceiling (pooled across all fly-frames)
%     dates_dt      — datetime per experiment
%     screen_type   — 'VNC', 'VNC2', or 'VNC3'
%     is_control    — logical, true for YNA_K_162984
%     fly_oc_per_exp — cell array, each entry is a 1xN vector of per-fly % onceiling

listfile = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/expdirlist_allVNC_production_8663.txt';
outfile = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/allVNC_onceiling_resnet_v2_summary.mat';

fid = fopen(listfile, 'r');
expdirs = textscan(fid, '%s'); expdirs = expdirs{1};
fclose(fid);

n = numel(expdirs);
pct_oc = nan(n, 1);
dates_dt = NaT(n, 1);
screen_type = cell(n, 1);
expnames = cell(n, 1);
is_control = false(n, 1);
fly_oc_per_exp = cell(n, 1);

for i = 1:n
    [~, ename] = fileparts(expdirs{i});
    expnames{i} = ename;

    % Parse date
    tokens = regexp(ename, '_(\d{8})T\d{6}$', 'tokens');
    if ~isempty(tokens)
        dates_dt(i) = datetime(tokens{1}{1}, 'InputFormat', 'yyyyMMdd');
    end

    % Screen type
    if startsWith(ename, 'VNC3_'), screen_type{i} = 'VNC3';
    elseif startsWith(ename, 'VNC2_'), screen_type{i} = 'VNC2';
    elseif startsWith(ename, 'VNC_'), screen_type{i} = 'VNC';
    end

    % Control line
    is_control(i) = contains(ename, 'YNA_K_162984');

    % Load scores and compute per-fly and per-experiment onceiling %
    s = load(fullfile(expdirs{i}, 'scores_onceiling_resnet_v2.mat'));
    nflies = numel(s.allScores.scores);
    fly_pcts = nan(1, nflies);
    all_sc = [];
    for f = 1:nflies
        sc = s.allScores.scores{f};
        fly_pcts(f) = 100 * mean(sc > 0);
        all_sc = [all_sc, sc(:)']; %#ok<AGROW>
    end
    fly_oc_per_exp{i} = fly_pcts;
    pct_oc(i) = 100 * mean(all_sc > 0);

    if mod(i, 500) == 0, fprintf('  %d/%d\n', i, n); end
end
fprintf('Done. %d experiments.\n', n);

save(outfile, 'expnames', 'expdirs', 'pct_oc', 'dates_dt', ...
    'screen_type', 'is_control', 'fly_oc_per_exp', '-v7.3');
fprintf('Saved to %s\n', outfile);

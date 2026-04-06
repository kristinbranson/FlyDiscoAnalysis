%% script_export_linemeans_csv.m
% Export normalized line means and z-scores as CSV files for collaborators.
%
% Outputs:
%   linemeans_normalized_20260331.csv — normalized means per line per stat
%   linemeans_zscore_20260331.csv — z-scores per line per stat

modpath;

%% load data
S = load('/groups/branson/bransonlab/flydisco_linelevel_VNC/CollectedVNC23PerFrameStats20260331.mat', ...
    'linestats', 'statfns', 'controlmean', 'controlstd', 'nlines');

outdir = '/groups/branson/home/robiea/Code_versioned/locomotion_analysis/Alice';
nlines = numel(S.linestats.line_names);
nstats = numel(S.statfns);

%% compute z-scores
normmeans = nan(nlines, nstats);
zscores = nan(nlines, nstats);
for si = 1:nstats
    fn = S.statfns{si};
    normmeans(:, si) = S.linestats.normmeans.(fn)';
    zscores(:, si) = (S.linestats.normmeans.(fn)' - S.controlmean(si)) / S.controlstd(si);
end

%% write normalized means CSV
outfile = fullfile(outdir, 'linemeans_normalized_20260402.csv');
fid = fopen(outfile, 'w');

% header
fprintf(fid, 'line');
for si = 1:nstats
    fprintf(fid, ',%s', S.statfns{si});
end
fprintf(fid, '\n');

% data
for li = 1:nlines
    fprintf(fid, '%s', S.linestats.line_names{li});
    for si = 1:nstats
        fprintf(fid, ',%.6g', normmeans(li, si));
    end
    fprintf(fid, '\n');
end
fclose(fid);
fprintf('Saved %s (%d lines x %d stats)\n', outfile, nlines, nstats);

%% write z-score CSV
outfile = fullfile(outdir, 'linemeans_zscore_20260402.csv');
fid = fopen(outfile, 'w');

% header
fprintf(fid, 'line');
for si = 1:nstats
    fprintf(fid, ',%s', S.statfns{si});
end
fprintf(fid, '\n');

% data
for li = 1:nlines
    fprintf(fid, '%s', S.linestats.line_names{li});
    for si = 1:nstats
        fprintf(fid, ',%.6g', zscores(li, si));
    end
    fprintf(fid, '\n');
end
fclose(fid);
fprintf('Saved %s (%d lines x %d stats)\n', outfile, nlines, nstats);

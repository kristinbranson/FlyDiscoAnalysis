%% script_submit_pvalue_jobs.m
% Submit ComputePValueBySampling jobs to the cluster.
% Run this from FlyDiscoAnalysis on a cluster login node after running
% ScriptComputePValuesVNC_20260331.m through the save(outmatfile,...) cell.
%
% Usage:
%   ssh -X login2
%   bsub -XF -Is -n1 -W 48:00 /bin/bash
%   cd /groups/branson/home/robiea/Code_versioned/FlyDiscoAnalysis
%   matlab
%   >> script_submit_pvalue_jobs

modpath;

%% parameters — update these to match your ScriptComputePValuesVNC script
outmatfile = '/groups/branson/home/robiea/Code_versioned/locomotion_analysis/Alice/ComputePValueBySamplingData20260402.mat';
resultsdir = '/groups/branson/home/robiea/Code_versioned/locomotion_analysis/Alice/ComputePValueBySamplingResults20260402';
forcecompute = false; % set to true to recompute all stats, even if result files exist

%% load and verify outmatfile has required fields
fprintf('Loading %s...\n', outmatfile);
D = load(outmatfile);

required_fields = {'statfns','nstats','setstats','allstats','allnframestotal',...
    'minnframes','minnexps','setiscontrol','nsets','maxnexps','exp2setidx',...
    'nlines','linestats','set2lineidx','nsamples','controlmean','controlstd','iscircstat'};

missing = setdiff(required_fields, fieldnames(D));
if ~isempty(missing)
    error('outmatfile is missing required fields: %s', strjoin(missing, ', '));
end
fprintf('All %d required fields present. nstats=%d, nlines=%d\n', ...
    numel(required_fields), D.nstats, D.nlines);

nstats = D.nstats;

%% create results directory
if ~exist(resultsdir,'dir'), mkdir(resultsdir); end

%% submit jobs (skip stats with existing results unless forcecompute)
n_submitted = 0;
n_skipped = 0;
for stati = 1:nstats
    resultfile = fullfile(resultsdir, sprintf('PvaluesForStat%03d.mat', stati));
    if ~forcecompute && exist(resultfile, 'file')
        n_skipped = n_skipped + 1;
        continue;
    end

    jobname = sprintf('pval_%d', stati);
    logfile = fullfile(resultsdir, sprintf('ComputePValuesBySampling_P%d.log', stati));

    cmd = sprintf(['bsub -J %s -o %s -n 1 -W 4:00 ', ...
        'matlab -nodisplay -nosplash -singleCompThread -r ', ...
        '"cd(''/groups/branson/home/robiea/Code_versioned/FlyDiscoAnalysis''); modpath; ', ...
        'cd(''/groups/branson/home/robiea/Code_versioned/locomotion_analysis/Alice''); ', ...
        'ComputePValueBySampling(%d, ''%s'', ''%s''); exit"'], ...
        jobname, logfile, stati, outmatfile, resultsdir);

    system(cmd);
    n_submitted = n_submitted + 1;
    if mod(n_submitted, 100) == 0
        fprintf('Submitted %d jobs so far (%d skipped)\n', n_submitted, n_skipped);
    end
end

fprintf('Submitted %d jobs, skipped %d with existing results. Use bjobs to check status.\n', n_submitted, n_skipped);

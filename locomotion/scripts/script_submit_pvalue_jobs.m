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

%% parameters — must match ScriptComputePValuesVNC_20260331.m
outmatfile = '/groups/branson/home/robiea/Code_versioned/locomotion_analysis/Alice/ComputePValueBySamplingData20260331.mat';
resultsdir = '/groups/branson/home/robiea/Code_versioned/locomotion_analysis/Alice/ComputePValueBySamplingResults20260331';

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

%% submit jobs
for stati = 1:nstats
    jobname = sprintf('pval_%d', stati);
    logfile = fullfile(resultsdir, sprintf('ComputePValuesBySampling_P%d.log', stati));

    cmd = sprintf(['bsub -J %s -o %s -n 1 -W 4:00 ', ...
        'matlab -nodisplay -nosplash -singleCompThread -r ', ...
        '"cd(''/groups/branson/home/robiea/Code_versioned/FlyDiscoAnalysis''); modpath; ', ...
        'cd(''/groups/branson/home/robiea/Code_versioned/locomotion_analysis/Alice''); ', ...
        'ComputePValueBySampling(%d, ''%s'', ''%s''); exit"'], ...
        jobname, logfile, stati, outmatfile, resultsdir);

    system(cmd);
    if mod(stati, 100) == 0
        fprintf('Submitted %d/%d jobs\n', stati, nstats);
    end
end

fprintf('All %d jobs submitted. Use bjobs to check status.\n', nstats);

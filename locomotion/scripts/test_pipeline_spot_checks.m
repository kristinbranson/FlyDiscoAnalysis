%% test_pipeline_spot_checks.m
% Spot checks for onfloor filtering pipeline before launching on full dataset.
%
% Tests:
%   1. Low-ceiling VNC2 experiment (filtering should be near no-op)
%   2. VNC experiment (different screen type)
%   3. VNC3 experiment (different screen type)
%   4. Visual sanity: AEPy filtered vs unfiltered on high-ceiling experiment
%   5. Graceful failure: missing score files
%
% Prerequisites:
%   - 20260326 settings directories must have dolocomotionmetrics,force
%   - Production experiments must have scores_nottracking.mat and
%     scores_onceiling_resnet_v2.mat

modpath;

settings_folder_path = '/groups/branson/home/robiea/Code_versioned/BransonFlyDiscoSettings/settings';
test_root = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20260327_testingOnFloorPipeline';
if ~isfolder(test_root), mkdir(test_root); end

% Pipeline parameters for local run
cluster_billing_account_name = 'branson';
user_name_for_configuration_purposes = 'bransonlab';
do_use_bqueue = false;
do_actually_submit_jobs = false;
do_try = true;
ssh_host_name = 'login2.int.janelia.org';

results = struct('name', {}, 'status', {}, 'details', {});

%% ========================================================================
%% TEST 1: Low-ceiling VNC2 (0.2% ceiling — filtering near no-op)
%% ========================================================================
fprintf('\n============================================================\n');
fprintf('TEST 1: Low-ceiling VNC2 experiment\n');
fprintf('============================================================\n');

source_exp = '/groups/branson/bransonlab/flydisco_data/VNC2_YNA_K_162984_RigA_20220525T094753';
protocol = '20260326_flybubble_LED_VNC2';
[~, ename] = fileparts(source_exp);
test_expdir = fullfile(test_root, ename);

try
    % Create hard copy
    if ~isfolder(test_expdir)
        fprintf('Creating hard copy of %s...\n', ename);
        SymbolicCopyExperimentDirectory(source_exp, test_root, 'dosoftlink', false);
    end

    % Run full pipeline
    fprintf('Running pipeline with protocol %s...\n', protocol);
    run_transfero_FlyDiscoPipeline_wrapper_on_experiment_list({test_expdir}, ...
        cluster_billing_account_name, user_name_for_configuration_purposes, ...
        do_use_bqueue, do_actually_submit_jobs, do_try, ssh_host_name, ...
        'settingsdir', settings_folder_path, 'analysis_protocol', protocol);

    % Check output
    outfile = fullfile(test_expdir, 'locostatsperexp_onfloor.mat');
    if exist(outfile, 'file')
        S = load(outfile);
        stats = S.locostatsperexp;
        fns = fieldnames(stats);
        fprintf('PASS: locostatsperexp_onfloor.mat produced with %d fields\n', numel(fns));

        % Check walk Z count — should be close to total walks for low-ceiling
        walkfield = 'velmag_ctr__walk__LEDoff__all';
        if isfield(stats, walkfield)
            fprintf('  velmag_ctr walk LEDoff Z = %d\n', stats.(walkfield).Z);
        end
        results(end+1) = struct('name', 'Test1_LowCeiling_VNC2', 'status', 'PASS', 'details', sprintf('%d fields', numel(fns)));
    else
        fprintf('FAIL: locostatsperexp_onfloor.mat not found\n');
        results(end+1) = struct('name', 'Test1_LowCeiling_VNC2', 'status', 'FAIL', 'details', 'output file missing');
    end
catch ME
    fprintf('FAIL: %s\n', ME.message);
    results(end+1) = struct('name', 'Test1_LowCeiling_VNC2', 'status', 'FAIL', 'details', ME.message);
end

%% ========================================================================
%% TEST 2: VNC experiment (different screen type)
%% ========================================================================
fprintf('\n============================================================\n');
fprintf('TEST 2: VNC experiment\n');
fprintf('============================================================\n');

source_exp = '/groups/branson/bransonlab/flydisco_data/VNC_YNA_K_162984_RigA_20210406T150414';
protocol = '20260326_flybubble_LED_VNC';
[~, ename] = fileparts(source_exp);
test_expdir = fullfile(test_root, ename);

try
    if ~isfolder(test_expdir)
        fprintf('Creating hard copy of %s...\n', ename);
        SymbolicCopyExperimentDirectory(source_exp, test_root, 'dosoftlink', false);
    end

    fprintf('Running pipeline with protocol %s...\n', protocol);
    run_transfero_FlyDiscoPipeline_wrapper_on_experiment_list({test_expdir}, ...
        cluster_billing_account_name, user_name_for_configuration_purposes, ...
        do_use_bqueue, do_actually_submit_jobs, do_try, ssh_host_name, ...
        'settingsdir', settings_folder_path, 'analysis_protocol', protocol);

    outfile = fullfile(test_expdir, 'locostatsperexp_onfloor.mat');
    if exist(outfile, 'file')
        S = load(outfile);
        stats = S.locostatsperexp;
        fns = fieldnames(stats);
        fprintf('PASS: locostatsperexp_onfloor.mat produced with %d fields\n', numel(fns));
        results(end+1) = struct('name', 'Test2_VNC', 'status', 'PASS', 'details', sprintf('%d fields', numel(fns)));
    else
        fprintf('FAIL: locostatsperexp_onfloor.mat not found\n');
        results(end+1) = struct('name', 'Test2_VNC', 'status', 'FAIL', 'details', 'output file missing');
    end
catch ME
    fprintf('FAIL: %s\n', ME.message);
    results(end+1) = struct('name', 'Test2_VNC', 'status', 'FAIL', 'details', ME.message);
end

%% ========================================================================
%% TEST 3: VNC3 experiment (different screen type)
%% ========================================================================
fprintf('\n============================================================\n');
fprintf('TEST 3: VNC3 experiment\n');
fprintf('============================================================\n');

source_exp = '/groups/branson/bransonlab/flydisco_data/VNC3_YNA_K_162984_RigA_20240806T115247';
protocol = '20260326_flybubble_LED_VNC3';
[~, ename] = fileparts(source_exp);
test_expdir = fullfile(test_root, ename);

try
    if ~isfolder(test_expdir)
        fprintf('Creating hard copy of %s...\n', ename);
        SymbolicCopyExperimentDirectory(source_exp, test_root, 'dosoftlink', false);
    end

    fprintf('Running pipeline with protocol %s...\n', protocol);
    run_transfero_FlyDiscoPipeline_wrapper_on_experiment_list({test_expdir}, ...
        cluster_billing_account_name, user_name_for_configuration_purposes, ...
        do_use_bqueue, do_actually_submit_jobs, do_try, ssh_host_name, ...
        'settingsdir', settings_folder_path, 'analysis_protocol', protocol);

    outfile = fullfile(test_expdir, 'locostatsperexp_onfloor.mat');
    if exist(outfile, 'file')
        S = load(outfile);
        stats = S.locostatsperexp;
        fns = fieldnames(stats);
        fprintf('PASS: locostatsperexp_onfloor.mat produced with %d fields\n', numel(fns));
        results(end+1) = struct('name', 'Test3_VNC3', 'status', 'PASS', 'details', sprintf('%d fields', numel(fns)));
    else
        fprintf('FAIL: locostatsperexp_onfloor.mat not found\n');
        results(end+1) = struct('name', 'Test3_VNC3', 'status', 'FAIL', 'details', 'output file missing');
    end
catch ME
    fprintf('FAIL: %s\n', ME.message);
    results(end+1) = struct('name', 'Test3_VNC3', 'status', 'FAIL', 'details', ME.message);
end

%% ========================================================================
%% TEST 4: Visual sanity — AEPy filtered vs unfiltered
%% ========================================================================
fprintf('\n============================================================\n');
fprintf('TEST 4: AEPy visual sanity check (high-ceiling experiment)\n');
fprintf('============================================================\n');

highceiling_dir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20260325_testingOnFloorfiltering/VNC2_YNA_K_162984_RigB_20230913T100609';

try
    orig = load(fullfile(highceiling_dir, 'locostatsperexp_original.mat'));
    filt = load(fullfile(highceiling_dir, 'locostatsperexp_onfloor_reference.mat'));
    if isfield(orig, 'locostatsperexp'), orig_stats = orig.locostatsperexp; else, orig_stats = orig; end
    if isfield(filt, 'locostatsperexp'), filt_stats = filt.locostatsperexp; else, filt_stats = filt; end

    % Compare AEPy stats
    aepy_field = 'AEPy__step__LEDoff__all';
    fprintf('  AEPy (LEDoff, all limbs):\n');
    fprintf('    Unfiltered: mean=%.4f, std=%.4f, Z=%d\n', ...
        orig_stats.(aepy_field).mean, orig_stats.(aepy_field).std, orig_stats.(aepy_field).Z);
    fprintf('    Filtered:   mean=%.4f, std=%.4f, Z=%d\n', ...
        filt_stats.(aepy_field).mean, filt_stats.(aepy_field).std, filt_stats.(aepy_field).Z);
    fprintf('    Z reduction: %d steps removed (%.1f%%)\n', ...
        orig_stats.(aepy_field).Z - filt_stats.(aepy_field).Z, ...
        100*(1 - filt_stats.(aepy_field).Z / orig_stats.(aepy_field).Z));

    % Per-limb comparison
    fprintf('\n  Per-limb AEPy means (LEDoff):\n');
    fprintf('  %-8s %12s %12s %8s\n', 'Limb', 'Unfiltered', 'Filtered', 'Delta');
    for limb = 1:6
        fn = sprintf('AEPy__step__LEDoff__limb%d', limb);
        om = orig_stats.(fn).mean;
        fm = filt_stats.(fn).mean;
        fprintf('  limb%-3d %12.4f %12.4f %8.4f\n', limb, om, fm, fm - om);
    end

    % Plot: per-limb AEPy comparison
    figure('Position', [100 100 800 400]);
    limb_labels = {'RF','RM','RH','LH','LM','LF'};
    orig_means = zeros(1,6);
    filt_means = zeros(1,6);
    for limb = 1:6
        fn = sprintf('AEPy__step__LEDoff__limb%d', limb);
        orig_means(limb) = orig_stats.(fn).mean;
        filt_means(limb) = filt_stats.(fn).mean;
    end
    bar_data = [orig_means; filt_means]';
    bar(bar_data);
    set(gca, 'XTickLabel', limb_labels, 'TickLabelInterpreter', 'none');
    legend({'Unfiltered', 'Onfloor filtered'}, 'Location', 'best');
    ylabel('AEPy (mm)', 'Interpreter', 'none');
    title('AEPy per limb: unfiltered vs onfloor filtered (high-ceiling exp)', 'Interpreter', 'none');

    plotdir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Locomotion_analysis/claudeplots_20260327_onfloor_spotchecks';
    if ~isfolder(plotdir), mkdir(plotdir); end
    saveas(gcf, fullfile(plotdir, 'AEPy_filtered_vs_unfiltered.png'));
    fprintf('\n  Plot saved to %s\n', plotdir);

    results(end+1) = struct('name', 'Test4_AEPy_Visual', 'status', 'PASS', 'details', 'plot saved');
catch ME
    fprintf('FAIL: %s\n', ME.message);
    results(end+1) = struct('name', 'Test4_AEPy_Visual', 'status', 'FAIL', 'details', ME.message);
end

%% ========================================================================
%% TEST 5: Graceful failure — missing score files
%% ========================================================================
fprintf('\n============================================================\n');
fprintf('TEST 5: Graceful failure with missing score files\n');
fprintf('============================================================\n');

missing_root = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20260327_testingOnFloorPipeline_missingscores';
if ~isfolder(missing_root), mkdir(missing_root); end

source_exp = '/groups/branson/bransonlab/flydisco_data/VNC2_YNA_K_162984_RigA_20220525T094753';
protocol = '20260326_flybubble_LED_VNC2';
[~, ename] = fileparts(source_exp);
missing_expdir = fullfile(missing_root, ename);

try
    % Create hard copy
    if ~isfolder(missing_expdir)
        fprintf('Creating hard copy for missing-scores test...\n');
        SymbolicCopyExperimentDirectory(source_exp, missing_root, 'dosoftlink', false);
    end

    % Rename score file to simulate missing
    scorefile = fullfile(missing_expdir, 'scores_nottracking.mat');
    bakfile = fullfile(missing_expdir, 'scores_nottracking.mat.bak');
    if exist(scorefile, 'file') && ~exist(bakfile, 'file')
        movefile(scorefile, bakfile);
        fprintf('Renamed scores_nottracking.mat -> .bak\n');
    end

    % Run pipeline with do_try=true — should fail gracefully
    fprintf('Running pipeline (expecting loco stage to fail gracefully)...\n');
    run_transfero_FlyDiscoPipeline_wrapper_on_experiment_list({missing_expdir}, ...
        cluster_billing_account_name, user_name_for_configuration_purposes, ...
        do_use_bqueue, do_actually_submit_jobs, do_try, ssh_host_name, ...
        'settingsdir', settings_folder_path, 'analysis_protocol', protocol);

    % Check: locostatsperexp_onfloor.mat should NOT exist
    outfile = fullfile(missing_expdir, 'locostatsperexp_onfloor.mat');
    if ~exist(outfile, 'file')
        fprintf('PASS: locostatsperexp_onfloor.mat correctly not produced\n');

        % Check ACC ran and logged the failure
        accfile = fullfile(missing_expdir, 'automatic_checks_complete_info.mat');
        if exist(accfile, 'file')
            acc = load(accfile);
            fprintf('  ACC file exists — pipeline continued past failed loco stage\n');
        end
        results(end+1) = struct('name', 'Test5_MissingScores', 'status', 'PASS', 'details', 'loco stage failed gracefully');
    else
        fprintf('FAIL: locostatsperexp_onfloor.mat was produced despite missing scores\n');
        results(end+1) = struct('name', 'Test5_MissingScores', 'status', 'FAIL', 'details', 'output produced unexpectedly');
    end

    % Restore renamed file
    if exist(bakfile, 'file')
        movefile(bakfile, scorefile);
        fprintf('Restored scores_nottracking.mat\n');
    end
catch ME
    fprintf('FAIL: %s\n', ME.message);
    % Restore renamed file on error
    if exist(bakfile, 'file')
        movefile(bakfile, scorefile);
    end
    results(end+1) = struct('name', 'Test5_MissingScores', 'status', 'FAIL', 'details', ME.message);
end

%% ========================================================================
%% SUMMARY
%% ========================================================================
fprintf('\n============================================================\n');
fprintf('SUMMARY\n');
fprintf('============================================================\n');

n_pass = 0;
n_fail = 0;
for i = 1:numel(results)
    fprintf('  %-30s %s  %s\n', results(i).name, results(i).status, results(i).details);
    if strcmp(results(i).status, 'PASS')
        n_pass = n_pass + 1;
    else
        n_fail = n_fail + 1;
    end
end
fprintf('\n%d PASS, %d FAIL out of %d tests\n', n_pass, n_fail, numel(results));

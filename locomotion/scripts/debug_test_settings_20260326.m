%% Test 20260326 settings update locally
% Tests that the new jaaba_detect_params.txt (with fastcomputepffs,1)
% and added nottracking_aptv7_nowindowdata.jab classifier work correctly.
%
% Runs pipeline locally (no cluster) on symlinked copies of one VNC, VNC2,
% and VNC3 production experiment.

modpath;

%% Create writable test copies
source_expdirs = {
    '/groups/branson/bransonlab/flydisco_data/VNC_YNA_K_162984_RigA_20210406T150414'
    '/groups/branson/bransonlab/flydisco_data/VNC2_YNA_K_162984_RigA_20220419T110856'
    '/groups/branson/bransonlab/flydisco_data/VNC3_YNA_K_162984_RigA_20240514T113219'
};
analysis_protocols = {
    '20260326_flybubble_LED_VNC'
    '20260326_flybubble_LED_VNC2'
    '20260326_flybubble_LED_VNC3'
};

test_root = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20260326_test_settings_update';
if ~isfolder(test_root), mkdir(test_root); end

test_expdirs = cell(size(source_expdirs));
for i = 1:numel(source_expdirs)
    SymbolicCopyExperimentDirectory(source_expdirs{i}, test_root);
    [~, ename] = fileparts(source_expdirs{i});
    test_expdirs{i} = fullfile(test_root, ename);
    fprintf('Created test copy: %s\n', test_expdirs{i});
end

%% Run pipeline locally on each experiment with its matching settings
cluster_billing_account_name = 'branson';
user_name_for_configuration_purposes = 'bransonlab';
settings_folder_path = '/groups/branson/home/robiea/Code_versioned/BransonFlyDiscoSettings/settings';

do_use_bqueue = false;
do_actually_submit_jobs = false;
do_try = true;
ssh_host_name = 'login2.int.janelia.org';

for i = 1:numel(test_expdirs)
    fprintf('\n========================================\n');
    fprintf('Running pipeline on %s\n', test_expdirs{i});
    fprintf('Analysis protocol: %s\n', analysis_protocols{i});
    fprintf('========================================\n');

    optional_argument_list = {'settingsdir', settings_folder_path, ...
                              'analysis_protocol', analysis_protocols{i}};

    folder_path_from_experiment_index = test_expdirs(i);

    run_transfero_FlyDiscoPipeline_wrapper_on_experiment_list(folder_path_from_experiment_index, ...
        cluster_billing_account_name, ...
        user_name_for_configuration_purposes, ...
        do_use_bqueue, ...
        do_actually_submit_jobs, ...
        do_try, ...
        ssh_host_name, ...
        optional_argument_list{:});
end

%% Script to explore correlations between phase difference and other metrics
% Creates symbolic copies of Dec 5 (low phase) and Dec 6 (high phase) experiments,
% runs LimbBoutAnalyzer pipeline, and saves intermediate per-fly data for correlation analysis.

modpath

%% Configuration
rootoutputdir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20260211_exploringphaseissue';
settingsdir = '/groups/branson/home/robiea/Code_versioned/BransonFlyDiscoSettings/settings';
analysis_protocol = '20251009_flybubble_LED_VNC2';

% Source experiment directories (Dec 5 = low phase, Dec 6 = high phase)
% These are the control line experiments from the phase investigation
% Paths from script_flywise_phase_comparison.m
source_expdirs = {
    % Dec 5 (low phase day)
    '/groups/branson/bransonlab/flydisco_data/VNC2_YNA_K_162984_RigA_20231205T114519'
    '/groups/branson/bransonlab/flydisco_data/VNC2_YNA_K_162984_RigB_20231205T114631'
    '/groups/branson/bransonlab/flydisco_data/VNC2_YNA_K_162984_RigC_20231205T114709'
    '/groups/branson/bransonlab/flydisco_data/VNC2_YNA_K_162984_RigD_20231205T114743'
    % Dec 6 (high phase day)
    '/groups/branson/bransonlab/flydisco_data/VNC2_YNA_K_162984_RigA_20231206T125420'
    '/groups/branson/bransonlab/flydisco_data/VNC2_YNA_K_162984_RigB_20231206T125502'
    '/groups/branson/bransonlab/flydisco_data/VNC2_YNA_K_162984_RigC_20231206T125605'
    '/groups/branson/bransonlab/flydisco_data/VNC2_YNA_K_162984_RigD_20231206T125700'
};

%% Create output directory if needed
if ~exist(rootoutputdir, 'dir')
    mkdir(rootoutputdir);
    fprintf('Created output directory: %s\n', rootoutputdir);
end

%% Create symbolic copies
fprintf('\n=== Creating symbolic copies ===\n');
out_expdirs = cell(size(source_expdirs));
for i = 1:numel(source_expdirs)
    expdir = source_expdirs{i};
    [~, expname] = fileparts(expdir);
    fprintf('Processing %d/%d: %s\n', i, numel(source_expdirs), expname);

    out_expdirs{i} = SymbolicCopyExperimentDirectory(expdir, rootoutputdir);
    fprintf('  -> %s\n', out_expdirs{i});
end

%% Run LimbBoutAnalyzer pipeline on each experiment
fprintf('\n=== Running LimbBoutAnalyzer pipeline ===\n');
for i = 2:numel(out_expdirs)
    expdir = out_expdirs{i};
    [~, expname] = fileparts(expdir);
    fprintf('\n--- Processing %d/%d: %s ---\n', i, numel(out_expdirs), expname);

    % Check if already processed
    resultsfile = fullfile(expdir, 'locomotionmetricsswingstanceboutstats.mat');
    if exist(resultsfile, 'file')
        fprintf('  Results already exist, skipping.\n');
        continue;
    end

    try
        % Initialize trx
        fprintf('  Initializing trx...\n');
        trx = FBATrx('analysis_protocol', analysis_protocol, 'settingsdir', settingsdir, ...
            'datalocparamsfilestr', 'dataloc_params.txt');
        trx.AddExpDir(expdir, 'dooverwrite', false, 'openmovie', false);

        % Load APT data
        aptfile = trx.dataloc_params.apttrkfilestr;
        aptdata = TrkFile.load(fullfile(expdir, aptfile));

        % Read stage params
        stageparamsfile = fullfile(trx.settingsdir, trx.analysis_protocol, trx.dataloc_params.locomotionmetricsparamsfilestr);
        stage_params = ReadParams(stageparamsfile);
        legtip_landmarknums = stage_params.legtip_landmarknums;

        % Load or compute tips_velmag
        tips_velmag_file = fullfile(expdir, 'tips_velmag.mat');
        if exist(tips_velmag_file, 'file')
            load(tips_velmag_file, 'tips_velmag');
        else
            error('tips_velmag.mat not found - needs to be computed first');
        end

        % Load or compute tips_pos_body
        tips_pos_body_file = fullfile(expdir, 'tips_pos_body.mat');
        if exist(tips_pos_body_file, 'file')
            load(tips_pos_body_file, 'tips_pos_body');
        else
            error('tips_pos_body.mat not found - needs to be computed first');
        end

        % Load or compute groundcontact
        groundcontact_file = fullfile(expdir, 'groundcontact.mat');
        if exist(groundcontact_file, 'file')
            load(groundcontact_file, 'groundcontact');
        else
            % Compute groundcontact
            gc_threshold_low = stage_params.gc_threshold_low;
            gc_threshold_high = stage_params.gc_threshold_high;
            pairs = stage_params.pairs;
            minimum_bout = stage_params.minimum_bout_groundcontact;
            [groundcontact] = compute_groundcontact(tips_velmag, 'pairs', pairs, ...
                'gc_threshold_low', gc_threshold_low, 'gc_threshold_high', gc_threshold_high, ...
                'minimum_bout', minimum_bout);
            save(groundcontact_file, 'groundcontact');
        end

        % Load walking scores and LED indicator
        [~, walking_scores] = LoadScoresFromFile(trx, 'scores_Walk2', 1);
        indicatordata = trx.getIndicatorLED(1);
        digitalindicator = indicatordata.indicatordigital;

        % Initialize LimbBoutAnalyzer
        fprintf('  Initializing LimbBoutAnalyzer...\n');
        loco_analyzer = LimbBoutAnalyzer(trx, aptdata, tips_pos_body, legtip_landmarknums, ...
            groundcontact, digitalindicator, walking_scores, ...
            'phase_methods', {'phasediff_hilbert'}, ...
            'expdir', expdir);

        % Run analysis
        fprintf('  Running analyzeBoutAndStimConditions...\n');
        tic;
        loco_analyzer.analyzeBoutAndStimConditions();
        t1 = toc;
        fprintf('    Done in %.2f s\n', t1);

        fprintf('  Running analyzeWalkAndStimConditions...\n');
        tic;
        loco_analyzer.analyzeWalkAndStimConditions();
        t2 = toc;
        fprintf('    Done in %.2f s\n', t2);

        % Save results
        fprintf('  Saving results...\n');
        loco_analyzer.saveResults();

        fprintf('  Total time: %.2f s\n', t1 + t2);

    catch ME
        fprintf('  ERROR: %s\n', ME.message);
        fprintf('  Stack trace:\n');
        for k = 1:numel(ME.stack)
            fprintf('    %s:%d\n', ME.stack(k).name, ME.stack(k).line);
        end
    end
end

fprintf('\n=== Done ===\n');
fprintf('Results saved to: %s\n', rootoutputdir);

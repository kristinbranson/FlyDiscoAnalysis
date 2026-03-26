% AR 20250812 with Claude.ai query using my code base to create class

classdef LimbBoutAnalyzer < handle

    properties
        % input data
        trx                    % FBATrx object with per-fly trajectories: trx(fly)
        aptdata               % aptdata from TrkFile.load - leg tip tracking per fly
        tips_pos_body         % tip positions  in body reference frame (cell array per fly)
        legtip_landmarknums   % keypoint numbers for leg tips
        groundcontact         % ground contact data from compute_groundcontact (trajectory format)
        limbBoutData          % bout data - computed from groundcontact in constructor
        digitalindicator      % LED stimulus indicator modified to be trajectory format (per-fly)
        walking_scores        % walking scores - trajectory format (per-fly)

        % parameters
        states = {'swing', 'stance', 'step'} %types of bouts
        pairs = [1,6; 2,5; 3,4] % limb pair idx into legtip_landmarks
        pairNames = {'Front', 'Mid', 'Rear'}
        binedges = 5:2:34 % histogram bin edges

        % Results storage
        restrictedBoutData = containers.Map()  % Store multiple restricted datasets
        boutMetrics = containers.Map()         % Store multiple bout metrics
        validFrames = containers.Map()   % store computed digtal signals for restricting analysis frames in walk features
        walkMetrics = containers.Map()         % Store muttiple walk metric datasets
        locoStatsPerExp = struct()              % Flat struct of per-experiment stats; LED condition encoded in field names (LEDon/LEDoff)
        walkStruct = struct()                  % Per-condition walk/step struct-of-arrays; fields: walkStruct.OFF.walk_struct, walkStruct.OFF.step_struct

        % Settings
        phase_methods = {'phaselag', 'phasediff_interp', 'phasediff_hilbert', 'phasediff_hilbert_global'}
        debug = false
        expdir = ''
        nflies = 0

        % Onfloor filtering
        do_onfloor_filtering = false        % If true, run onfloor filtering path
        onceiling_scores = {}               % Trajectory-format onceiling labels {fly}(1,nframes)
        nottracking_scores = {}             % Trajectory-format nottracking labels {fly}(1,nframes)
        frac_onfloor_threshold = 0.70       % Keep walks with frac_onfloor >= this threshold
        is_onfloor_perframe = {}             % Per-frame onfloor label {fly}(1,nframes): ~(onceiling | nottracking)
        walkonfloor_digital = {}            % Per-fly digital signal {fly}(1,nframes): true for frames in walks with frac_onfloor >= frac_onfloor_threshold
        walkFracOnfloor = containers.Map()  % Keyed by condition, double array (1 x nwalks)
        walkValidity = containers.Map()     % Keyed by condition, logical array (1 x nwalks)

    end

    methods
        function obj = LimbBoutAnalyzer(trx, aptdata, tips_pos_body, legtip_landmarknums, groundcontact, digitalindicator, walking_scores, varargin)
            % Constructor - automatically computes limb bout data and converts digitalindicator
            % Usage: analyzer = LimbBoutAnalyzer(trx, aptdata, tips_pos_body, legtip_landmarknums, groundcontact, digitalindicator, walking_scores, 'debug', true)
            %
            % Inputs:
            %   trx - FBATrx object with trx(fly) for each fly
            %   aptdata - from TrkFile.load - leg tip data per fly
            %   tips_pos_body - cell array {fly} with tip positions
            %   legtip_landmarknums - landmark numbers for leg tips
            %   limbBoutData - bout data (trajectory or movie format)
            %   digitalindicator - stimulus indicator (movie timeline, single vector)

            % Set required properties
            obj.trx = trx;
            obj.aptdata = aptdata;
            obj.tips_pos_body = tips_pos_body;
            obj.legtip_landmarknums = legtip_landmarknums;
            obj.walking_scores = walking_scores;
            obj.nflies = obj.trx.nflies;

            % Compute limb bout data from ground contact
            % fprintf('Computing limb bout data from ground contact...\n');
            obj.limbBoutData = limbSwingStanceStep(groundcontact);
            % fprintf('Limb bout data computation complete\n');

            % Convert digitalindicator to trajectory format if needed
            obj.digitalindicator = obj.convertDigitalIndicatorToTrajectory(digitalindicator);

            % Validate that all trajectory data has consistent number of flies
            obj.validateTrajectoryData();


            % Parse optional arguments
            p = inputParser;
            addParameter(p, 'debug', false, @islogical);
            addParameter(p, 'expdir', '', @ischar);
            addParameter(p, 'binedges', 5:2:34, @isnumeric);
            addParameter(p, 'phase_methods', {'phaselag', 'phasediff_interp', 'phasediff_hilbert', 'phasediff_hilbert_global'}, @iscell);
            addParameter(p, 'do_onfloor_filtering', false, @islogical);
            addParameter(p, 'onceiling_scores', {}, @iscell);
            addParameter(p, 'nottracking_scores', {}, @iscell);
            addParameter(p, 'frac_onfloor_threshold', 0.70, @isnumeric);
            parse(p, varargin{:});

            obj.debug = p.Results.debug;
            obj.expdir = p.Results.expdir;
            obj.binedges = p.Results.binedges;
            obj.phase_methods = p.Results.phase_methods;
            obj.do_onfloor_filtering = p.Results.do_onfloor_filtering;
            obj.onceiling_scores = p.Results.onceiling_scores;
            obj.nottracking_scores = p.Results.nottracking_scores;
            obj.frac_onfloor_threshold = p.Results.frac_onfloor_threshold;           

            % Initialize containers
            obj.restrictedBoutData = containers.Map();
            obj.boutMetrics = containers.Map();
            obj.validFrames = containers.Map();
            obj.walkMetrics = containers.Map();
            obj.locoStatsPerExp = struct();
            obj.walkFracOnfloor = containers.Map();
            obj.walkValidity = containers.Map();

            
            if obj.debug
            obj.displayDataStructure();
            end
        end

        function validateTrajectoryData(obj)
            % Validate that all trajectory format data is consistent

            % Check digitalindicator (should now be trajectory format)
            if iscell(obj.digitalindicator)
                assert(length(obj.digitalindicator) == obj.nflies, ...
                    'digitalindicator must have %d flies', obj.nflies);
            else
                error('digitalindicator conversion failed - should be cell array in trajectory format');
            end

            % Check walking_scores
            if iscell(obj.walking_scores)
                assert(length(obj.walking_scores) == obj.nflies, ...
                    'walking_scores must have %d flies', obj.nflies);
            else
                error('walking_scores must be cell array in trajectory format: walking_scores{fly}');
            end

            % Check limbBoutData
            assert(isstruct(obj.limbBoutData) && length(obj.limbBoutData) == obj.nflies, ...
                'limbBoutData must be struct array with %d flies', obj.nflies);

            % Check onfloor filtering scores (only if filtering is enabled)
            if obj.do_onfloor_filtering
                if isempty(obj.onceiling_scores) || isempty(obj.nottracking_scores)
                    error('do_onfloor_filtering is true but onceiling_scores or nottracking_scores are empty');
                end
                assert(iscell(obj.onceiling_scores) && length(obj.onceiling_scores) == obj.nflies, ...
                    'onceiling_scores must be cell array with %d flies', obj.nflies);
                assert(iscell(obj.nottracking_scores) && length(obj.nottracking_scores) == obj.nflies, ...
                    'nottracking_scores must be cell array with %d flies', obj.nflies);
            end
        end

        function createValidFrames_LED(obj)
            obj.validFrames('led_on_traj') = obj.digitalindicator;

            obj.validFrames('led_off_traj') = obj.invertTrajectoryIndicator(obj.digitalindicator);
        end


        function computeOnfloorSignal(obj)
            % Compute walkonfloor_digital: trajectory-format digital signal
            % where true = frame is in a walk with frac_onfloor >= threshold.
            %
            % For each fly:
            %   1. Compute is_onfloor_perframe = ~(onceiling | nottracking)
            %   2. Detect walk bouts from walking_scores
            %   3. For each walk, compute frac_onfloor = mean(is_onfloor_perframe(t0:t1))
            %   4. Set walkonfloor_digital(t0:t1) = true only for walks with frac_onfloor >= threshold
            %
            % Also stores is_onfloor_perframe for use by buildWalkStructForCondition
            % to compute frac_onfloor for condition-specific walks.

            if isempty(obj.onceiling_scores) || isempty(obj.nottracking_scores)
                error('LimbBoutAnalyzer:computeOnfloorSignal', ...
                    'onceiling_scores and nottracking_scores must be set before calling computeOnfloorSignal.');
            end

            fprintf('Computing walkonfloor_digital (threshold = %.0f%%)...\n', obj.frac_onfloor_threshold * 100);

            obj.is_onfloor_perframe = cell(1, obj.nflies);
            obj.walkonfloor_digital = cell(1, obj.nflies);
            total_walks = 0;
            valid_walks = 0;

            for fly = 1:obj.nflies
                nframes = numel(obj.walking_scores{fly});

                % Per-frame: true if neither classifier flags the frame
                % NaN == 1 is false, so NaN is treated as not-flagged
                obj.is_onfloor_perframe{fly} = ~(obj.onceiling_scores{fly} == 1 | obj.nottracking_scores{fly} == 1);

                % Detect walk bouts from walking_scores
                [walk_t0s, walk_t1s] = detect_bouts(obj.walking_scores{fly});

                % walkonfloor_digital: true only for frames in walks with frac_onfloor >= threshold
                obj.walkonfloor_digital{fly} = false(1, nframes);

                for w = 1:numel(walk_t0s)
                    t0 = walk_t0s(w);
                    t1 = min(walk_t1s(w), nframes);
                    if t0 < 1, continue; end

                    frac_onfloor = mean(obj.is_onfloor_perframe{fly}(t0:t1));
                    total_walks = total_walks + 1;

                    if frac_onfloor >= obj.frac_onfloor_threshold
                        obj.walkonfloor_digital{fly}(t0:t1) = true;
                        valid_walks = valid_walks + 1;
                    end
                end
            end

            fprintf('  %d / %d walks pass onfloor threshold (%.1f%%)\n', ...
                valid_walks, total_walks, 100 * valid_walks / max(1, total_walks));
        end


        function restrictBoutsbyCondition(obj,condition_name,boutdata,digital_signal_traj)
            % Restrict bout data using trajectory format digital signal
            % uses restrictedLimbBoutData - restrict to 'bouts' fully
            % contained within the digital signal bouts. 

            % Usage: analyzer.restrictBoutDataByCondition('stimON', digitalindicator_traj)

            % Apply restriction in trajectory format
            restricted_data_traj = restrictedLimbBoutData(digital_signal_traj, boutdata, 'trajectory');

            % Store trajectory format
            obj.restrictedBoutData([condition_name '_traj']) = restricted_data_traj;

        end

        function computeBoutMetricsByCondition(obj, condition_name)
            % Compute bout metrics using your computeboutmetrics2 function

            if ~obj.restrictedBoutData.isKey(condition_name)
                error('Condition "%s" not found. Run restrictBoutDataByCondition first.', condition_name);
            end

            boutstruct = obj.restrictedBoutData(condition_name);
            bout_metrics = computeboutmetrics2(obj.trx, obj.aptdata, obj.tips_pos_body, obj.legtip_landmarknums, boutstruct);
            obj.boutMetrics(condition_name) = bout_metrics;
        end

        function computeWalkMetricsforValidFrames(obj, validframes_name)

            digital_signal = obj.validFrames(validframes_name);
            walk_metrics = computeWalkMetrics(obj,digital_signal,obj.phase_methods);
            obj.walkMetrics(validframes_name) = walk_metrics;
        end
  
        function analyzeBoutAndStimConditions(obj)
            % Analyze combinations of walking and stimulus conditions
            % Uses both stored walking_scores and digitalindicator (both trajectory format)
            
            fprintf('Analyzing bout + stimulus conditions (trajectory format)...\n');
            
            % Restrict to walking periods first
            if obj.debug
            fprintf('Restricting to walking periods...\n');
            end

            % Now restrict limbBoutData data by walking scores as condition
            restrictBoutsbyCondition(obj,'walkingBoutData',obj.limbBoutData,obj.walking_scores)
            
            % Now restrict walking data by stimulus conditions
            if obj.debug
            fprintf('Applying stimulus restrictions to walking data...\n');
            end            
            restrictBoutsbyCondition(obj,'walking_stimON',obj.restrictedBoutData('walkingBoutData_traj'),obj.digitalindicator)            
            
            % stim OFF
            stim_off_indicator = obj.invertTrajectoryIndicator(obj.digitalindicator);
            % stimOFF_walkingData = restrictedLimbBoutData(stim_off_indicator, walkingBoutData, 'trajectory');
            restrictBoutsbyCondition(obj,'walking_stimOFF',obj.restrictedBoutData('walkingBoutData_traj'),stim_off_indicator)
            
            % Compute metrics
            if obj.debug
            fprintf('Computing bout metrics for walking conditions...\n');
            end
            obj.computeBoutMetricsByCondition('walking_stimON_traj');
            obj.computeBoutMetricsByCondition('walking_stimOFF_traj');
            
            if obj.debug
            fprintf('Walking + stimulus analysis complete!\n');
            end
        end


        function analyzeBoutAndStimConditions_onfloor(obj)
            % Analyze bout + stimulus conditions, restricted to onfloor walks.
            % Same restriction chain as analyzeBoutAndStimConditions but with
            % walkonfloor_digital inserted: walking → onfloor → LED.
            % Stores results under distinct keys with '_onfloor' suffix.

            if isempty(obj.walkonfloor_digital)
                obj.computeOnfloorSignal();
            end

            fprintf('Analyzing bout + stimulus conditions (onfloor filtered)...\n');

            % Restrict to walking periods first (reuses existing walkingBoutData if available)
            if ~obj.restrictedBoutData.isKey('walkingBoutData_traj')
                restrictBoutsbyCondition(obj, 'walkingBoutData', obj.limbBoutData, obj.walking_scores);
            end

            % Restrict walking bouts to onfloor walks
            restrictBoutsbyCondition(obj, 'walkingOnfloorBoutData', obj.restrictedBoutData('walkingBoutData_traj'), obj.walkonfloor_digital);

            % Restrict onfloor walking bouts by stimulus conditions
            restrictBoutsbyCondition(obj, 'walking_stimON_onfloor', obj.restrictedBoutData('walkingOnfloorBoutData_traj'), obj.digitalindicator);

            stim_off_indicator = obj.invertTrajectoryIndicator(obj.digitalindicator);
            restrictBoutsbyCondition(obj, 'walking_stimOFF_onfloor', obj.restrictedBoutData('walkingOnfloorBoutData_traj'), stim_off_indicator);

            % Compute bout metrics
            obj.computeBoutMetricsByCondition('walking_stimON_onfloor_traj');
            obj.computeBoutMetricsByCondition('walking_stimOFF_onfloor_traj');
        end


        function analyzeWalkAndStimConditions(obj)
            if obj.debug
            fprintf('Analyzing walk + stimulus conditions ...\n');
            end 
            % create valid frame signals for LED on and LED off
            obj.createValidFrames_LED();           

            % Compute metrics
            obj.computeWalkMetricsforValidFrames('led_on_traj');
            obj.computeWalkMetricsforValidFrames('led_off_traj');

        end


        function analyzeWalkAndStimConditions_onfloor(obj)
            % Compute walk metrics for LED on/off, restricted to onfloor walks.
            % Uses walkonfloor_digital AND'd with digitalindicator.
            % Stores results under distinct keys: 'led_on_onfloor_traj', 'led_off_onfloor_traj'

            if isempty(obj.walkonfloor_digital)
                obj.computeOnfloorSignal();
            end

            fprintf('Analyzing walk + stimulus conditions (onfloor filtered)...\n');

            % AND digitalindicator with walkonfloor_digital for each LED condition
            led_on_onfloor = cell(obj.nflies, 1);
            led_off_onfloor = cell(obj.nflies, 1);
            stim_off_indicator = obj.invertTrajectoryIndicator(obj.digitalindicator);
            for fly = 1:obj.nflies
                led_on_onfloor{fly} = obj.digitalindicator{fly} & obj.walkonfloor_digital{fly};
                led_off_onfloor{fly} = stim_off_indicator{fly} & obj.walkonfloor_digital{fly};
            end

            obj.validFrames('led_on_onfloor_traj') = led_on_onfloor;
            obj.validFrames('led_off_onfloor_traj') = led_off_onfloor;

            % Compute walk metrics on filtered valid frames
            obj.computeWalkMetricsforValidFrames('led_on_onfloor_traj');
            obj.computeWalkMetricsforValidFrames('led_off_onfloor_traj');
        end


        function saveBoutData(obj, condition_names)
            % Save bout data for specified conditions
            if nargin < 2
                condition_names = keys(obj.restrictedBoutData);
            end

            if ~isempty(obj.expdir)
                save_data = struct();
                for i = 1:length(condition_names)
                    condition = condition_names{i};
                    if obj.restrictedBoutData.isKey(condition)
                        % Use variable names matching your original code
                        if strcmp(condition, 'stimON')
                            save_data.stimON_walkingLimbBoutData = obj.restrictedBoutData(condition);
                        elseif strcmp(condition, 'stimOFF')
                            save_data.stimOFF_walkingLimbBoutData = obj.restrictedBoutData(condition);
                        else
                            save_data.(condition) = obj.restrictedBoutData(condition);
                        end
                    end
                end

                save(fullfile(obj.expdir, 'swingtstancestep.mat'), '-struct', 'save_data');
                fprintf('Bout data saved to: %s\n', fullfile(obj.expdir, 'swingtstancestep.mat'));
            end
        end
        
        function saveResults(obj, filename)
            % Save analysis results in your format
            if nargin < 2
                filename = 'locomotionmetricsswingstanceboutstats.mat';
            end
            if ~isempty(obj.expdir)
                filepath = fullfile(obj.expdir, filename);
            else
                filepath = filename;
            end

            % Save in format matching your original code

            bout_metrics_ON = [];
            bout_metrics_OFF = [];

            % TO DO cleanup condition names everywhere at some point
            if obj.boutMetrics.isKey('walking_stimON_traj')
                bout_metrics_ON = obj.boutMetrics('walking_stimON_traj');
            end
            if obj.boutMetrics.isKey('walking_stimOFF_traj')
                bout_metrics_OFF = obj.boutMetrics('walking_stimOFF_traj');
            end
            % add walk metrics
            if obj.walkMetrics.isKey('led_on_traj')
                walk_metrics_ON = obj.walkMetrics('led_on_traj');
            end
            if obj.walkMetrics.isKey('led_off_traj')
                walk_metrics_OFF = obj.walkMetrics('led_off_traj');
            end
            save(filepath, 'bout_metrics_ON', 'bout_metrics_OFF','walk_metrics_ON','walk_metrics_OFF');
            fprintf('Results saved to: %s\n', filepath);
        end

        function locostatsperexp = computeStatsPerExp(obj, conditions, key_suffix)
            % Flatten bout_metrics and walk_metrics into a single per-experiment
            % summary struct. LED condition is encoded in each field name
            % (e.g. feature__swing__LEDon__all, feature__walk__LEDoff__all).
            %
            % Usage:
            %   stats = analyzer.computeStatsPerExp();                        % both ON and OFF, unfiltered
            %   stats = analyzer.computeStatsPerExp({'ON'});                  % ON only
            %   stats = analyzer.computeStatsPerExp({'ON','OFF'}, '_onfloor'); % onfloor-filtered keys
            %
            % key_suffix appended to map keys (e.g. '_onfloor' looks up
            % 'walking_stimON_onfloor_traj' and 'led_on_onfloor_traj').
            %
            % Based on LocomotionCombinePerFrameStats.m (Alice)

            if nargin < 2
                conditions = {'ON', 'OFF'};
            end
            if nargin < 3
                key_suffix = '';
            end

            % Map condition labels to stored keys and LED field name labels
            condMap = struct();
            condMap.ON.boutKey  = ['walking_stimON' key_suffix '_traj'];
            condMap.ON.walkKey  = ['led_on' key_suffix '_traj'];
            condMap.ON.label    = 'LEDon';
            condMap.OFF.boutKey = ['walking_stimOFF' key_suffix '_traj'];
            condMap.OFF.walkKey = ['led_off' key_suffix '_traj'];
            condMap.OFF.label   = 'LEDoff';

            statsperexp = struct;

            for c = 1:numel(conditions)
                cond = conditions{c};
                if ~isfield(condMap, cond)
                    error('Unknown condition "%s". Use "ON" or "OFF".', cond);
                end
                cm = condMap.(cond);

                % Bout metrics (swing/stance scalars, per-frame structs, velmag bins)
                if obj.boutMetrics.isKey(cm.boutKey)
                    bout_metrics = obj.boutMetrics(cm.boutKey);
                    statsperexp = obj.combineBoutMetrics(bout_metrics, cm.label, statsperexp);
                    statsperexp = obj.combineStepMetrics(bout_metrics, cm.label, statsperexp);
                else
                    fprintf('Warning: bout metrics key "%s" not found, skipping bout/step stats for %s\n', cm.boutKey, cond);
                end

                % Walk metrics (per-frame features and phase)
                if obj.walkMetrics.isKey(cm.walkKey)
                    walk_metrics = obj.walkMetrics(cm.walkKey);
                    statsperexp = obj.combineWalkMetrics(walk_metrics, cm.label, statsperexp);
                    statsperexp = obj.combinePhaseMetrics(walk_metrics, cm.label, statsperexp);
                else
                    fprintf('Warning: walk metrics key "%s" not found, skipping walk/phase stats for %s\n', cm.walkKey, cond);
                end
            end

            obj.locoStatsPerExp = statsperexp;
            locostatsperexp = statsperexp;
            fprintf('Per-experiment stats computed for conditions: %s\n', strjoin(conditions, ', '));
        end

        function saveStatsPerExp(obj, filename)
            % Save per-experiment stats to .mat file
            % Usage: analyzer.saveStatsPerExp('locostatsperexp.mat')

            if nargin < 2
                filename = 'locostatsperexp.mat';
            end
            if ~isempty(obj.expdir)
                filepath = fullfile(obj.expdir, filename);
            else
                filepath = filename;
            end

            locostatsperexp = obj.locoStatsPerExp; %#ok<PROP>
            save(filepath, 'locostatsperexp');
            fprintf('Per-experiment stats saved to: %s\n', filepath);
        end

        function buildWalkStruct(obj, conditions, key_suffix)
            % Build flat per-walk and per-step struct-of-arrays from in-memory
            % walk_metrics and bout_metrics. Same logic as build_walk_struct.m
            % but operates on data already held in the analyzer.
            %
            % Usage:
            %   analyzer.buildWalkStruct();                         % both ON and OFF, unfiltered
            %   analyzer.buildWalkStruct({'OFF'});                  % OFF only
            %   analyzer.buildWalkStruct({'ON','OFF'}, '_onfloor'); % onfloor-filtered keys

            if nargin < 2
                conditions = {'ON', 'OFF'};
            end
            if nargin < 3
                key_suffix = '';
            end

            condMap = struct();
            condMap.ON.boutKey  = ['walking_stimON' key_suffix '_traj'];
            condMap.ON.walkKey  = ['led_on' key_suffix '_traj'];
            condMap.OFF.boutKey = ['walking_stimOFF' key_suffix '_traj'];
            condMap.OFF.walkKey = ['led_off' key_suffix '_traj'];

            for c = 1:numel(conditions)
                cond = conditions{c};
                if ~isfield(condMap, cond)
                    error('Unknown condition "%s". Use "ON" or "OFF".', cond);
                end
                cm = condMap.(cond);

                if ~obj.boutMetrics.isKey(cm.boutKey)
                    error('Bout metrics key "%s" not found. Run analyzeBoutAndStimConditions first.', cm.boutKey);
                end
                if ~obj.walkMetrics.isKey(cm.walkKey)
                    error('Walk metrics key "%s" not found. Run analyzeWalkAndStimConditions first.', cm.walkKey);
                end

                bout_metrics = obj.boutMetrics(cm.boutKey);
                walk_metrics = obj.walkMetrics(cm.walkKey);

                [walk_struct, step_struct] = obj.buildWalkStructForCondition(walk_metrics, bout_metrics);
                obj.walkStruct.(cond).walk_struct = walk_struct;
                obj.walkStruct.(cond).step_struct = step_struct;
            end

            fprintf('Walk struct built for conditions: %s\n', strjoin(conditions, ', '));
        end

        function saveWalkStruct(obj, filename)
            % Save walk struct to .mat file with variables named
            % walk_struct_OFF, step_struct_OFF, etc.
            %
            % Usage: analyzer.saveWalkStruct('locomotion_walkstruct.mat')

            if nargin < 2
                filename = 'locomotion_walkstruct.mat';
            end
            if ~isempty(obj.expdir)
                filepath = fullfile(obj.expdir, filename);
            else
                filepath = filename;
            end

            conds = fieldnames(obj.walkStruct);
            save_data = struct();
            for c = 1:numel(conds)
                cond = conds{c};
                save_data.(['walk_struct_' cond]) = obj.walkStruct.(cond).walk_struct;
                save_data.(['step_struct_' cond]) = obj.walkStruct.(cond).step_struct;
            end

            save(filepath, '-struct', 'save_data');
            fprintf('Walk struct saved to: %s\n', filepath);
        end

        function analyzeCustomCondition(obj, condition_name, custom_digital_signal_traj)
            % Claude addition - haven't tested. Might be useful for other
            % ways to chop up data by stim. 
            % Analyze data with custom digital indicator (trajectory format)
            % Usage: analyzer.analyzeCustomCondition('custom_condition', custom_signal_traj, true)
            %   where custom_signal_traj is cell array: custom_signal_traj{fly}

            fprintf('Analyzing condition: %s\n', condition_name);
            obj.restrictBoutsbyCondition(condition_name, custom_digital_signal_traj);            

        end

        function data = getConditionData(obj, condition_name, data_type)
            % Get data for a specific condition
            % Usage: data = analyzer.getConditionData('walking_stimON', 'restricted') % or 'metrics'
            
            if strcmp(data_type, 'restricted')
                if obj.restrictedBoutData.isKey(condition_name)
                    data = obj.restrictedBoutData(condition_name);
                else
                    error('Condition "%s" not found in restricted bout data', condition_name);
                end
            elseif strcmp(data_type, 'metrics')
                if obj.boutMetrics.isKey(condition_name)
                    data = obj.boutMetrics(condition_name);
                else
                    error('Condition "%s" not found in bout metrics', condition_name);
                end
            else
                error('data_type must be "restricted" or "metrics"');
            end
        end

        function boutmetrics_list = listBoutMetrics(obj)
            % List all analyzed conditions
            boutmetrics_list = keys(obj.boutMetrics);
            fprintf('Bout Metrics for contions:\n');
            for i = 1:length(boutmetrics_list)
                fprintf('  %s\n', boutmetrics_list{i});
            end
        end

        function conditions_list = listRestrictionConditions(obj)
            % List all analyzed conditions
            conditions_list = keys(obj.restrictedBoutData);
            fprintf('Restricted conditions:\n');
            for i = 1:length(conditions_list)
                fprintf('  %s\n', conditions_list{i});
            end
        end


        function summary = getConditionSummary(obj, condition_name)
            % Get summary statistics for a specific condition
            if ~obj.boutMetrics.isKey(condition_name)
                error('Condition "%s" not found', condition_name);
            end
            
            bout_metrics = obj.boutMetrics(condition_name);
            summary = struct();
            summary.condition = condition_name;
            summary.nflies = obj.nflies;
            
            for state_idx = 1:2  % swing and stance only
                state_name = obj.states{state_idx};
                
                if isfield(bout_metrics.allflies.all_limbs, state_name) && ...
                   isfield(bout_metrics.allflies.all_limbs.(state_name), 'durations_time')
                    
                    data = bout_metrics.allflies.all_limbs.(state_name).durations_time;
                    summary.(state_name).mean = mean(data);
                    summary.(state_name).std = std(data);
                    summary.(state_name).count = length(data);
                    summary.(state_name).median = median(data);
                end
            end
        end        

        function displayDataStructure(obj)
            % Display information about the loaded data structure
            fprintf('\n=== LimbBoutAnalyzer Data Structure ===\n');
            fprintf('Number of flies: %d\n', obj.nflies);
            fprintf('Data format: trajectory (all data per-fly)\n');

            % Show per-fly stimulus and walking info
            fprintf('\nPer-fly data summary:\n');
            for fly = 1:min(3, obj.nflies)  % Show first 3 flies
                stim_frames = sum(obj.digitalindicator{fly});
                walking_frames = sum(obj.walking_scores{fly});
                total_frames = length(obj.digitalindicator{fly});

                fprintf('  Fly %d: %d frames, %d stim (%.1f%%), %d walking (%.1f%%)\n', ...
                    fly, total_frames, stim_frames, 100*stim_frames/total_frames, ...
                    walking_frames, 100*walking_frames/total_frames);
            end
            if obj.nflies > 3
                fprintf('  ... and %d more flies\n', obj.nflies - 3);
            end

            % Show bout data structure
            if ~isempty(obj.limbBoutData)
                fprintf('\nBout data structure:\n');
                fprintf('  Flies: %d\n', length(obj.limbBoutData));
                if length(obj.limbBoutData) >= 1
                    fprintf('  Limbs per fly: %d\n', length(obj.limbBoutData(1).perlimb));
                    if length(obj.limbBoutData(1).perlimb) >= 1
                        states_available = fieldnames(obj.limbBoutData(1).perlimb(1));
                        fprintf('  States available: %s\n', strjoin(states_available, ', '));
                    end
                end
            end

            % Show restricted bout data conditions
            fprintf('\nRestricted bout data conditions:\n');
            if obj.restrictedBoutData.Count > 0
                restricted_conditions = keys(obj.restrictedBoutData);
                fprintf('  Number of conditions: %d\n', obj.restrictedBoutData.Count);
                fprintf('  Conditions: %s\n', strjoin(restricted_conditions, ', '));
            else
                fprintf('  No conditions analyzed yet\n');
            end

            % Show bout metrics conditions
            fprintf('\nBout metrics conditions:\n');
            if obj.boutMetrics.Count > 0
                metrics_conditions = keys(obj.boutMetrics);
                fprintf('  Number of conditions: %d\n', obj.boutMetrics.Count);
                fprintf('  Conditions: %s\n', strjoin(metrics_conditions, ', '));
            else
                fprintf('  No metrics computed yet\n');
            end

            % Show analysis suggestions
            fprintf('\nAnalysis suggestions:\n');
            if obj.restrictedBoutData.Count == 0 && obj.boutMetrics.Count == 0
                fprintf('  Run analyzer.analyzeBoutAndStimConditions() for walking + stimulus analysis (with metrics)\n');
                fprintf('  Run analyzer.analyzeStimConditions() for basic stimulus restriction (no metrics)\n');
            elseif obj.restrictedBoutData.Count > 0 && obj.boutMetrics.Count == 0
                fprintf('  Run analyzer.computeBoutMetricsByCondition(condition_name) to compute metrics\n');
                fprintf('  Note: Only walking+stimulus conditions typically need metrics\n');
            else
                fprintf('  Run analyzer.listConditions() to see all available conditions\n');
                fprintf('  Run analyzer.plotConditionComparison(cond1, cond2) to compare conditions\n');
                fprintf('  Run analyzer.getConditionSummary(condition_name) for detailed statistics\n');
            end

            fprintf('=======================================\n\n');
        end

    end

    methods (Access = private)
        function digitalindicator_traj = convertDigitalIndicatorToTrajectory(obj, digitalindicator)

            % Convert digitalindicator to trajectory format (auto-detects input format)

            if iscell(digitalindicator)
                % Already in trajectory format
                if length(digitalindicator) == obj.nflies
                    digitalindicator_traj = digitalindicator;
                    fprintf('Digital indicator already in trajectory format\n');
                else
                    error('Digital indicator cell array must have %d flies', obj.nflies);
                end
            elseif isvector(digitalindicator)
                % Movie format - convert to trajectory format
                fprintf('Converting digital indicator from movie to trajectory format...\n');
                digitalindicator_traj = cell(obj.nflies, 1);

                for fly = 1:obj.nflies
                    % Get this fly's frames in the movie using your method
                    fly_frames = obj.getFlyFramesInMovie(fly);

                    % Extract this fly's stimulus data
                    digitalindicator_traj{fly} = digitalindicator(fly_frames);
                end
                fprintf('Conversion complete\n');
            else
                error('digitalindicator must be either vector (movie format) or cell array (trajectory format)');
            end
        end

        function fly_frames = getFlyFramesInMovie(obj, fly)
            % Get the movie frames that correspond to this fly's trajectory
            % This depends on how trx stores the mapping between trajectory and movie frames
            fly_frames = obj.trx.firstframes(fly):obj.trx.endframes(fly);
        end

        function inverted_indicator = invertTrajectoryIndicator(obj, trajectory_indicator)
            % Invert trajectory format indicator for all flies
            inverted_indicator = cell(obj.nflies, 1);
            for fly = 1:obj.nflies
                inverted_indicator{fly} = ~trajectory_indicator{fly};
            end
        end

        function limbname = getLimbName(~, numlimb, ll)
            % Get limb name string based on number of limbs in grouping
            if numlimb == 1
                limbname = 'all';
            elseif numlimb == 3
                limbname = sprintf('pair%s', num2str(ll));
            elseif numlimb == 6
                limbname = sprintf('limb%s', num2str(ll));
            end
        end

        function [walk_struct, step_struct] = buildWalkStructForCondition(obj, walk_metrics, bout_metrics) %#ok<INUSL>
            % Build flat per-walk and per-step struct-of-arrays from
            % walk_metrics and bout_metrics held in memory.
            % Core logic ported from build_walk_struct.m.

            pw = walk_metrics.perwalk;
            nwalks = numel(pw);
            nlimb = 6;
            legnames = {'RF','RM','RH','LH','LM','LF'};

            %% Field name definitions
            movement_fields = {'velmag_ctr','absdv_ctr','absdu_ctr','absdtheta', ...
                'forward_vel','backward_vel','left_vel','right_vel', ...
                'left_dtheta','right_dtheta','CoM_stability'};

            phase_pairs = {'LM_LH','RM_RH','LF_LM','RF_RM','RM_LH','LM_RH','LF_RM','RF_LM'};
            abs_phase_pairs = {'absRF_LF','absRM_LM','absRH_LH'};
            abs_phase_outnames = {'absphase_RF_LF','absphase_RM_LM','absphase_RH_LH'};
            phase_groups = {'ipsi_post_2','ipsi_ant_2','tripods_4','ipsi_P2A_4'};

            phase_group_pairs = struct();
            phase_group_pairs.ipsi_post_2 = {{'LM_LH','RM_RH'}};
            phase_group_pairs.ipsi_ant_2 = {{'LF_LM','RF_RM'}};
            phase_group_pairs.tripods_4 = {{'RF_LM','LM_RH','LF_RM','RM_LH'}};
            phase_group_pairs.ipsi_P2A_4 = {{'LM_LH','RM_RH','LF_LM','RF_RM'}};
            abs_group_pairs = {'absRF_LF','absRM_LM','absRH_LH'};

            step_geom_2row = {'AEP','AEP_BL','PEP','PEP_BL'};
            step_geom_scalar = {'amplitude_BL','amplitude_px','distance_BL','distance_px', ...
                'length_BL','length_px','step_direction','speed_BLpers','speed_pxpers'};

            tipspeed_bodyref = {'mean_tips_speed_bodyref','std_tips_speed_bodyref','min_speed_bodyref','max_speed_bodyref'};
            tipspeed_globalref = {'mean_tips_speed_globalref','std_tips_speed_globalref','min_speed_globalref','max_speed_globalref'};
            tipspeed_out_bodyref = {'mean','std','min','max'};
            tipspeed_out_globalref = {'mean','std','min','max'};

            %% ===== BUILD WALK_STRUCT =====

            walk_struct = struct();

            % --- Metadata ---
            walk_struct.fly = [pw.fly];
            walk_struct.walk_t0 = [pw.walk_t0];
            walk_struct.walk_t1 = [pw.walk_t1];
            walk_struct.walk_duration = walk_struct.walk_t1 - walk_struct.walk_t0 + 1;
            walk_struct.nframes = walk_struct.walk_duration;

            % --- Onfloor filtering ---
            % Compute frac_onfloor per walk from is_onfloor_perframe
            walk_struct.frac_onfloor = nan(1, nwalks);
            walk_struct.walk_valid = true(1, nwalks);  % default: all valid
            if ~isempty(obj.is_onfloor_perframe)
                for w = 1:nwalks
                    fly = pw(w).fly;
                    t0 = pw(w).walk_t0;
                    t1 = min(pw(w).walk_t1, numel(obj.is_onfloor_perframe{fly}));
                    if t0 >= 1 && t0 <= numel(obj.is_onfloor_perframe{fly})
                        walk_struct.frac_onfloor(w) = mean(obj.is_onfloor_perframe{fly}(t0:t1));
                    end
                end
                walk_struct.walk_valid = walk_struct.frac_onfloor >= obj.frac_onfloor_threshold;
            end

            % --- Movement means ---
            for m = 1:numel(movement_fields)
                fn = movement_fields{m};
                vals = nan(1, nwalks);
                for w = 1:nwalks
                    vals(w) = pw(w).(fn).mean;
                end
                walk_struct.(fn) = vals;
            end

            % --- Phase groups (circular) ---
            for g = 1:numel(phase_groups)
                gname = phase_groups{g};
                vals = nan(1, nwalks);
                for w = 1:nwalks
                    vals(w) = pw(w).phasediff_hilbert.(gname).mean;
                end
                walk_struct.(['phase_' gname]) = vals;
            end

            % --- Abs phase group (non-circular) ---
            vals = nan(1, nwalks);
            for w = 1:nwalks
                vals(w) = pw(w).phasediff_hilbert.abscontra_L2R_3.mean;
            end
            walk_struct.absphase_contra_L2R_3 = vals;

            % --- Phase pairs (circular) ---
            for p = 1:numel(phase_pairs)
                pname = phase_pairs{p};
                vals = nan(1, nwalks);
                for w = 1:nwalks
                    vals(w) = pw(w).phasediff_hilbert.(pname).mean;
                end
                walk_struct.(['phase_' pname]) = vals;
            end

            % --- Abs phase pairs (non-circular) ---
            for p = 1:numel(abs_phase_pairs)
                pname = abs_phase_pairs{p};
                vals = nan(1, nwalks);
                for w = 1:nwalks
                    vals(w) = pw(w).phasediff_hilbert.(pname).mean;
                end
                walk_struct.(abs_phase_outnames{p}) = vals;
            end

            % --- Phase trim per limb ---
            for limb = 1:nlimb
                ln = legnames{limb};
                sidx_vals = nan(1, nwalks);
                eidx_vals = nan(1, nwalks);
                for w = 1:nwalks
                    pd = pw(w).phasediff_hilbert.phasedata(limb,:);
                    valid = find(~isnan(pd));
                    if ~isempty(valid)
                        sidx_vals(w) = valid(1);
                        eidx_vals(w) = valid(end);
                    end
                end
                walk_struct.(['phase_sidx_' ln]) = sidx_vals;
                walk_struct.(['phase_eidx_' ln]) = eidx_vals;
            end

            % --- Initialize per-limb walk fields ---
            for limb = 1:nlimb
                ln = legnames{limb};
                for f = 1:numel(step_geom_2row)
                    fname = step_geom_2row{f};
                    walk_struct.(['step_' fname 'x_' ln]) = nan(1, nwalks);
                    walk_struct.(['step_' fname 'y_' ln]) = nan(1, nwalks);
                end
                for f = 1:numel(step_geom_scalar)
                    walk_struct.(['step_' step_geom_scalar{f} '_' ln]) = nan(1, nwalks);
                end
                walk_struct.(['step_instataeous_frequency_steps_' ln]) = nan(1, nwalks);
                walk_struct.(['step_duration_' ln]) = nan(1, nwalks);
                walk_struct.(['swing_duration_' ln]) = nan(1, nwalks);
                walk_struct.(['stance_duration_' ln]) = nan(1, nwalks);
                for stat = tipspeed_out_bodyref
                    walk_struct.(['swing_tips_speed_bodyref_' stat{1} '_' ln]) = nan(1, nwalks);
                    walk_struct.(['stance_tips_speed_bodyref_' stat{1} '_' ln]) = nan(1, nwalks);
                end
                for stat = tipspeed_out_globalref
                    walk_struct.(['swing_tips_speed_globalref_' stat{1} '_' ln]) = nan(1, nwalks);
                    walk_struct.(['stance_tips_speed_globalref_' stat{1} '_' ln]) = nan(1, nwalks);
                end
                walk_struct.(['step_nsteps_' ln]) = zeros(1, nwalks);
            end

            %% ===== BUILD INDEX MAPS AND COUNT STEPS =====

            nflies = numel(bout_metrics.perfly);
            stance_map = cell(nflies, nlimb);
            swing_map = cell(nflies, nlimb);
            for fly = 1:nflies
                for limb = 1:nlimb
                    stance_starts = bout_metrics.perfly(fly).perlimb(limb).stance.start_indices;
                    swing_ends = bout_metrics.perfly(fly).perlimb(limb).swing.end_indices;
                    sm = containers.Map('KeyType','int64','ValueType','int64');
                    for k = 1:numel(stance_starts)
                        sm(int64(stance_starts(k))) = int64(k);
                    end
                    stance_map{fly, limb} = sm;
                    wm = containers.Map('KeyType','int64','ValueType','int64');
                    for k = 1:numel(swing_ends)
                        wm(int64(swing_ends(k))) = int64(k);
                    end
                    swing_map{fly, limb} = wm;
                end
            end

            total_steps = 0;
            for w = 1:nwalks
                fly = pw(w).fly;
                wt0 = pw(w).walk_t0;
                wt1 = pw(w).walk_t1;
                for limb = 1:nlimb
                    step_data = bout_metrics.perfly(fly).perlimb(limb).step;
                    step_t0s = step_data.stepfeatures.start_indices;
                    step_t1s = step_data.stepfeatures.end_indices;
                    total_steps = total_steps + sum(step_t0s >= wt0 & step_t1s <= wt1);
                end
            end

            %% ===== PRE-ALLOCATE STEP_STRUCT =====

            step_struct = struct();
            step_struct.walk_idx = zeros(1, total_steps);
            step_struct.fly = zeros(1, total_steps);
            step_struct.limb = zeros(1, total_steps);
            step_struct.step_t0 = nan(1, total_steps);
            step_struct.step_t1 = nan(1, total_steps);

            for f = 1:numel(step_geom_2row)
                fname = step_geom_2row{f};
                step_struct.([fname 'x']) = nan(1, total_steps);
                step_struct.([fname 'y']) = nan(1, total_steps);
            end
            for f = 1:numel(step_geom_scalar)
                step_struct.(step_geom_scalar{f}) = nan(1, total_steps);
            end

            step_struct.step_duration = nan(1, total_steps);
            step_struct.instataeous_frequency_steps = nan(1, total_steps);

            step_struct.stance_t0 = nan(1, total_steps);
            step_struct.stance_t1 = nan(1, total_steps);
            step_struct.stance_duration = nan(1, total_steps);
            step_struct.swing_t0 = nan(1, total_steps);
            step_struct.swing_t1 = nan(1, total_steps);
            step_struct.swing_duration = nan(1, total_steps);

            for stat = tipspeed_out_bodyref
                step_struct.(['swing_tips_speed_bodyref_' stat{1}]) = nan(1, total_steps);
                step_struct.(['stance_tips_speed_bodyref_' stat{1}]) = nan(1, total_steps);
            end
            for stat = tipspeed_out_globalref
                step_struct.(['swing_tips_speed_globalref_' stat{1}]) = nan(1, total_steps);
                step_struct.(['stance_tips_speed_globalref_' stat{1}]) = nan(1, total_steps);
            end

            for m = 1:numel(movement_fields)
                step_struct.(movement_fields{m}) = nan(1, total_steps);
            end

            for p = 1:numel(phase_pairs)
                step_struct.(['phase_' phase_pairs{p}]) = nan(1, total_steps);
            end
            for g = 1:numel(phase_groups)
                step_struct.(['phase_' phase_groups{g}]) = nan(1, total_steps);
            end
            for p = 1:numel(abs_phase_outnames)
                step_struct.(abs_phase_outnames{p}) = nan(1, total_steps);
            end
            step_struct.absphase_contra_L2R_3 = nan(1, total_steps);
            step_struct.phase_valid = false(1, total_steps);

            %% ===== MAIN LOOP: POPULATE WALK PER-LIMB FIELDS AND STEP_STRUCT =====

            step_ct = 0;

            for w = 1:nwalks
                fly = pw(w).fly;
                wt0 = pw(w).walk_t0;
                wt1 = pw(w).walk_t1;

                for limb = 1:nlimb
                    ln = legnames{limb};

                    step_data = bout_metrics.perfly(fly).perlimb(limb).step;
                    swing_data = bout_metrics.perfly(fly).perlimb(limb).swing;
                    stance_data = bout_metrics.perfly(fly).perlimb(limb).stance;

                    step_t0s = step_data.stepfeatures.start_indices;
                    step_t1s = step_data.stepfeatures.end_indices;

                    match_mask = step_t0s >= wt0 & step_t1s <= wt1;
                    match_idx = find(match_mask);

                    if isempty(match_idx)
                        continue;
                    end

                    nmatched = numel(match_idx);
                    walk_struct.(['step_nsteps_' ln])(w) = nmatched;

                    % --- Step geometry 2-row fields ---
                    for f = 1:numel(step_geom_2row)
                        fname = step_geom_2row{f};
                        geom_vals_x = nan(1, nmatched);
                        geom_vals_y = nan(1, nmatched);
                        for si = 1:nmatched
                            i = match_idx(si);
                            st_key = int64(step_t0s(i));
                            if stance_map{fly, limb}.isKey(st_key)
                                sidx = stance_map{fly, limb}(st_key);
                                geom_vals_x(si) = step_data.stepfeatures.(fname)(1, sidx);
                                geom_vals_y(si) = step_data.stepfeatures.(fname)(2, sidx);
                            end
                        end
                        walk_struct.(['step_' fname 'x_' ln])(w) = mean(geom_vals_x, 'omitnan');
                        walk_struct.(['step_' fname 'y_' ln])(w) = mean(geom_vals_y, 'omitnan');
                    end

                    % --- Step geometry scalar fields ---
                    for f = 1:numel(step_geom_scalar)
                        fn = step_geom_scalar{f};
                        walk_struct.(['step_' fn '_' ln])(w) = mean(step_data.stepfeatures.(fn)(match_idx), 'omitnan');
                    end
                    walk_struct.(['step_instataeous_frequency_steps_' ln])(w) = mean(step_data.stepfeatures.instataeous_frequency_steps(match_idx), 'omitnan');

                    % --- Step duration ---
                    walk_struct.(['step_duration_' ln])(w) = mean(step_data.stepfeatures.durations_time(match_idx), 'omitnan');

                    % --- Stance and swing durations and tip speeds ---
                    stance_durs = nan(1, nmatched);
                    swing_durs = nan(1, nmatched);
                    swing_tipspeed = struct();
                    stance_tipspeed = struct();
                    for ts = 1:numel(tipspeed_bodyref)
                        swing_tipspeed.(tipspeed_bodyref{ts}) = nan(1, nmatched);
                        stance_tipspeed.(tipspeed_bodyref{ts}) = nan(1, nmatched);
                    end
                    for ts = 1:numel(tipspeed_globalref)
                        swing_tipspeed.(tipspeed_globalref{ts}) = nan(1, nmatched);
                        stance_tipspeed.(tipspeed_globalref{ts}) = nan(1, nmatched);
                    end

                    for si = 1:nmatched
                        i = match_idx(si);
                        st_key = int64(step_t0s(i));
                        if stance_map{fly, limb}.isKey(st_key)
                            sidx = stance_map{fly, limb}(st_key);
                            stance_durs(si) = stance_data.durations_time(sidx);
                            for ts = 1:numel(tipspeed_bodyref)
                                stance_tipspeed.(tipspeed_bodyref{ts})(si) = stance_data.(tipspeed_bodyref{ts})(sidx);
                            end
                            for ts = 1:numel(tipspeed_globalref)
                                stance_tipspeed.(tipspeed_globalref{ts})(si) = stance_data.(tipspeed_globalref{ts})(sidx);
                            end
                        end
                        sw_key = int64(step_t1s(i));
                        if swing_map{fly, limb}.isKey(sw_key)
                            swidx = swing_map{fly, limb}(sw_key);
                            swing_durs(si) = swing_data.durations_time(swidx);
                            for ts = 1:numel(tipspeed_bodyref)
                                swing_tipspeed.(tipspeed_bodyref{ts})(si) = swing_data.(tipspeed_bodyref{ts})(swidx);
                            end
                            for ts = 1:numel(tipspeed_globalref)
                                swing_tipspeed.(tipspeed_globalref{ts})(si) = swing_data.(tipspeed_globalref{ts})(swidx);
                            end
                        end
                    end

                    walk_struct.(['stance_duration_' ln])(w) = mean(stance_durs, 'omitnan');
                    walk_struct.(['swing_duration_' ln])(w) = mean(swing_durs, 'omitnan');

                    for ts = 1:numel(tipspeed_out_bodyref)
                        walk_struct.(['swing_tips_speed_bodyref_' tipspeed_out_bodyref{ts} '_' ln])(w) = mean(swing_tipspeed.(tipspeed_bodyref{ts}), 'omitnan');
                        walk_struct.(['stance_tips_speed_bodyref_' tipspeed_out_bodyref{ts} '_' ln])(w) = mean(stance_tipspeed.(tipspeed_bodyref{ts}), 'omitnan');
                    end
                    for ts = 1:numel(tipspeed_out_globalref)
                        walk_struct.(['swing_tips_speed_globalref_' tipspeed_out_globalref{ts} '_' ln])(w) = mean(swing_tipspeed.(tipspeed_globalref{ts}), 'omitnan');
                        walk_struct.(['stance_tips_speed_globalref_' tipspeed_out_globalref{ts} '_' ln])(w) = mean(stance_tipspeed.(tipspeed_globalref{ts}), 'omitnan');
                    end

                    % --- Populate step_struct for each matched step ---
                    for si = 1:nmatched
                        step_ct = step_ct + 1;
                        i = match_idx(si);

                        step_struct.walk_idx(step_ct) = w;
                        step_struct.fly(step_ct) = fly;
                        step_struct.limb(step_ct) = limb;
                        step_struct.step_t0(step_ct) = step_t0s(i);
                        step_struct.step_t1(step_ct) = step_t1s(i);

                        step_struct.step_duration(step_ct) = step_data.stepfeatures.durations_time(i);
                        step_struct.instataeous_frequency_steps(step_ct) = step_data.stepfeatures.instataeous_frequency_steps(i);

                        st_key = int64(step_t0s(i));
                        if stance_map{fly, limb}.isKey(st_key)
                            sidx = stance_map{fly, limb}(st_key);
                            for f = 1:numel(step_geom_2row)
                                fname = step_geom_2row{f};
                                step_struct.([fname 'x'])(step_ct) = step_data.stepfeatures.(fname)(1, sidx);
                                step_struct.([fname 'y'])(step_ct) = step_data.stepfeatures.(fname)(2, sidx);
                            end
                            step_struct.stance_t0(step_ct) = stance_data.start_indices(sidx);
                            step_struct.stance_t1(step_ct) = stance_data.end_indices(sidx);
                            step_struct.stance_duration(step_ct) = stance_data.durations_time(sidx);
                            for ts = 1:numel(tipspeed_bodyref)
                                step_struct.(['stance_tips_speed_bodyref_' tipspeed_out_bodyref{ts}])(step_ct) = stance_data.(tipspeed_bodyref{ts})(sidx);
                            end
                            for ts = 1:numel(tipspeed_globalref)
                                step_struct.(['stance_tips_speed_globalref_' tipspeed_out_globalref{ts}])(step_ct) = stance_data.(tipspeed_globalref{ts})(sidx);
                            end
                        end

                        for f = 1:numel(step_geom_scalar)
                            fn = step_geom_scalar{f};
                            step_struct.(fn)(step_ct) = step_data.stepfeatures.(fn)(i);
                        end

                        sw_key = int64(step_t1s(i));
                        if swing_map{fly, limb}.isKey(sw_key)
                            swidx = swing_map{fly, limb}(sw_key);
                            step_struct.swing_t0(step_ct) = swing_data.start_indices(swidx);
                            step_struct.swing_t1(step_ct) = swing_data.end_indices(swidx);
                            step_struct.swing_duration(step_ct) = swing_data.durations_time(swidx);
                            for ts = 1:numel(tipspeed_bodyref)
                                step_struct.(['swing_tips_speed_bodyref_' tipspeed_out_bodyref{ts}])(step_ct) = swing_data.(tipspeed_bodyref{ts})(swidx);
                            end
                            for ts = 1:numel(tipspeed_globalref)
                                step_struct.(['swing_tips_speed_globalref_' tipspeed_out_globalref{ts}])(step_ct) = swing_data.(tipspeed_globalref{ts})(swidx);
                            end
                        end

                        for m = 1:numel(movement_fields)
                            fn = movement_fields{m};
                            step_struct.(fn)(step_ct) = step_data.(fn).mean(i);
                        end

                        % Phase during step
                        rel_t0 = step_t0s(i) - wt0 + 1;
                        rel_t1 = step_t1s(i) - wt0 + 1;
                        nwf = size(pw(w).phasediff_hilbert.phasedata, 2);
                        rel_t0 = max(1, rel_t0);
                        rel_t1 = min(nwf, rel_t1);

                        for p = 1:numel(phase_pairs)
                            pname = phase_pairs{p};
                            pdata = pw(w).phasediff_hilbert.(pname).data(rel_t0:rel_t1);
                            valid_pd = pdata(~isnan(pdata));
                            if ~isempty(valid_pd)
                                step_struct.(['phase_' pname])(step_ct) = circ_mean(valid_pd(:));
                            end
                        end

                        for p = 1:numel(abs_phase_pairs)
                            pname = abs_phase_pairs{p};
                            pdata = pw(w).phasediff_hilbert.(pname).data(rel_t0:rel_t1);
                            valid_pd = pdata(~isnan(pdata));
                            if ~isempty(valid_pd)
                                step_struct.(abs_phase_outnames{p})(step_ct) = mean(valid_pd);
                            end
                        end

                        for g = 1:numel(phase_groups)
                            gname = phase_groups{g};
                            pairs_in_group = phase_group_pairs.(gname){1};
                            all_data = [];
                            for pg = 1:numel(pairs_in_group)
                                pdata = pw(w).phasediff_hilbert.(pairs_in_group{pg}).data(rel_t0:rel_t1);
                                all_data = [all_data, pdata]; %#ok<AGROW>
                            end
                            valid_pd = all_data(~isnan(all_data));
                            if ~isempty(valid_pd)
                                step_struct.(['phase_' gname])(step_ct) = circ_mean(valid_pd(:));
                            end
                        end

                        all_data = [];
                        for pg = 1:numel(abs_group_pairs)
                            pdata = pw(w).phasediff_hilbert.(abs_group_pairs{pg}).data(rel_t0:rel_t1);
                            all_data = [all_data, pdata]; %#ok<AGROW>
                        end
                        valid_pd = all_data(~isnan(all_data));
                        if ~isempty(valid_pd)
                            step_struct.absphase_contra_L2R_3(step_ct) = mean(valid_pd);
                        end

                        phase_sidx = walk_struct.(['phase_sidx_' ln])(w);
                        phase_eidx = walk_struct.(['phase_eidx_' ln])(w);
                        if ~isnan(phase_sidx) && ~isnan(phase_eidx)
                            step_struct.phase_valid(step_ct) = (rel_t0 >= phase_sidx) && (rel_t1 <= phase_eidx);
                        end
                    end
                end

                if mod(w, 200) == 0
                    fprintf('  Processed %d/%d walks\n', w, nwalks);
                end
            end

            % Inherit walk_valid into step_struct via walk_idx
            valid_idx = step_struct.walk_idx(1:step_ct);
            step_struct.walk_valid = false(1, total_steps);
            step_struct.walk_valid(1:step_ct) = walk_struct.walk_valid(valid_idx);

            fprintf('  Built walk_struct: %d walks (%d valid), step_struct: %d steps\n', ...
                nwalks, sum(walk_struct.walk_valid), step_ct);
        end

        function statsperexp = combineBoutMetrics(obj, bout_metrics, led_label, statsperexp)
            % Flatten swing/stance scalar metrics, per-frame struct metrics,
            % and velmag-conditioned durations into statsperexp.

            flies = 'allflies';
            limbs = {'all_limbs', 'perlimb', 'pairs'};
            state = {'swing', 'stance'};
            ignorelist = {'fly', 'start_indices', 'end_indices', ...
                'Nboutspervelmagbin', 'velmagbincenters', ...
                'meanboutdurationsofvelmagbins', 'stdboutdurationsofvelmagbins'};

            for l = 1:numel(limbs)
                numlimb = numel(bout_metrics.(flies).(limbs{l}));

                for ll = 1:numlimb
                    limbname = obj.getLimbName(numlimb, ll);

                    for s = 1:numel(state)
                        featurefields = fields(bout_metrics.(flies).(limbs{l})(ll).(state{s}));
                        for f = 1:numel(featurefields)
                            if ismember(featurefields{f}, ignorelist)
                                continue;
                            end

                            fieldval = bout_metrics.(flies).(limbs{l})(ll).(state{s}).(featurefields{f});
                            funname = sprintf('%s__%s__%s__%s', featurefields{f}, state{s}, led_label, limbname);
                            currstruct = struct;

                            if ~isstruct(fieldval)
                                % Scalar per-bout data
                                currdata = [];
                                if ~isempty(fieldval)
                                    currdata = fieldval;
                                end
                                currstruct.mean = mean(currdata, 'omitnan');
                                currstruct.std = std(currdata, 'omitnan');
                                currstruct.Z = nnz(~isnan(currdata));
                            else
                                % Per-frame feature struct — use .mean field
                                currdata = [];
                                if ~isempty(fieldval)
                                    currdata = fieldval.mean;
                                end
                                currstruct.mean = mean(currdata, 'omitnan');
                                currstruct.std = std(currdata, 'omitnan');
                                currstruct.Z = nnz(~isnan(currdata));
                            end
                            statsperexp.(funname) = currstruct;
                        end

                        % Velmag-conditioned durations
                        bins = [3, 10];
                        for b = 1:numel(bins)
                            funname = sprintf('durations_time__%s__%s__%s__bin%svelmag', state{s}, led_label, limbname, num2str(bins(b)));
                            currstruct = struct;
                            if ~isfield(bout_metrics.(flies).(limbs{l})(ll).(state{s}), 'stdboutdurationsofvelmagbins')
                                currstruct.mean = NaN;
                                currstruct.std = NaN;
                                currstruct.Z = NaN;
                            else
                                currstruct.mean = bout_metrics.(flies).(limbs{l})(ll).(state{s}).meanboutdurationsofvelmagbins(bins(b));
                                currstruct.std = bout_metrics.(flies).(limbs{l})(ll).(state{s}).stdboutdurationsofvelmagbins(bins(b));
                                currstruct.Z = bout_metrics.(flies).(limbs{l})(ll).(state{s}).Nboutspervelmagbin(bins(b));
                            end
                            statsperexp.(funname) = currstruct;
                        end
                    end
                end
            end
        end

        function statsperexp = combineStepMetrics(obj, bout_metrics, led_label, statsperexp)
            % Flatten step per-frame features, step features, and x/y step
            % features into statsperexp.

            flies = 'allflies';
            limbs = {'all_limbs', 'perlimb', 'pairs'};

            perframe_features_step = {'velmag_ctr', 'absdv_ctr', 'absdu_ctr', 'absdtheta', ...
                'left_vel', 'right_vel', 'forward_vel', 'backward_vel', ...
                'right_dtheta', 'left_dtheta', 'CoM_stability'};

            stepfeatures = {'durations_frames', 'durations_time', 'instataeous_frequency_steps', ...
                'amplitude_px', 'amplitude_BL', 'step_direction', 'distance_px', 'distance_BL', ...
                'speed_pxpers', 'speed_BLpers', 'length_px', 'length_BL'};

            xystepfeatures = {'AEP', 'AEP_BL', 'PEP', 'PEP_BL'};
            xyname = {'x', 'y'};

            for l = 1:numel(limbs)
                numlimb = numel(bout_metrics.(flies).(limbs{l}));

                for ll = 1:numlimb
                    limbname = obj.getLimbName(numlimb, ll);

                    % Per-frame features during steps
                    for pff = 1:numel(perframe_features_step)
                        funname = sprintf('%s__step__%s__%s', perframe_features_step{pff}, led_label, limbname);
                        currstruct = struct;
                        currdata = [];
                        if ~isempty(bout_metrics.(flies).(limbs{l})(ll).step.(perframe_features_step{pff}))
                            currdata = bout_metrics.(flies).(limbs{l})(ll).step.(perframe_features_step{pff}).mean;
                        end
                        currstruct.mean = mean(currdata, 'omitnan');
                        currstruct.std = std(currdata, 'omitnan');
                        currstruct.Z = nnz(~isnan(currdata));
                        statsperexp.(funname) = currstruct;
                    end

                    % Step-level features
                    for sf = 1:numel(stepfeatures)
                        funname = sprintf('%s__step__%s__%s', stepfeatures{sf}, led_label, limbname);
                        currstruct = struct;
                        currdata = [];
                        if ~isempty(bout_metrics.(flies).(limbs{l})(ll).step.stepfeatures.(stepfeatures{sf}))
                            currdata = bout_metrics.(flies).(limbs{l})(ll).step.stepfeatures.(stepfeatures{sf});
                        end
                        currstruct.mean = mean(currdata, 'omitnan');
                        currstruct.std = std(currdata, 'omitnan');
                        currstruct.Z = nnz(~isnan(currdata));
                        statsperexp.(funname) = currstruct;
                    end

                    % X,Y step features (AEP, PEP)
                    for xysf = 1:numel(xystepfeatures)
                        for xy = 1:2
                            funname_str = [xystepfeatures{xysf}, xyname{xy}];
                            funname = sprintf('%s__step__%s__%s', funname_str, led_label, limbname);
                            currstruct = struct;
                            currdata = [];
                            if ~isempty(bout_metrics.(flies).(limbs{l})(ll).step.stepfeatures.(xystepfeatures{xysf}))
                                currdata = bout_metrics.(flies).(limbs{l})(ll).step.stepfeatures.(xystepfeatures{xysf})(xy,:);
                            end
                            currstruct.mean = mean(currdata, 'omitnan');
                            currstruct.std = std(currdata, 'omitnan');
                            currstruct.Z = nnz(~isnan(currdata));
                            statsperexp.(funname) = currstruct;
                        end
                    end
                end
            end
        end

        function statsperexp = combineWalkMetrics(~, walk_metrics, led_label, statsperexp)
            % Flatten walk per-frame features into statsperexp.

            perframe_features = {'velmag_ctr', 'absdv_ctr', 'absdu_ctr', 'absdtheta', ...
                'left_vel', 'right_vel', 'forward_vel', 'backward_vel', ...
                'right_dtheta', 'left_dtheta', 'CoM_stability'};

            for pff = 1:numel(perframe_features)
                funname = sprintf('%s__walk__%s__all', perframe_features{pff}, led_label);
                currstruct = struct;
                currdata = [];
                if ~isempty(fields(walk_metrics.perexp))
                    currdata = walk_metrics.perexp.(perframe_features{pff}).frm_data_exp;
                end
                currstruct.mean = mean(currdata, 'omitnan');
                currstruct.std = std(currdata, 'omitnan');
                currstruct.Z = nnz(~isnan(currdata));
                statsperexp.(funname) = currstruct;
            end
        end

        function statsperexp = combinePhaseMetrics(~, walk_metrics, led_label, statsperexp)
            % Flatten phase difference metrics into statsperexp.
            % Uses circ_mean/circ_std for signed phase, regular mean/std for absolute phase.

            feature = 'phasediff_hilbert';
            stat = 'frm_data_exp';

            % Signed phase — use circular statistics
            subfeature_signed = {'tripods_4', 'ipsi_P2A_4', 'ipsi_ant_2', 'ipsi_post_2'};

            for sf = 1:numel(subfeature_signed)
                funname = sprintf('%s__walk__%s__%s', feature, led_label, subfeature_signed{sf});
                currstruct = struct;
                currdata = [];
                if ~isempty(fields(walk_metrics.perexp))
                    currdata = walk_metrics.perexp.(feature).(subfeature_signed{sf}).(stat);
                end
                valid = currdata(~isnan(currdata));
                currstruct.mean = circ_mean(valid);
                currstruct.std = circ_std(valid);
                currstruct.Z = nnz(~isnan(currdata));
                statsperexp.(funname) = currstruct;
            end

            % Absolute phase — use regular statistics
            subfeature_abs = {'abscontra_L2R_3', 'absRF_LF', 'absRM_LM', 'absRH_LH'};

            for sf = 1:numel(subfeature_abs)
                funname = sprintf('abs%s__walk__%s__%s', feature, led_label, subfeature_abs{sf});
                currstruct = struct;
                currdata = [];
                if ~isempty(fields(walk_metrics.perexp))
                    currdata = walk_metrics.perexp.(feature).(subfeature_abs{sf}).(stat);
                end
                currstruct.mean = mean(currdata, 'omitnan');
                currstruct.std = std(currdata, 'omitnan');
                currstruct.Z = nnz(~isnan(currdata));
                statsperexp.(funname) = currstruct;
            end
        end

    end

end

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

        % Settings
        debug = false
        expdir = ''
        nflies = 0

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
            parse(p, varargin{:});

            obj.debug = p.Results.debug;
            obj.expdir = p.Results.expdir;
            obj.binedges = p.Results.binedges;           

            % Initialize containers
            obj.restrictedBoutData = containers.Map();
            obj.boutMetrics = containers.Map();
            obj.validFrames = containers.Map();
            obj.walkMetrics = containers.Map();
            obj.locoStatsPerExp = struct();

            
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
        end

        function createValidFrames_LED(obj)
            obj.validFrames('led_on_traj') = obj.digitalindicator;

            obj.validFrames('led_off_traj') = obj.invertTrajectoryIndicator(obj.digitalindicator);
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
            walk_metrics = computeWalkMetrics(obj,digital_signal);
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

        function locostatsperexp = computeStatsPerExp(obj, conditions)
            % Flatten bout_metrics and walk_metrics into a single per-experiment
            % summary struct. LED condition is encoded in each field name
            % (e.g. feature__swing__LEDon__all, feature__walk__LEDoff__all).
            %
            % Usage:
            %   stats = analyzer.computeStatsPerExp();          % both ON and OFF
            %   stats = analyzer.computeStatsPerExp({'ON'});    % ON only
            %
            % Based on LocomotionCombinePerFrameStats.m (Alice)

            if nargin < 2
                conditions = {'ON', 'OFF'};
            end

            % Map condition labels to stored keys and LED field name labels
            condMap = struct();
            condMap.ON.boutKey  = 'walking_stimON_traj';
            condMap.ON.walkKey  = 'led_on_traj';
            condMap.ON.label    = 'LEDon';
            condMap.OFF.boutKey = 'walking_stimOFF_traj';
            condMap.OFF.walkKey = 'led_off_traj';
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

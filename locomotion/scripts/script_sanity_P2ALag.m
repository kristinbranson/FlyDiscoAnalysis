%% script_sanity_P2ALag.m
% Sanity-check plots for the posterior-to-anterior (P2A) onset lag metrics
% (Pliftoff2Aliftoff_lag, Ptouchdown2Aliftoff_lag) on the test experiment.
%
% Runs the locomotion walk-metrics pipeline (no onfloor filtering, for
% simplicity), builds walk_struct for ON/OFF, and plots the distribution of
% per-walk lag values for all 7 sub-fields (4 pairs + H_to_M, M_to_F, all),
% overlaying LED ON vs OFF. Shared (fixed) x- and y-limits per metric so the
% panels are directly comparable.
%
% Expected sanity checks:
%   - lags non-negative (forward-only pairing window), in [0, ~step period]
%   - touchdown->liftoff lag (Rule 2) shorter than liftoff->liftoff lag
%   - distributions roughly unimodal; H_to_M / M_to_F resemble their pairs

addpath('/groups/branson/home/robiea/Code_versioned/FlyDiscoAnalysis');
modpath;

%% Config
settingsdir = '/groups/branson/home/robiea/Code_versioned/BransonFlyDiscoSettings/settings';
analysis_protocol = '20260326_flybubble_LED_VNC2';
expdir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20260112_testing_locomotioncomputeperframestats/VNC2_JRC_SS57983_RigD_20230913T120134';

plotdir = sprintf('/groups/branson/home/robiea/Projects_data/FlyDisco/Locomotion_analysis/claudeplots_%s_P2ALag', datestr(now,'yyyymmdd'));
saveplots = true;
if saveplots && ~exist(plotdir,'dir'), mkdir(plotdir); end

metrics = {'Pliftoff2Aliftoff_lag', 'Ptouchdown2Aliftoff_lag'};
subfields = {'RH_to_RM','RM_to_RF','LH_to_LM','LM_to_LF','H_to_M','M_to_F','all'};

%% Run pipeline (no onfloor filtering)
fprintf('Initializing trx...\n');
trx = FBATrx('analysis_protocol',analysis_protocol,'settingsdir',settingsdir,...
    'datalocparamsfilestr','dataloc_params.txt');
trx.AddExpDir(expdir,'dooverwrite',false,'openmovie',false);

aptfile = trx.dataloc_params.apttrkfilestr;
aptdata = TrkFile.load(fullfile(expdir,aptfile));

stageparamsfile = fullfile(trx.settingsdir,trx.analysis_protocol,trx.dataloc_params.locomotionmetricsparamsfilestr);
stage_params = ReadParams(stageparamsfile);
legtip_landmarknums = stage_params.legtip_landmarknums;

load(fullfile(expdir,'tips_velmag.mat'),'tips_velmag');
load(fullfile(expdir,'tips_pos_body.mat'),'tips_pos_body');

groundcontact = compute_groundcontact(tips_velmag, ...
    'pairs', stage_params.pairs, ...
    'gc_threshold_low', stage_params.gc_threshold_low, ...
    'gc_threshold_high', stage_params.gc_threshold_high, ...
    'minimum_bout', stage_params.minimum_bout_groundcontact);

[~,walking_scores] = LoadScoresFromFile(trx,'scores_Walk2',1);
indicatordata = trx.getIndicatorLED(1);
digitalindicator = indicatordata.indicatordigital;

fprintf('Building LimbBoutAnalyzer and running walk metrics...\n');
loco_analyzer = LimbBoutAnalyzer(trx, aptdata, tips_pos_body, legtip_landmarknums, ...
    groundcontact, digitalindicator, walking_scores, ...
    'phase_methods', {'phasediff_hilbert'}, 'expdir', expdir);
loco_analyzer.analyzeWalkAndStimConditions();
loco_analyzer.analyzeBoutAndStimConditions();
loco_analyzer.buildWalkStruct({'ON','OFF'});

ws_on  = loco_analyzer.walkStruct.ON.walk_struct;
ws_off = loco_analyzer.walkStruct.OFF.walk_struct;

%% Plot one figure per metric
for mi = 1:numel(metrics)
    mn = metrics{mi};

    % pass 1: shared x-edges (0 to 99th pct over all sub-fields/conditions)
    allvals = [];
    for si = 1:numel(subfields)
        fn = [mn '_' subfields{si}];
        allvals = [allvals, ws_off.(fn), ws_on.(fn)]; %#ok<AGROW>
    end
    allvals = allvals(~isnan(allvals));
    if isempty(allvals)
        fprintf('  %s: no valid data, skipping plot\n', mn);
        continue
    end
    xmax = prctile(allvals, 99);
    edges = linspace(0, xmax, 40);

    % pass 1b: shared y-limit (max probability over all panels/conditions)
    ymax = 0;
    for si = 1:numel(subfields)
        fn = [mn '_' subfields{si}];
        for d = {ws_off.(fn), ws_on.(fn)}
            v = d{1}; v = v(~isnan(v));
            if ~isempty(v)
                h = histcounts(v, edges, 'Normalization','probability');
                ymax = max(ymax, max(h));
            end
        end
    end
    if ymax == 0, ymax = 1; end

    % pass 2: plot
    fig = figure('Position',[100 100 1400 600],'Color','w');
    tl = tiledlayout(2,4,'TileSpacing','compact','Padding','compact');
    for si = 1:numel(subfields)
        fn = [mn '_' subfields{si}];
        voff = ws_off.(fn); voff = voff(~isnan(voff));
        von  = ws_on.(fn);  von  = von(~isnan(von));

        nexttile; hold on
        histogram(voff, edges, 'Normalization','probability', ...
            'DisplayStyle','stairs','EdgeColor',[0 0 0],'LineWidth',1.5);
        histogram(von, edges, 'Normalization','probability', ...
            'DisplayStyle','stairs','EdgeColor',[0.85 0.2 0.2],'LineWidth',1.5);
        if ~isempty(voff)
            xline(mean(voff),'-','Color',[0 0 0],'LineWidth',1);
        end
        if ~isempty(von)
            xline(mean(von),'-','Color',[0.85 0.2 0.2],'LineWidth',1);
        end
        xlim([0 xmax]); ylim([0 ymax*1.05]);
        title(sprintf('%s (n_{off}=%d, n_{on}=%d)', subfields{si}, numel(voff), numel(von)), ...
            'Interpreter','none');
        if si == 1, legend({'LED off','LED on'},'Location','northeast','Box','off'); end
        set(gca,'TickLabelInterpreter','none');
    end
    xlabel(tl, 'per-walk lag (ms)');
    ylabel(tl, 'probability');
    title(tl, sprintf('%s — %s', mn, 'VNC2\_JRC\_SS57983\_RigD\_20230913T120134'), 'Interpreter','none');

    if saveplots
        outfile = fullfile(plotdir, sprintf('sanity_%s.png', mn));
        exportgraphics(fig, outfile, 'Resolution',150);
        fprintf('  saved %s\n', outfile);
    end
end

fprintf('Done. Plots in %s\n', plotdir);

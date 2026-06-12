%% script_explore_P2Alag_window.m
% Empirically choose the P2A pairing window. For every floor walk, find the
% NEAREST anterior swing onset to each reference onset within a WIDE +/-1
% reference period (uncapped), signed, expressed as a FRACTION of the
% reference period (period-independent). Histogram per metric x pair so we
% can see where lags actually fall and pick a backward/forward cap ratio
% (candidate: -0.25 .. +0.75) without the cap biasing the distribution.
%
% Positive fraction = anterior follows posterior (normal wave).

addpath('/groups/branson/home/robiea/Code_versioned/FlyDiscoAnalysis');
modpath;

%% Config (freshly-processed control copy, floor walks)
settingsdir = '/groups/branson/home/robiea/Code_versioned/BransonFlyDiscoSettings/settings';
analysis_protocol = '20260326_flybubble_LED_VNC2';
expdir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20260612_testingLag/VNC2_YNA_K_162984_RigA_20220419T110856';
plotdir = sprintf('/groups/branson/home/robiea/Projects_data/FlyDisco/Locomotion_analysis/claudeplots_%s_P2ALag_verify', datestr(now,'yyyymmdd'));
if ~exist(plotdir,'dir'), mkdir(plotdir); end

pairs = {[3 2],[2 1]};  pairnames = {'RH_to_RM','RM_to_RF'};
refevents = {'swing','stance'};
metricnames = {'Pliftoff2Aliftoff','Ptouchdown2Aliftoff'};
cap_pre = 0.25; cap_post = 0.75;   % candidate caps (fraction of period)

%% Build analyzer (onfloor) and get floor walks
trx = FBATrx('analysis_protocol',analysis_protocol,'settingsdir',settingsdir,'datalocparamsfilestr','dataloc_params.txt');
trx.AddExpDir(expdir,'dooverwrite',false,'openmovie',false);
aptdata = TrkFile.load(fullfile(expdir,trx.dataloc_params.apttrkfilestr));
stage_params = ReadParams(fullfile(trx.settingsdir,trx.analysis_protocol,trx.dataloc_params.locomotionmetricsparamsfilestr));
load(fullfile(expdir,'tips_velmag.mat'),'tips_velmag');
load(fullfile(expdir,'tips_pos_body.mat'),'tips_pos_body');
groundcontact = compute_groundcontact(tips_velmag,'pairs',stage_params.pairs, ...
    'gc_threshold_low',stage_params.gc_threshold_low,'gc_threshold_high',stage_params.gc_threshold_high, ...
    'minimum_bout',stage_params.minimum_bout_groundcontact);
[~,walking_scores]    = LoadScoresFromFile(trx,'scores_Walk2',1);
[~,onceiling_scores]  = LoadScoresFromFile(trx,'scores_onceiling_resnet_v2',1);
[~,nottracking_scores]= LoadScoresFromFile(trx,'scores_nottracking',1);
digitalindicator = trx.getIndicatorLED(1).indicatordigital;
loco = LimbBoutAnalyzer(trx, aptdata, tips_pos_body, stage_params.legtip_landmarknums, ...
    groundcontact, digitalindicator, walking_scores, 'phase_methods',{'phasediff_hilbert'}, ...
    'expdir',expdir,'do_onfloor_filtering',true,'onceiling_scores',onceiling_scores, ...
    'nottracking_scores',nottracking_scores,'frac_onfloor_threshold',stage_params.frac_onfloor_threshold);
loco.analyzeBoutAndStimConditions_onfloor();
loco.analyzeWalkAndStimConditions_onfloor();
pw = loco.walkMetrics('led_off_onfloor_traj').perwalk;
all_ts = trx.movie_timestamps{1};
fprintf('%d floor walks\n', numel(pw));

%% Accumulate uncapped nearest signed fraction-of-period per metric x pair
frac = struct();  % frac.(metric).(pair) = vector
for m = 1:2, for p = 1:2, frac.(metricnames{m}).(pairnames{p}) = []; end, end

for w = 1:numel(pw)
    fly = pw(w).fly; t0 = pw(w).walk_t0; t1 = pw(w).walk_t1;
    bd = loco.limbBoutData(fly);
    for m = 1:2
        for p = 1:2
            refl = pairs{p}(1); antl = pairs{p}(2);
            rb = bd.perlimb(refl).(refevents{m});
            idx = rb.start_indices >= t0 & rb.end_indices <= t1;
            ron = sort(rb.start_indices(idx));
            if numel(ron) < 2, continue; end
            P = median(diff(ron));                       % frames
            aon = sort(bd.perlimb(antl).swing.start_indices(:)');
            f = nan(1,numel(ron));
            for s = 1:numel(ron)
                r = ron(s);
                cand = aon(aon >= r - P & aon <= r + P);  % wide +/- 1 period
                if isempty(cand), continue; end
                [~,k] = min(abs(cand - r));
                f(s) = (cand(k) - r) / P;                 % signed fraction of period
            end
            frac.(metricnames{m}).(pairnames{p}) = [frac.(metricnames{m}).(pairnames{p}), f(~isnan(f))]; %#ok<AGROW>
        end
    end
end

%% Histograms: 2 metrics x 2 pairs, fraction of period
edges = -1:0.04:1;
fig = figure('Position',[80 80 1300 800],'Color','w','Name','P2A lag window exploration','NumberTitle','off');
tl = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
for m = 1:2
    for p = 1:2
        nexttile; hold on;
        d = frac.(metricnames{m}).(pairnames{p});
        histogram(d, edges, 'Normalization','probability','FaceColor',[0.3 0.5 0.8],'EdgeColor','none');
        xline(0,'k-'); xline(-0.5,'k:'); xline(0.5,'k:');
        xline(-cap_pre,'r--','LineWidth',1.5); xline(cap_post,'r--','LineWidth',1.5);
        within = mean(d >= -cap_pre & d < cap_post) * 100;
        title(sprintf('%s  %s   (n=%d, median=%.2f, %.0f%% in [-%.2f,%.2f))', ...
            metricnames{m}, pairnames{p}, numel(d), median(d), within, cap_pre, cap_post), ...
            'Interpreter','none','FontSize',10);
        xlabel('signed lag / reference period'); ylabel('probability'); xlim([-1 1]);
        set(gca,'TickLabelInterpreter','none');
    end
end
title(tl, sprintf('P2A nearest-onset signed lag (fraction of period), all floor walks\nred dashed = candidate cap [-%.2f, %.2f];  black dotted = +/-0.5', cap_pre, cap_post), ...
    'FontWeight','bold','Interpreter','none');
outfile = fullfile(plotdir,'P2Alag_window_histograms.png');
exportgraphics(fig, outfile, 'Resolution',150);
fprintf('saved %s\n', outfile);

%% Console summary: how much falls outside candidate cap, and at edges
fprintf('\nfraction-of-period summary (candidate cap [-%.2f, %.2f)):\n', cap_pre, cap_post);
for m = 1:2
    for p = 1:2
        d = frac.(metricnames{m}).(pairnames{p});
        fprintf('  %-20s %-9s n=%5d  median=%+.2f  pre-edge(<-%.2f)=%.1f%%  post-edge(>=%.2f)=%.1f%%  beyond|0.5|=%.1f%%\n', ...
            metricnames{m}, pairnames{p}, numel(d), median(d), cap_pre, mean(d < -cap_pre)*100, ...
            cap_post, mean(d >= cap_post)*100, mean(abs(d) > 0.5)*100);
    end
end

%% script_plot_P2ALag_distributions.m
% Distributions of the FINAL P2A lag metrics (as produced by the pipeline with
% the chosen pairings: liftoff=forward, touchdown=signed +/-0.5P), pooled over
% per-event values across all floor walks of the test control. Units = ms.
%
%   row 1: Pliftoff2Aliftoff_lag  (forward, >= 0)   -> metachronal wave
%   row 2: Ptouchdown2Aliftoff_lag (signed, ~0)     -> Cruse Rule 2 latency

addpath('/groups/branson/home/robiea/Code_versioned/FlyDiscoAnalysis');
modpath;

settingsdir = '/groups/branson/home/robiea/Code_versioned/BransonFlyDiscoSettings/settings';
analysis_protocol = '20260326_flybubble_LED_VNC2';
expdir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20260612_testingLag/VNC2_YNA_K_162984_RigA_20220419T110856';
plotdir = sprintf('/groups/branson/home/robiea/Projects_data/FlyDisco/Locomotion_analysis/claudeplots_%s_P2ALag_verify', datestr(now,'yyyymmdd'));
if ~exist(plotdir,'dir'), mkdir(plotdir); end

metrics = {'Pliftoff2Aliftoff_lag','Ptouchdown2Aliftoff_lag'};
metrictitles = {'P liftoff \rightarrow A liftoff  (forward)','P touchdown \rightarrow A liftoff  (signed)'};
subfields = {'RH_to_RM','RM_to_RF','LH_to_LM','LM_to_LF','all'};

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
fprintf('%d floor walks\n', numel(pw));

%% Pool per-event .data (ms) across all floor walks
D = struct();
for mi = 1:2
    for si = 1:numel(subfields)
        c = cell(1,numel(pw));
        for w = 1:numel(pw), c{w} = pw(w).(metrics{mi}).(subfields{si}).data; end
        d = horzcat(c{:}); D.(metrics{mi}).(subfields{si}) = d(~isnan(d));
    end
end

%% Fixed x-range per metric (row): forward = [0, p99]; signed = symmetric +/-p99(|.|)
xl = cell(1,2);
dall1 = D.(metrics{1}).all;  xl{1} = [0, prctile(dall1,99)];
dall2 = D.(metrics{2}).all;  xx = prctile(abs(dall2),99);  xl{2} = [-xx, xx];

%% Plot: 2 metrics x 5 subfields
fig = figure('Position',[60 60 1600 700],'Color','w','Name','P2A lag distributions (new pairing)','NumberTitle','off');
tl = tiledlayout(2,numel(subfields),'TileSpacing','compact','Padding','compact');
for mi = 1:2
    edges = linspace(xl{mi}(1), xl{mi}(2), 50);
    for si = 1:numel(subfields)
        nexttile; hold on;
        d = D.(metrics{mi}).(subfields{si});
        histogram(d, edges, 'Normalization','probability','FaceColor',[0.3 0.5 0.8],'EdgeColor','none');
        if mi == 2, xline(0,'k-'); end
        xline(mean(d),'r-','LineWidth',1.5);
        title(sprintf('%s  (n=%d)\nmean=%.1f  median=%.1f ms', subfields{si}, numel(d), mean(d), median(d)), ...
            'Interpreter','none','FontSize',9);
        xlim(xl{mi});
        if si == 1, ylabel(metrictitles{mi},'Interpreter','tex','FontWeight','bold'); end
        if mi == 2, xlabel('signed lag (ms)'); else, xlabel('lag (ms)'); end
        set(gca,'TickLabelInterpreter','none');
    end
end
title(tl, sprintf('P2A lag distributions (new pairing), per-event pooled over %d floor walks\nred = mean', numel(pw)), ...
    'FontWeight','bold','Interpreter','none');
outfile = fullfile(plotdir,'P2Alag_distributions_newpairing.png');
exportgraphics(fig, outfile, 'Resolution',150);
fprintf('saved %s\n', outfile);

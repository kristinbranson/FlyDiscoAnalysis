%% script_verify_P2ALag_examples.m
% Per-example verification of the P2A onset-lag metrics on the test experiment.
%
% For a handful of example walks (clean / missing-match / terminal edge bout)
% it draws the right-side ipsilateral triple (RH->RM->RF) swing/stance
% timeline and, for BOTH metrics, connects each reference onset to its matched
% anterior liftoff, annotated with the lag computed INDEPENDENTLY here (not by
% computeP2ALag) versus the pipeline's stored per-event value. A console table
% reports indep-vs-code agreement per event.
%
%   top panel    : Pliftoff2Aliftoff_lag  (posterior liftoff  -> anterior liftoff)
%   bottom panel : Ptouchdown2Aliftoff_lag (posterior touchdown -> anterior liftoff)
%
% Right triple limb indices (legorder {RF,RM,RH,LH,LM,LF}): RH=3, RM=2, RF=1.
% Pairs: RH_to_RM (3->2), RM_to_RF (2->1).

addpath('/groups/branson/home/robiea/Code_versioned/FlyDiscoAnalysis');
modpath;

%% Config
settingsdir = '/groups/branson/home/robiea/Code_versioned/BransonFlyDiscoSettings/settings';
analysis_protocol = '20260326_flybubble_LED_VNC2';
% freshly-processed control (symbolic copy + current locomotion stage, 2026-06-12);
% has onceiling/nottracking scores so we can filter to floor-only walks
expdir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20260612_testingLag/VNC2_YNA_K_162984_RigA_20220419T110856';
plotdir = sprintf('/groups/branson/home/robiea/Projects_data/FlyDisco/Locomotion_analysis/claudeplots_%s_P2ALag_verify', datestr(now,'yyyymmdd'));
saveplots = true;
if saveplots && ~exist(plotdir,'dir'), mkdir(plotdir); end

% right triple, plotted posterior(bottom)->anterior(top)
triple   = [3 2 1];            % RH, RM, RF -> y = 1,2,3
tripnames = {'RH','RM','RF'};
pairs    = {[3 2],[2 1]};      % RH->RM, RM->RF (ref -> ant)
pairnames = {'RH_to_RM','RM_to_RF'};
metrics  = {'Pliftoff2Aliftoff_lag','Ptouchdown2Aliftoff_lag'};
refevents = {'swing','stance'};       % reference event per metric
pairings  = {'forward','signed'};     % pairing rule per metric (must match computeP2ALag)
metrictitles = {'P liftoff \rightarrow A liftoff  (forward)','P touchdown \rightarrow A liftoff  (signed)'};

%% Run pipeline (no onfloor filtering)
fprintf('Initializing trx and running walk metrics...\n');
trx = FBATrx('analysis_protocol',analysis_protocol,'settingsdir',settingsdir,...
    'datalocparamsfilestr','dataloc_params.txt');
trx.AddExpDir(expdir,'dooverwrite',false,'openmovie',false);
aptdata = TrkFile.load(fullfile(expdir,trx.dataloc_params.apttrkfilestr));
trxs = load(fullfile(expdir,'registered_trx.mat')); trx_s = trxs.trx;   % pixel x/y/theta + firstframe
[readframe, ~, movie_fid] = get_readframe_fcn(fullfile(expdir, trx.dataloc_params.moviefilestr));
stage_params = ReadParams(fullfile(trx.settingsdir,trx.analysis_protocol,trx.dataloc_params.locomotionmetricsparamsfilestr));
load(fullfile(expdir,'tips_velmag.mat'),'tips_velmag');
load(fullfile(expdir,'tips_pos_body.mat'),'tips_pos_body');
groundcontact = compute_groundcontact(tips_velmag,'pairs',stage_params.pairs, ...
    'gc_threshold_low',stage_params.gc_threshold_low,'gc_threshold_high',stage_params.gc_threshold_high, ...
    'minimum_bout',stage_params.minimum_bout_groundcontact);
[~,walking_scores] = LoadScoresFromFile(trx,'scores_Walk2',1);
[~,onceiling_scores]  = LoadScoresFromFile(trx,'scores_onceiling_resnet_v2',1);
[~,nottracking_scores] = LoadScoresFromFile(trx,'scores_nottracking',1);
digitalindicator = trx.getIndicatorLED(1).indicatordigital;
loco = LimbBoutAnalyzer(trx, aptdata, tips_pos_body, stage_params.legtip_landmarknums, ...
    groundcontact, digitalindicator, walking_scores, 'phase_methods',{'phasediff_hilbert'}, ...
    'expdir',expdir, 'do_onfloor_filtering',true, ...
    'onceiling_scores',onceiling_scores, 'nottracking_scores',nottracking_scores, ...
    'frac_onfloor_threshold',stage_params.frac_onfloor_threshold);
loco.analyzeBoutAndStimConditions_onfloor();
loco.analyzeWalkAndStimConditions_onfloor();
wm = loco.walkMetrics('led_off_onfloor_traj');
pw = wm.perwalk;
% per-walk fraction of frames on the floor (for selecting fully-floor examples)
onfloor_frac = @(w) mean(loco.is_onfloor_perframe{pw(w).fly}( ...
    pw(w).walk_t0 : min(pw(w).walk_t1, numel(loco.is_onfloor_perframe{pw(w).fly})) ));

all_ts = trx.movie_timestamps{1};

%% Select example walks: clean / missing-match / terminal edge
nan_in = @(w) any(isnan(pw(w).Pliftoff2Aliftoff_lag.RH_to_RM.data)) || ...
              any(isnan(pw(w).Pliftoff2Aliftoff_lag.RM_to_RF.data));

ex = struct('w',{},'label',{});
ONFLOOR = 0.99;   % require examples to be essentially fully on the floor
% clean = richest no-NaN floor walk, capped to <=12 RH_to_RM events for legibility
best_clean = 0; best_n = -1;
for w = 1:numel(pw)
    if ~nan_in(w) && onfloor_frac(w) >= ONFLOOR
        n1 = pw(w).Pliftoff2Aliftoff_lag.RH_to_RM.n;
        nn = n1 + pw(w).Pliftoff2Aliftoff_lag.RM_to_RF.n;
        if n1 <= 12 && nn > best_n, best_n = nn; best_clean = w; end
    end
end
if best_clean > 0, ex(end+1) = struct('w',best_clean,'label','clean floor walk'); end %#ok<SAGROW>
for w = 1:numel(pw)
    if nan_in(w) && onfloor_frac(w) >= ONFLOOR && ...
            (pw(w).Pliftoff2Aliftoff_lag.RH_to_RM.n + pw(w).Pliftoff2Aliftoff_lag.RM_to_RF.n) >= 2
        ex(end+1) = struct('w',w,'label','floor walk with missing match'); break; %#ok<SAGROW>
    end
end
% longest fully-floor walk (most events) as a third example
best_long = 0; best_ln = -1;
for w = 1:numel(pw)
    if onfloor_frac(w) >= ONFLOOR
        nn = pw(w).Pliftoff2Aliftoff_lag.RH_to_RM.n + pw(w).Pliftoff2Aliftoff_lag.RM_to_RF.n;
        if nn > best_ln, best_ln = nn; best_long = w; end
    end
end
if best_long > 0, ex(end+1) = struct('w',best_long,'label','long floor walk'); end %#ok<SAGROW>

isclose = @(a,b) (isnan(a)&&isnan(b)) || abs(a-b) < 1e-9;

% close P2A figures from a prior run of this script (leave other figures alone)
close(findobj('Type','figure','-regexp','Name','^P2A example'));

%% Build a figure per example
nmismatch_total = 0;
for e = 1:numel(ex)
    w = ex(e).w;
    fly = pw(w).fly;
    walk_t0 = pw(w).walk_t0; walk_t1 = pw(w).walk_t1;
    cf_ts = all_ts(trx.firstframes(fly):trx.endframes(fly));
    bd = loco.limbBoutData(fly);
    tms = @(f) (cf_ts(min(f,numel(cf_ts))) - cf_ts(walk_t0)) * 1000;  % frame -> ms from walk start
    % context window: show a margin before/after the classified walk bout
    margin = max(20, round(0.2*(walk_t1 - walk_t0)));
    disp_t0 = max(1, walk_t0 - margin);
    disp_t1 = min(walk_t1 + margin, numel(cf_ts));

    fprintf('\n===== Example %d: %s | fly %d, walk %d (frames %d-%d) =====\n', ...
        e, ex(e).label, fly, w, walk_t0, walk_t1);

    fig = figure('Position',[80 80 1500 920],'Color','w', ...
        'Name',sprintf('P2A example %d: %s  (fly %d, walk %d)', e, ex(e).label, fly, w),'NumberTitle','off');
    tl = tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

    % --- per-frame speed & turning over the walk (time-aligned) ---
    axpf = nexttile; hold(axpf,'on');
    velmag = trx.GetPerFrameData('velmag_ctr', fly);
    absdth = trx.GetPerFrameData('absdtheta', fly);
    fr = disp_t0:min(disp_t1, min(numel(velmag),numel(absdth)));
    xx = (cf_ts(fr) - cf_ts(walk_t0))*1000;
    yyaxis(axpf,'left');  plot(axpf, xx, velmag(fr), '-','LineWidth',1.2); ylabel(axpf,'velmag (mm/s)');
    yyaxis(axpf,'right'); plot(axpf, xx, absdth(fr), '-','LineWidth',1.2); ylabel(axpf,'absdtheta (rad/s)');
    xline(axpf, 0, '--', 'walk start','Color',[0.15 0.15 0.15],'LineWidth',1.2,'LabelVerticalAlignment','bottom','Interpreter','none');
    xline(axpf, tms(walk_t1), '--', 'walk end','Color',[0.15 0.15 0.15],'LineWidth',1.2,'LabelVerticalAlignment','bottom','Interpreter','none');
    xlim(axpf,[tms(disp_t0)-2 tms(disp_t1)+2]);
    title(axpf, sprintf('speed & turning  (walk mean: velmag=%.1f mm/s, absdtheta=%.2f rad/s)', ...
        pw(w).velmag_ctr.mean, pw(w).absdtheta.mean), 'Interpreter','none');
    grid(axpf,'on'); set(axpf,'TickLabelInterpreter','none');

    for m = 1:numel(metrics)
        ax = nexttile; hold(ax,'on');

        % --- swing + stance occupancy for the triple over the context window
        %     [disp_t0, disp_t1]; markers drawn for any onset in-window ---
        for i = 1:3
            limb = triple(i); y = i;
            % stance (on ground) as grey background
            st = bd.perlimb(limb).stance;
            for k = 1:numel(st.start_indices)
                s = st.start_indices(k); en = st.end_indices(k);
                if en <= disp_t0 || s >= disp_t1, continue; end
                plot(ax,[tms(max(s,disp_t0)) tms(min(en,disp_t1))],[y y],'-','LineWidth',7,'Color',[0.82 0.82 0.82]);
            end
            % swing (in air) as blue, with liftoff/touchdown markers
            sw = bd.perlimb(limb).swing;
            for k = 1:numel(sw.start_indices)
                s = sw.start_indices(k); en = sw.end_indices(k);
                if en <= disp_t0 || s >= disp_t1, continue; end
                plot(ax,[tms(max(s,disp_t0)) tms(min(en,disp_t1))],[y y],'-','LineWidth',7,'Color',[0.55 0.75 0.95]);
                if s >= disp_t0 && s <= disp_t1
                    plot(ax,tms(s),y,'^','MarkerSize',8,'MarkerFaceColor',[0.1 0.4 0.8],'MarkerEdgeColor','k'); % liftoff
                end
                if en >= disp_t0 && en <= disp_t1
                    plot(ax,tms(en),y,'v','MarkerSize',8,'MarkerFaceColor',[0.85 0.5 0.1],'MarkerEdgeColor','k'); % touchdown
                end
            end
        end
        % classified walk bout boundaries
        xline(ax, 0, '--', 'Color',[0.15 0.15 0.15], 'LineWidth',1.2);
        xline(ax, tms(walk_t1), '--', 'Color',[0.15 0.15 0.15], 'LineWidth',1.2);

        % --- connectors per pair: indep recompute vs pipeline .data ---
        fprintf('  [%s]\n', metrics{m});
        for p = 1:numel(pairs)
            refl = pairs{p}(1); antl = pairs{p}(2);
            yref = find(triple==refl); yant = find(triple==antl);
            [lag, reff, matchf] = indep_p2a_lag(bd, refl, antl, walk_t0, walk_t1, cf_ts, refevents{m}, pairings{m});
            code = pw(w).(metrics{m}).(pairnames{p}).data;
            n = min(numel(lag), numel(code));
            fprintf('    %s: %d events\n', pairnames{p}, n);
            for s = 1:n
                ok = isclose(lag(s), code(s));
                if ~ok, nmismatch_total = nmismatch_total + 1; end
                fprintf('      ev%d ref@%d match@%s  indep=%s  code=%s  %s\n', s, reff(s), ...
                    fnum(matchf(s)), fnum(lag(s)), fnum(code(s)), tern(ok,'OK','MISMATCH'));
                % draw connector (only if matched)
                if ~isnan(matchf(s))
                    xr = tms(reff(s)); xa = tms(matchf(s));
                    col = tern(ok,[0 0.6 0],[0.85 0 0]);
                    plot(ax,[xr xa],[yref yant],'-','LineWidth',1.5,'Color',col);
                    plot(ax,xr,yref,'o','MarkerSize',6,'MarkerFaceColor',col,'MarkerEdgeColor','k');
                    text(ax,(xr+xa)/2,(yref+yant)/2,sprintf(' %.1f/%.1f',lag(s),code(s)), ...
                        'FontSize',8,'Color',col,'FontWeight','bold','Interpreter','none');
                else
                    xr = tms(reff(s));
                    plot(ax,xr,yref,'x','MarkerSize',10,'Color',[0.85 0 0],'LineWidth',2);
                    text(ax,xr,yref+0.18,' NaN','FontSize',8,'Color',[0.85 0 0],'Interpreter','none');
                end
            end
        end

        set(ax,'YTick',1:3,'YTickLabel',tripnames,'TickLabelInterpreter','none');
        ylim(ax,[0.5 3.5]); xlim(ax,[tms(disp_t0)-2 tms(disp_t1)+2]);
        ylabel(ax,'leg (posterior \rightarrow anterior)');
        if m == numel(metrics), xlabel(ax,'time from walk start (ms)'); end
        title(ax, metrictitles{m}, 'Interpreter','tex','FontSize',12,'FontWeight','bold');
        grid(ax,'on');
    end

    title(tl, sprintf('EXAMPLE %d:  %s\nfly %d,  walk %d  (frames %d-%d)', ...
        e, ex(e).label, fly, w, walk_t0, walk_t1), 'FontSize',14,'FontWeight','bold','Interpreter','none');
    annotation(fig,'textbox',[0.02 0.0 0.96 0.028], ...
        'String',['blue = swing,  grey = stance      |      up-triangle = liftoff,  down-triangle = touchdown' ...
        '      |      connector label = independent / code lag (ms)'], ...
        'EdgeColor','none','HorizontalAlignment','center','FontSize',9,'Interpreter','none');

    if saveplots
        outfile = fullfile(plotdir, sprintf('verify_ex%d_%s.png', e, matlab.lang.makeValidName(ex(e).label)));
        exportgraphics(fig, outfile, 'Resolution',150);
        fprintf('  saved %s\n', outfile);
    end

    % global lab-frame movie animation of the walk bout (actual image + APT tips)
    giffile = fullfile(plotdir, sprintf('walk_ex%d_%s.gif', e, matlab.lang.makeValidName(ex(e).label)));
    render_walk_gif_global(readframe, trx_s, aptdata, groundcontact, fly, walk_t0, walk_t1, giffile);
    fprintf('  saved %s\n', giffile);
end

if exist('movie_fid','var') && movie_fid > 0, fclose(movie_fid); end

fprintf('\n===== TOTAL indep-vs-code mismatches across all examples: %d =====\n', nmismatch_total);
fprintf('Figures in %s\n', plotdir);

% ---------------------------------------------------------------------------
function [lag, refframes, matchframes] = indep_p2a_lag(bd, ref_limb, ant_limb, walk_t0, walk_t1, ts, ref_event, pairing)
% Independent re-implementation of the P2A pairing (no findClosestSteps),
% branching on pairing to match computeP2ALag:
%   'forward' : first anterior onset in [ref, next_ref); lag >= 0 (n-1 events)
%   'signed'  : nearest anterior onset in [ref-0.5P, ref+0.5P], signed (n events)
refb = bd.perlimb(ref_limb).(ref_event);
idx = find(refb.start_indices >= walk_t0 & refb.end_indices <= walk_t1);
ref_on = sort(refb.start_indices(idx));
ant_on = sort(bd.perlimb(ant_limb).swing.start_indices(:)');
nts = numel(ts);
n = numel(ref_on);

if strcmp(pairing, 'forward')
    lag = nan(1, max(0,n-1)); refframes = nan(1, max(0,n-1)); matchframes = nan(1, max(0,n-1));
    for s = 1:n-1
        w0 = ref_on(s); w1 = ref_on(s+1);   % window [w0, w1) = one reference period
        cand = ant_on(ant_on >= w0 & ant_on < w1);
        refframes(s) = w0;
        if ~isempty(cand)
            mfr = cand(1);                  % closest to w0 (all candidates >= w0)
            matchframes(s) = mfr;
            if mfr >= 1 && mfr <= nts && w0 >= 1 && w0 <= nts
                lag(s) = (ts(mfr) - ts(w0)) * 1000;
            end
        end
    end
else  % 'signed'
    lag = nan(1,n); refframes = nan(1,n); matchframes = nan(1,n);
    if n < 2, return; end
    halfwin = 0.5 * median(diff(ref_on));
    for s = 1:n
        r = ref_on(s); refframes(s) = r;
        cand = ant_on(ant_on >= r - halfwin & ant_on <= r + halfwin);
        if isempty(cand), continue; end
        [~, k] = min(abs(cand - r));        % nearest (prev or next)
        mfr = cand(k); matchframes(s) = mfr;
        if mfr >= 1 && mfr <= nts && r >= 1 && r <= nts
            lag(s) = (ts(mfr) - ts(r)) * 1000;   % signed
        end
    end
end
end

function s = fnum(x)
if isnan(x), s = 'NaN'; else, s = sprintf('%.4g', x); end
end

function out = tern(c,a,b)
if c, out = a; else, out = b; end
end

function render_walk_gif_global(readframe, trx_s, aptdata, gc, fly, t0, t1, outfile)
% Global lab-frame animation over a walk bout: the actual movie image in a
% FIXED window covering the fly's path, with APT leg tips overlaid
% (orange = stance, blue = swing) and a yellow body-center trajectory trail,
% so global translation and turning are visible.
legtip_kps = [12 13 14 15 16 17];          % RF RM RH LH LM LF -> limbs 1..6
ff = trx_s(fly).firstframe;
nT = trx_s(fly).endframe - ff + 1;
gcf_ = gc{fly}; if size(gcf_,1) ~= 6 && size(gcf_,2) == 6, gcf_ = gcf_'; end
t1 = min([t1, nT, size(gcf_,2)]);
tidxs = t0:t1;

% fixed window = bbox of the fly's path over the walk + pad (pixels)
xs = trx_s(fly).x(tidxs); ys = trx_s(fly).y(tidxs);
pad = 70;
win_x = [floor(min(xs))-pad, ceil(max(xs))+pad];
win_y = [floor(min(ys))-pad, ceil(max(ys))+pad];

apt = aptdata.pTrk{fly};                   % 21 x 2 x T
% fixed grayscale range (from first frame) so the background doesn't flicker
im0 = readframe(ff + t0 - 1);
clim = double([min(im0(:)) max(im0(:))]); if clim(1) >= clim(2), clim = [0 255]; end
fh = figure('Position',[100 100 540 560],'Color','w');
ax = axes(fh);
for ii = 1:numel(tidxs)
    tidx = tidxs(ii);
    im = readframe(ff + tidx - 1);
    cla(ax);
    imagesc(ax, im); colormap(ax, gray); set(ax,'CLim',clim); hold(ax,'on');
    plot(ax, xs(1:ii), ys(1:ii), '-', 'Color',[1 1 0], 'LineWidth',1.5);     % path trail
    plot(ax, xs(ii), ys(ii), 'y+', 'MarkerSize',10, 'LineWidth',1.5);
    for L = 1:6
        kp = legtip_kps(L);
        kx = apt(kp,1,tidx); ky = apt(kp,2,tidx);
        if isnan(kx), continue; end
        if gcf_(L,tidx) == 1, c = [1 0.55 0.1]; else, c = [0.2 0.5 1]; end
        plot(ax, kx, ky, 'o','MarkerFaceColor',c,'MarkerEdgeColor','k','MarkerSize',4);
    end
    xlim(ax, win_x); ylim(ax, win_y); set(ax,'DataAspectRatio',[1 1 1],'XTick',[],'YTick',[]);
    title(ax, sprintf('frame %d   (orange=stance, blue=swing, yellow=path)', ff+tidx-1),'FontSize',10);
    drawnow;
    [A,map] = rgb2ind(frame2im(getframe(fh)),128);
    if ii == 1
        imwrite(A,map,outfile,'gif','LoopCount',Inf,'DelayTime',0.06);
    else
        imwrite(A,map,outfile,'gif','WriteMode','append','DelayTime',0.06);
    end
end
close(fh);
end

%% setpath
modpath

%% set parameters
settingsdir = '/groups/branson/home/robiea/Code_versioned/BransonFlyDiscoSettings/settings';
% screen-matched dated protocols with do_onfloor_filtering,1 + save_swingstancebouts,1
proto_by_screen = containers.Map( ...
    {'VNC','VNC2','VNC3'}, ...
    {'20260622_flybubble_LED_VNC','20260622_flybubble_LED_VNC2','20260622_flybubble_LED_VNC3'});

%% expdir
% balanced YNA_K_162984 control sample (2 per screen x rig) for gait-vs-speed
L = load('/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20260619_velmag_speedbins/plots/gaitspeed_explist.mat');
explist = L.explist;
expscreen = L.expscreen;

%% run code
% force a full recompute so onfloor+save run cleanly (also clears read-only
% symlinked perframe files). Rerun ALL 24 since prior runs were non-onfloor.
for i = 1:numel(explist)
    expdir = explist{i};
    proto = proto_by_screen(expscreen{i});
    fprintf('\n=== [%d/%d] %s (%s) ===\n', i, numel(explist), expdir, proto);
    params = {'settingsdir',settingsdir,'analysis_protocol',proto,'forcecompute',true,'debug',false};
    try
        FlyDiscoComputeLocomotionMetrics(expdir,params{:})
    catch ME
        fprintf(2,'FAILED: %s\n', ME.message);
    end
end
fprintf('\nALL DONE\n');

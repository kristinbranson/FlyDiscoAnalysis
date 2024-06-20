function FlyDiscoComputeAptPerFrameFeatures(expdir, varargin)
% This function implements the computeaptperframefeatures stage.

% varargin will be a sequence of key-value pairs, including at least the keys
% 'settingsdir', 'analysis_protocol', and 'forcecompute'.  It will also
% contain any additional stage-specific key-value pairs passed as the last
% argument to FlyDiscoPipelineStage().


version = '0.1';
% addpath /groups/branson/home/robiea/Code_versioned/locomotion_analysis/Alice
% Parse the optional arguments
[settingsdir, analysis_protocol,datalocparamsfilestr,forcecompute] = ...
    myparse(varargin,...
    'settingsdir', default_settings_folder_path(), ...
    'analysis_protocol', 'current', ...
    'datalocparamsfilestr','dataloc_params.txt',...
    'forcecompute', false) ;

% starting stage messgage
logfid = 1;
timestamp = datestr(now,'yyyymmddTHHMMSS');
real_analysis_protocol = GetRealAnalysisProtocol(analysis_protocol,settingsdir);

fprintf(logfid,'\n\n***\nRunning FlyDiscoComputeAptPerFrameFeatures version %s analysis_protocol %s (linked to %s) at %s\n',version,analysis_protocol,real_analysis_protocol,timestamp);



% make list of special perframe features being computed
aptpfflist = {'nfeet_ground'};

% Initialize trx class object
fprintf('Initializing trx...\n');
trx = FBATrx('analysis_protocol',analysis_protocol,'settingsdir',settingsdir,...
    'datalocparamsfilestr',datalocparamsfilestr);
trx.AddExpDir(expdir,'dooverwrite',false,'openmovie',false);

% if force compute is true, delete aptPFF if they exist
%should check if files exist first
% TO DO ADD  tips_velmag to list for deleting if forcecompute 
if forcecompute,
    for i = 1:numel(aptpfflist)
        curr_aptpff = fullfile(expdir,trx.dataloc_params.perframedir,[aptpfflist{i} '.mat']);
        if exist(curr_aptpff,'file')
            fprintf('Deleting per-frame data file %s\n',aptpfflist{i});
            delete(curr_aptpff);
        end
    end
end



% read in parameters
aptperframefeaturesparamsfile = fullfile(trx.settingsdir,trx.analysis_protocol,trx.dataloc_params.aptperframefeaturesparamsfilestr);
aptperframefeatures_params = ReadParams(aptperframefeaturesparamsfile);


% get leg tip velocities

% % specify leg tip points
legtip_landmarknums = aptperframefeatures_params.legtip_landmarknums;

% if perframe features exist load from them

% load velmag of leg tips from JAABA apt perframe features
% transform to format ground contact expects
% nfly x 6 leg tips x frames cell array
perframestr = 'apt_view1_global_velmag_';
tips_velmag =cell(1,trx.nflies);
perframeload_success = ones(1,numel(legtip_landmarknums));

for i = 1:numel(legtip_landmarknums)
    fn = [perframestr,num2str(legtip_landmarknums(i))];
    if exist(fullfile(expdir,trx.dataloc_params.perframedir,[fn,'.mat']),'file') % currently trx class can't compute apt features on the fly
        for fly = 1:trx.nflies
            tips_velmag{fly}(i,:) = trx(fly).(fn);
        end
    else
        perframeload_success(i) = 0;
    end
end

if ~all(perframeload_success) && exist(fullfile(expdir,"tips_velmag.mat"),'file') 

    % if tips_velmag already exists load them
    load(fullfile(expdir,"tips_velmag.mat"),'tips_velmag');

else
    % compute and save if doesn't already exist
    % need apt data 
    aptfile = trx.dataloc_params.apttrkfilestr;
    aptdata = TrkFile.load(fullfile(expdir,aptfile));
    ts = trx.movie_timestamps{:};
    dt = diff(ts)';
    pxpermm = trx.pxpermm;
  
        % leg tip velocity
    [tips_velmag] = compute_legtipvelmag(aptdata,dt,pxpermm,legtip_landmarknums);

    save(fullfile(expdir,'tips_velmag.mat'),'tips_velmag');
end

% run groundcontact detections
gc_threshold_low = aptperframefeatures_params.gc_threshold_low;
gc_threshold_high = aptperframefeatures_params.gc_threshold_high;
[groundcontact] = compute_groundcontact(tips_velmag,'gc_threshold_low',gc_threshold_low,'gc_threshold_high',gc_threshold_high);

%compute number of feet on the ground
data = cell(1,trx.nflies);
for fly = 1:trx.nflies
    data{fly} = sum(groundcontact{fly},1);
end
units = parseunits('unit');
save(fullfile(expdir,trx.dataloc_params.perframedir,'nfeet_ground.mat'),'data','units');

%compute swing stance bouts
[perfly_limbboutdata] = limbSwingStance(groundcontact);

% filter swing stance for during walking
% how to load score data from trx file??
% ScoreFile = fullfile(expdir,"scores_Walk2.mat");
% load(ScoreFile,'allScores');
% digitalscores = allScores.postprocessed;
[~,digitalscores] = LoadScoresFromFile(trx,'scores_Walk2',1);
[walkingSwingStance] = restrictedSwingStance(digitalscores,perfly_limbboutdata,'trajectory');


% convert bout start and end indices to movie reference frame
[walkingSwingStance] = convert_boutdata_to_movieformat(trx,walkingSwingStance);

% seperate walking stwing stance for stimon and stimoff data
% how to load stim data from trx file??
% indicatorfile = fullfile(expdir,"indicatordata.mat");
% load(indicatorfile,'indicatorLED');
% digitalindicator = indicatorLED.indicatordigital;
indicatordata = trx.getIndicatorLED(1);
digitalindicator = indicatordata.indicatordigital;
[stimON_walkingSwingStance] = restrictedSwingStance(digitalindicator,walkingSwingStance,'movie');
[stimOFF_walkingSwingStance] = restrictedSwingStance(~digitalindicator,walkingSwingStance,'movie');


% save(fullfile(expdir,'swingstance.mat'),'stimON_walkingSwingStance','stimOFF_walkingSwingStance');


% compute durations for bouts
[bout_metrics_ON] = computeboutmetrics2(trx,stimON_walkingSwingStance);
[bout_metrics_OFF] = computeboutmetrics2(trx,stimOFF_walkingSwingStance);



% If something goes wrong, throw an error.  Any values returned from this
% function are ignored.
try
    save(fullfile(expdir,trx.dataloc_params.aptperframeboutstatsfilestr),'bout_metrics_ON','bout_metrics_OFF')
catch ME
    warning('FlyDiscoComputeAptPerFrameFeatures:save',...
        'Could not save information to file %s: %s',filename,getReport(ME));
end

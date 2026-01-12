function [] = LocomotionCombinePerFrameStats(expdir,varargin)

[settingsdir, analysis_protocol,datalocparamsfilestr,forcecompute,~,debug] = ...
    myparse(varargin,...
    'settingsdir', default_settings_folder_path(), ...
    'analysis_protocol', 'current', ...
    'datalocparamsfilestr','dataloc_params.txt',...
    'forcecompute', false, ...
    'do_run', [],...
    'debug', []) ;

 %% parameters
  datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
  dataloc_params = ReadParams(datalocparamsfile);
  locomotiondatafile = fullfile(expdir,dataloc_params.locomotionmetricsswingstanceboutstatsfilestr);
  

if ~exists(locomotiondatafile)
    error(' locomotion data file %s does nott exist',locomotiondatafile)
end

load(locomotiondatafile);

%% create locomotionstats_perframefeatures from my data structs in each expdir
% need a experiment level mean, std, and Z for each stat. 

% statsperexp
% struct with field names = each stat
% meanmean
% stdmean
% Z

statsperexp = struct;
% bout_metrics_ON (leaving OFF aside for now). 
flies = 'allflies';
limbs = {'all_limbs','perlimb','pairs'};
state = {'swing','stance'};
ignorelist = {'fly','start_indices','end_indices','Nboutspervelmagbin','velmagbincenters','meanboutdurationsofvelmagbins','stdboutdurationsofvelmagbins'};


for l = 1:numel(limbs)
    numlimb = numel(bout_metrics_ON.(flies).(limbs{l}));

    for ll = 1:numlimb
        if numlimb ==1
            limbname = 'all';
        elseif numlimb ==3
            limbname = sprintf('pair%s',num2str(ll));
        elseif numlimb== 6
            limbname = sprintf('limb%s',num2str(ll));
        end

        for s = 1:numel(state)
            featurefields = fields(bout_metrics_ON.(flies).(limbs{l})(ll).(state{s}));
            for f = 1:numel(featurefields)
                if ~ismember(featurefields{f},ignorelist) && ~isstruct(bout_metrics_ON.(flies).(limbs{l})(ll).(state{s}).(featurefields{f}))
                    % Note mean = mean (data for all bouts in experiments)
                    funname = sprintf('%s__%s__LEDon__%s',featurefields{f},state{s},limbname);
                    currstruct = struct;
                    currdata = bout_metrics_ON.(flies).(limbs{l})(ll).(state{s}).(featurefields{f});
                    currstruct.mean = mean(currdata);
                    currstruct.std = std(currdata);
                    currstruct.Z = numel(currdata); % number of bouts because 1 data per bout
                    statsperexp.(funname) = currstruct;
                elseif ~ismember(featurefields{f},ignorelist) && isstruct(bout_metrics_ON.(flies).(limbs{l})(ll).(state{s}).(featurefields{f}))
                    % Note: mean = mean (means of perframe features for
                    % each bout in the experiemnt). NOT per fly means
                    % Z = number of bouts
                    funname = sprintf('%s__%s__LEDon__%s',featurefields{f},state{s},limbname);
                    currstruct = struct;
                    currdata = bout_metrics_ON.(flies).(limbs{l})(ll).(state{s}).(featurefields{f}).mean;
                    currstruct.mean = mean(currdata);
                    currstruct.std = std(currdata);
                    currstruct.nbouts = numel(currdata);
                    currstruct.Z = sum(bout_metrics_ON.(flies).(limbs{l})(ll).(state{s}).(featurefields{f}).n);
                    statsperexp.(funname) = currstruct;
                end
            end
        end
    end
end

% walk_metrics_ON






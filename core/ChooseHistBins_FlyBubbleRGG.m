% Computes general histogram bin edges and centers
% for use by PlotPerFrameHistsGroup.
% Saves the output file, hist_perframebins.mat, to the analysis_protocol folder. 

 
%% set up path
modpath
%% USER defined parameters
% User need to modify:
% 1. settingsdir, analysis protocol, control line name
% 2. behaviornames list
% 3. full metadata data struct (either load metadata structure or pull parameters for
% getExperimentDirsFlyDisco) 

settingsdir = '/groups/branson/home/robiea/Code_versioned/FlyDiscoAnalysis/settings-internal';

analysis_protocol = '20240402_flybubble_LED_VNC2';
datalocparamsfilestr = 'dataloc_params.txt';
params = {'analysis_protocol',analysis_protocol,...
  'datalocparamsfilestr',datalocparamsfilestr,...
  'settingsdir',settingsdir};

control_line = 'YNA_K_162984';
% any behavior defined as a frame condition in statframeconditions.txt and
% used with the special case 'duration' in hist_perframefeatures.txt
behaviornames = {'walk'};


%% data pull or load experimental data metadata file

savefile = ['/groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/dataforHistBins_',datestr(now,'yyyymmdd')];
if exist(savefile,file)
    load savefile
    experiments_all = expdirstruct;
else
rootdatadir = '/groups/branson/bransonlab/flydisco_data';
[expdirstruct] = getExperimentDirsFlyDisco(rootdatadir,'metadatafile','Metadata.xml','expdirname','VNC2*','autocheckin',true,'autocheckcomplete',true,'manualcheck',true,'TrajNum',false);
% filter out fails
pidx = find([expdirstruct.automated_pf] == 'P');
mpidx = find([expdirstruct.manual_fail] =='U');
idx = pidx & mpidx;

ogexpdir = expdirstruct;
expdirstruct = expdirstruct(idx);

save(savefile,'expdirstruct');
experiments_all = expdirstruct;
end

% load('/groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/expdirs_VNC_20240711_VNC_passAandM_metadata.mat');
% experiments_all = expdirstruct;

%% choose experiments to histogram

weight_order = {'genotype','date','rig','bowl'};

% % 
ngal4 = 100;
ncontrols = 100;

[experiments,ngal4_chosen,ncontrols_chosen] = choose_expdirs(experiments_all,ngal4,ncontrols,'weight_order',weight_order,'control_line',control_line);

expdirs = {experiments.file_system_path};
nexpdirs = numel(expdirs);

%% parameters
datalocparamsfile = fullfile(settingsdir,analysis_protocol,'dataloc_params.txt');
dataloc_params = ReadParams(datalocparamsfile);

%% read parameters

% perframefnsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.perframefnsfilestr);
% perframefns = importdata(perframefnsfile);
prctile_store = 10;

aicparamsfile = fullfile(settingdir,analysis_protocol,dataloc_params.automaticchecksincomingparamsfilestr);
aic_parameters = ReadParms(aicparamsfile);
fps = aic_parameters.fps;


histparamsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.histparamsfilestr);

histperframefeaturesfile = fullfile(settingsdir,analysis_protocol,dataloc_params.histperframefeaturesfilestr);
hist_perframefeatures = ReadHistPerFrameFeatures2(histperframefeaturesfile);
hist_params = ReadParams(histparamsfile);

statframeconditionsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.statframeconditionfilestr);
frameconditiondict = ReadParams(statframeconditionsfile);
statflyconditionsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.statflyconditionfilestr);
flyconditiondict = ReadParams(statflyconditionsfile);



%% loop over per-frame parameters

% data structure we will store everything in
bins = struct;
perframefns = setdiff({hist_perframefeatures.field},{'duration'});

% add in behavior labels
for i = 1:numel(behaviornames),
  perframefns{end+1} = ['duration_',behaviornames{i}];
end
  
for fni = 1:numel(perframefns),
  
  fn = perframefns{fni};
  
  if strcmp(fn,'sex'),
    continue;
  end
  
  fprintf('Choosing bins for %s...\n',fn);
  bins.(fn) = struct;
  
  % store smallest and largest data to get lims for histogramming
  datasmall = [];
  datalarge = [];

  fprintf('Choosing lims...\n');
  % loop through experiments
  for expi = 1:numel(expdirs),
    
    expdir = expdirs{expi};
    %fprintf('exp(%d) = %s\n',expi,expdir);

    m = regexp(fn,'^duration_(.+)$','tokens','once');
    if ~isempty(m),
      m = m{1};
      s = regexp(frameconditiondict.(m){1},'^[^_]+_(.+)$','tokens','once');
      if isempty(s),
        error('Could not parse %s',m);
      end
      filename = fullfile(expdir,[s{1},'.mat']);
      if exist(filename,'file'),
        labeldata = load(filename);
        datacurr = (cell2mat(labeldata.t1s)-cell2mat(labeldata.t0s))/fps;
      else
        s{1} = regexprep(s{1},'label','score');
        filename = fullfile(expdir,[s{1},'.mat']);
        if ~exist(filename,'file'),
          error('Could not find score or label %s for %s',filename,m);
        end
        scoredata = load(filename);
        if isfield(scoredata.allScores,'postprocessed'),
          datacurr = [];
          for j = 1:numel(scoredata.allScores.postprocessed),
            scoredata.allScores.postprocessed{j}(isnan(scoredata.allScores.postprocessed{j})) = 0;
            [t0s,t1s] = get_interval_ends(scoredata.allScores.postprocessed{j});
            datacurr = [datacurr,(t1s-t0s)/fps]; %#ok<AGROW>
          end
        else
          datacurr = (cell2mat(scoredata.allScores.t1s)-cell2mat(scoredata.allScores.t0s))/fps;
        end
      end
    else
      filename = fullfile(expdir,dataloc_params.perframedir,[fn,'.mat']);
      if ~exist(filename,'file'), 
        warning('Could not find file %s, skipping',filename);
        continue;
      end
      % get per-frame data for these flies
      tmp = load(filename,'data');
      if ~iscell(tmp.data),
        warning('Data in %s is not a cell, skipping.',filename);
        continue;
      end
      datacurr = cell2mat(tmp.data);
    end
    thresh = prctile(datacurr,[prctile_store,100-prctile_store]);
    datasmall = [datasmall,datacurr(datacurr <= thresh(1))]; %#ok<AGROW>
    datalarge = [datalarge,datacurr(datacurr >= thresh(2))]; %#ok<AGROW>
    
  end
  
  % compute actual lim percentiles from the stored boundary data
  bins.(fn).lims = [nan,nan];

  mindata = min(datasmall);
  maxdata = max(datalarge);
  % set lower bound to 0 if this makes sense
  if mindata >= 0 && mindata / maxdata < 1e-3,
    bins.(fn).lims(1) = 0;
    fprintf('Choosing lb = 0 for %s\n',fn);
  % set lower bound to -pi if this makes sense
  elseif maxdata > -pi && abs(mindata + pi) / (maxdata + pi) < 1e-3,
    bins.(fn).lims(1) = -pi;
    fprintf('Choosing lb = -pi for %s\n',fn);
  % set lower bound to -pi/2 if this makes sense
  elseif maxdata > -pi/2 && abs(mindata + pi/2) / (maxdata + pi/2) < 1e-3,
    bins.(fn).lims(1) = -pi/2;
    fprintf('Choosing lb = -pi/2 for %s\n',fn);
  % set lower bound to -1 if this makes sense
  elseif maxdata > -1 && abs(mindata + 1) / (maxdata + 1) < 1e-3,
    bins.(fn).lims(1) = -1;
    fprintf('Choosing lb = -1 for %s\n',fn);
  else
    bins.(fn).lims(1) = prctile(datasmall,hist_params.lim_prctile(1)*prctile_store);
    fprintf('Choosing lb from percentile for %s\n',fn);
  end
  
  % set upper bound to pi if this makes sense
  if mindata < pi && abs(maxdata - pi) / (pi - mindata) < 1e-3,
    bins.(fn).lims(2) = pi;
    fprintf('Choosing ub = pi for %s\n',fn);
  % set upper bound to pi/2 if this makes sense
  elseif mindata < pi/2 && abs(maxdata - pi/2) / (pi/2 - mindata) < 1e-3,
    bins.(fn).lims(2) = pi/2;
    fprintf('Choosing ub = pi/2 for %s\n',fn);
  % set upper bound to 1 if this makes sense
  elseif mindata < 1 && abs(maxdata - 1) / (1 - mindata) < 1e-3,
    bins.(fn).lims(2) = 1;
    fprintf('Choosing ub = 1 for %s\n',fn);
  else    
    bins.(fn).lims(2) = prctile(datalarge,100-(100-hist_params.lim_prctile(2))*prctile_store);
    fprintf('Choosing ub from percentile for %s\n',fn);
  end

  % edges based on linear spacing
  [bins.(fn).edges_linear,bins.(fn).centers_linear] = ...
    SelectHistEdges(hist_params.nbins,bins.(fn).lims,'linear');
  
  % edges based on log spacing
  if mindata < 0,
    [bins.(fn).edges_log,bins.(fn).centers_log] = ...
      SelectHistEdges(hist_params.nbins,bins.(fn).lims,'logabs');
  else
    [bins.(fn).edges_log,bins.(fn).centers_log] = ...
      SelectHistEdges(hist_params.nbins,bins.(fn).lims,'log');
  end
  
end

%% save results

% savename = fullfile(settingsdir,analysis_protocol,dataloc_params.histperframebinsfilestr);
savename = fullfile(settingsdir,analysis_protocol,'histbins_20250401');
save(savename,'bins');
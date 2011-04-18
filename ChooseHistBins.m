
if ispc,
  addpath E:\Code\JCtrax\misc;
  addpath E:\Code\JCtrax\filehandling;
else
  addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
  addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
end

%% data locations

if ispc,
  settingsdir = 'E:\Code\FlyBowlAnalysis\settings';
  rootdir = 'E:\Data\FlyBowl\CtraxTest20110407';
else
  settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
  %rootdir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data';
  rootdir = '/groups/branson/bransonlab/tracking_data/olympiad/FlyBowl/CtraxTest20110407';
end

expdirs = dir(fullfile(rootdir,'*_*'));
expdirs = cellfun(@(s) fullfile(rootdir,s),{expdirs([expdirs.isdir]).name},'UniformOutput',false);

%% parameters

analysis_protocol = '20110407';
params = {'settingsdir',settingsdir,...
  'analysis_protocol',analysis_protocol};
datalocparamsfile = fullfile(settingsdir,analysis_protocol,'dataloc_params.txt');
dataloc_params = ReadParams(datalocparamsfile);

%% load data

% obj = Trx('settingsdir',settingsdir,'analysis_protocol',analysis_protocol);
% for i = 1:numel(expdirs),
%   obj.AddExpDir(expdirs{i});
% end

%% read parameters

perframefnsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.perframefnsfilestr);
perframefns = importdata(perframefnsfile);
prctile_store = 10;

histparamsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.histparamsfilestr);
hist_params = ReadParams(histparamsfile);

%% loop over per-frame parameters

% data structure we will store everything in
bins = struct;


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
    fprintf('exp(%d) = %s\n',expi,expdir);

    filename = fullfile(expdir,dataloc_params.perframedir,[fn,'.mat']);
    
    % get per-frame data for these flies
    tmp = load(filename,'data');
    datacurr = cell2mat(tmp.data);
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

savename = fullfile(settingsdir,analysis_protocol,dataloc_params.histperframebinsfilestr);
save(savename,'bins');
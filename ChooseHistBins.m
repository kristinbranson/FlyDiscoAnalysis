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
  analysis_protocol = '20110202_pc';
  dataloc_params = ReadParams(fullfile(settingsdir,analysis_protocol,'dataloc_params.txt'));
  expdir = 'pBDPGAL4U_TrpA_Rig1Plate10BowlA_20110202T105734';
  expdir_read = fullfile(dataloc_params.rootreaddir,expdir);
  
else
  settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
  protocol = '20110208';
  analysis_protocol = '20110208';
  [expdirs,expdir_reads,expdir_writes,experiments,rootreaddir,rootwritedir] = ...
    getExperimentDirs('protocol',protocol,'subreadfiles',{'perframe'});
end

%% number of experiments to analyze
maxnexpdirs_control = 10;
maxnexpdirs_gal4 = 10;
idx = cellfun(@isempty,regexp(expdirs,'^pBD.*$','once'));
ngal4 = nnz(idx);
ncontrol = nnz(~idx);
if ngal4 < maxnexpdirs_gal4,
  ncontrol = min(ncontrol,maxnexpdirs_control+maxnexpdirs_gal4-ngal4);
else
  ncontrol = maxnexpdirs_control;
end
if ncontrol < maxnexpdirs_control,
  ngal4 = min(ngal4,maxnexpdirs_gal4+maxnexpdirs_control-ncontrol);
else
  ngal4 = maxnexpdirs_gal4;
end

%% choose different experiments

% distance between experiments
D = zeros(numel(expdirs));
for i = 1:numel(expdirs)-1,
  exp1 = parseExpDir(expdirs{i});
  for j = i+1:numel(expdirs),
    exp2 = parseExpDir(expdirs{j});
    D(i,j) = double(~strcmp(expdirs{i},expdirs{j}))*100000 + ...
      double(~strcmp(exp1.line,exp2.line))*10000 + ...
      double(~strcmp(exp1.rig,exp2.rig))*1000 + ...
      double(~strcmp(exp1.bowl,exp2.bowl))*100 + ...
      double(~strcmp(exp1.date(1:8),exp2.date(1:8)))*10 + ...
      abs(str2double(exp1.date(10:11))-str2double(exp2.date(10:11)))*1;
    D(j,i) = D(i,j);
  end
end


% choose experiments with large average distance to experiments already
% chosen
Dgal4 = D(idx,idx);
idxgal4 = [1,nan(1,ngal4-1)];
sumd = Dgal4(idxgal4(1),:);
for i = 2:ngal4,
  [~,idxgal4(i)] = max(sumd);
  sumd = sumd+Dgal4(idxgal4(i),:);
end
idxgal41 = find(idx);
idxgal4 = sort(idxgal41(idxgal4));

Dcontrol = D(~idx,~idx);
idxcontrol = [1,nan(1,ncontrol-1)];
sumd = Dcontrol(idxcontrol(1),:);
for i = 2:ncontrol,
  [~,idxcontrol(i)] = max(sumd);
  sumd = sumd+Dcontrol(idxcontrol(i),:);
end
idxcontrol1 = find(~idx);
idxcontrol = sort(idxcontrol1(idxcontrol));

expdirs_gal4 = expdirs(idxgal4);
expdirs_control = expdirs(idxcontrol);
expdir_reads_gal4 = expdir_reads(idxgal4);
expdir_reads_control = expdir_reads(idxcontrol);
expdirs = [expdirs_control,expdirs_gal4];
expdir_reads = [expdir_reads_control,expdir_reads_gal4];

%% load data

obj = Trx('settingsdir',settingsdir,'analysis_protocol',analysis_protocol);
for i = 1:numel(expdir_reads),
  obj.AddExpDir(expdir_reads{i});
end

%% read parameters

perframefnsfile = fullfile(obj.settingsdir,obj.analysis_protocol,obj.dataloc_params.perframefnsfilestr);
perframefns = importdata(perframefnsfile);
prctile_store = 10;

histparamsfile = fullfile(obj.settingsdir,obj.analysis_protocol,obj.dataloc_params.histparamsfilestr);
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
  for expi = 1:obj.nexpdirs,
    
    fprintf('exp(%d) = %s\n',expi,obj.expdirs{expi});
    
    flies = obj.exp2flies{expi};
    
    % get per-frame data for these flies
    datacurr = cell2mat(obj(flies).(fn));
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

savename = fullfile(obj.settingsdir,obj.analysis_protocol,obj.dataloc_params.histperframebinsfilestr);
save(savename,'bins');
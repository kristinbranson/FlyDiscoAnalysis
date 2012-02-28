function [statsperfly,statsperexp] = ComputeHackHitCategoryStatistics_KB(expdir,varargin)

[analysis_protocol,settingsdir,datalocparamsfilestr,rootoutputdir] = ...
    myparse(varargin,...
    'analysis_protocol','current',...
    'settingsdir','/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings',...
    'datalocparamsfilestr','dataloc_params.txt',...
    'rootoutputdir',0);

%% read in the data locations
datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
dataloc_params = ReadParams(datalocparamsfile);

perframefns = {'velmag','dcenter','dnose2ell','dnose2ell',...
  'x','y','theta','a','b','timestamps',...
  'x_mm','y_mm','a_mm','b_mm','theta_mm','dt','sex'};

%% create a local version of the directory if we don't have write permissions
if ischar(rootoutputdir),
  
  % create the root directory
  if ~exist(rootoutputdir,'dir'),
    [success,msg] = mkdir(rootoutputdir);
    if ~success,
      error('Could not create root output directory %s: %s',rootoutputdir,msg);
    end
  end
  
  % create the experiment directory
  [~,basename] = fileparts(expdir);
  outexpdir = fullfile(rootoutputdir,basename);
  if ~exist(outexpdir,'dir'),
    [success,msg] = mkdir(outexpdir);
    if ~success,
      error('Could not create output directory %s: %s',outexpdir,msg);
    end
  end

  % soft link to the trx
  outtrxfilename = fullfile(outexpdir,dataloc_params.trxfilestr);
  if ~exist(outtrxfilename,'file'),
    intrxfilename = fullfile(expdir,dataloc_params.trxfilestr);
    if ~exist(intrxfilename,'file'),
      error('Input trx file %s does not exist',intrxfilename);
    end
  
    cmd = sprintf('ln -s %s %s',intrxfilename,outtrxfilename);
    system(cmd);
    
    if ~exist(outtrxfilename,'file'),
      error('Soft-linked trx file %s not created successfully',outtrxfilename);
    end
  end

  % soft link to the movie
  outmoviefilename = fullfile(outexpdir,dataloc_params.moviefilestr);
  if ~exist(outmoviefilename,'file'),
    inmoviefilename = fullfile(expdir,dataloc_params.moviefilestr);
    if ~exist(inmoviefilename,'file'),
      error('Input movie file %s does not exist',inmoviefilename);
    end
  
    cmd = sprintf('ln -s %s %s',inmoviefilename,outmoviefilename);
    system(cmd);
    
    if ~exist(outmoviefilename,'file'),
      error('Soft-linked movie file %s not created successfully',outmoviefilename);
    end
  end

  
  % create the per-frame directory
  outperframedir = fullfile(outexpdir,dataloc_params.perframedir);
  % soft links to everything in perframedir
  inperframedir = fullfile(expdir,dataloc_params.perframedir);
  if exist(inperframedir,'dir') && ~exist(outperframedir,'dir'),
    % do all the per-frame files we need exist?
    allthere = true;
    for i = 1:numel(perframefns),
      fn = fullfile(inperframedir,[perframefns{i},'.mat']);
      if ~exist(fn,'file'),
        allthere = false;
        break;
      end
    end
    % if they do, then just make a link to the perframedir
    if allthere,
      cmd = sprintf('ln -s %s %s',inperframedir,outperframedir);
      system(cmd);
      if ~exist(outperframedir,'dir'),
        error('Creating soft link to perframedir %s did not work',inperframedir);
      end
    end
  else
    allthere = false;
  end
  
  if ~allthere,
    if ~exist(outperframedir,'dir'),
      [success,msg] = mkdir(outperframedir);
      if ~success,
        error('Could not create output per-frame directory %s: %s',outperframedir,msg);
      end
    end
      
    fns = dir(fullfile(inperframedir,'*.mat'));
    fns = {fns.name};
    %perframefiles = cellfun(@(s) [s,'.mat'],perframefns,'UniformOutput',false);
    %fns_missing = setdiff(perframefiles,fns);
    for i = 1:numel(fns);
      infn = fullfile(inperframedir,fns{i});
      outfn = fullfile(outperframedir,fns{i});
      if exist(outfn,'file'),
        continue;
      end
      cmd = sprintf('ln -s %s %s',infn,outfn);
      system(cmd);
      if ~exist(outfn,'file'),
        warning('Did not successfully soft link to per-frame file %s',infn);
      end
    end
  end

else
  outexpdir = expdir;
end
  

%% read in parameters
stats_params = ReadParams(fullfile(settingsdir,analysis_protocol,dataloc_params.hackhitcategoryparamsfilestr));

%% create trx variable

trx = Trx('analysis_protocol',analysis_protocol,'settingsdir',settingsdir,...
    'datalocparamsfilestr',datalocparamsfilestr);

fprintf('Loading trajectories for %s...\n',expdir);

trx.AddExpDir(outexpdir,'dooverwrite',false);

%% initialize structures

statsperfly = struct;
statsperexp = struct;
%% compute fraction of time "jumping"

statsperfly.isjumping = struct;
statsperflycurr = struct;
behaviorlabel = 'jump';

for i = 1:trx.nflies,
    if trx(i).nframes < stats_params.min_nframes_jump,
        continue;
    end
    isjumping = double(trx(i).velmag >= stats_params.min_velmag_jump);
    doanalyze = true(size(isjumping));
    
    [statsperflycurr.Z(i),statsperflycurr.mean(i),...
        statsperflycurr.std(i),statsperflycurr.prctiles(:,i)] = ...
        ComputePerFrameStats(isjumping,doanalyze,...
        'prctiles_compute',[]);
    statsperflycurr.fracframesanalyzed(i) = nnz(doanalyze) / numel(doanalyze);
    if ~logical(nnz(isjumping)),
        continue;
    end
    [flies_all(i), t0s_all{i}, t1s_all{i}, names_all{i}, off_all(i)] = ...
        createlabelfile(behaviorlabel,isjumping,trx.firstframes(i),i);
end
statsperexp.isjumping = CombinePerFrameStats(statsperflycurr);
statsperfly.isjumping = statsperflycurr;
if exist('flies_all','var')
    idx = logical(flies_all);
    flies = flies_all(idx)';
    t0s = t0s_all(~cellfun('isempty', t0s_all));
    t1s = t1s_all(~cellfun('isempty', t1s_all));
    names = names_all(~cellfun('isempty', names_all));
    idx = logical(off_all);
    off = off_all(idx);
    timestamp = now;
    statsmatsavename = fullfile(outexpdir,'jump_labels.mat');
    save(statsmatsavename,'flies','t0s','t1s','names','off','timestamp');
end
%% compute mean speed while moving
%%%%% uses the same min parameter for moving as fraction of time moving %%%%

statsperfly.meanspeedmoving = struct;
statsperflycurr = struct;
for i = 1:trx.nflies,
    if trx(i).nframes < stats_params.min_nframes_meanspeedmoving,
        continue;
    end
    velmag = trx(i).velmag;
    doanalyze = velmag >= stats_params.min_velmag_moving & ...
        velmag  <= stats_params.max_velmag_moving;
    [statsperflycurr.Z(i),statsperflycurr.mean(i),...
        statsperflycurr.std(i),statsperflycurr.prctiles(:,i)] = ...
        ComputePerFrameStats(velmag,doanalyze);
    statsperflycurr.fracframesanalyzed(i) = 1;
end
statsperexp.meanspeedmoving = CombinePerFrameStats(statsperflycurr);
statsperfly.meanspeedmoving = statsperflycurr;

%% compute fraction of time moving
%%%%% uses the same min parameter for moving as mean speed while moving %%%%
statsperfly.ismoving = struct;
statsperflycurr = struct;
behaviorlabel = 'moving';
for i = 1:trx.nflies,
    if trx(i).nframes < stats_params.min_nframes_meanspeedmoving,
        continue;
    end
    velmag = trx(i).velmag;
    ismoving = double(velmag >= stats_params.min_velmag_moving);
    doanalyze = true(size(velmag));
    [statsperflycurr.Z(i),statsperflycurr.mean(i),...
        statsperflycurr.std(i),statsperflycurr.prctiles(:,i)] = ...
        ComputePerFrameStats(ismoving,doanalyze,...
        'prctiles_compute',[]);
    statsperflycurr.fracframesanalyzed(i) = nnz(doanalyze) / numel(doanalyze);
    if ~logical(nnz(ismoving)),
        continue;
    end
    [flies_all(i), t0s_all{i}, t1s_all{i}, names_all{i}, off_all(i)] = ...
        createlabelfile(behaviorlabel,ismoving,trx.firstframes(i),i);
end
statsperexp.ismoving = CombinePerFrameStats(statsperflycurr);
statsperfly.ismoving = statsperflycurr;
if exist('flies_all','var')
    idx = logical(flies_all);
    flies = flies_all(idx)';
    t0s = t0s_all(~cellfun('isempty', t0s_all));
    t1s = t1s_all(~cellfun('isempty', t1s_all));
    names = names_all(~cellfun('isempty', names_all));
    idx = logical(off_all);
    off = off_all(idx);
    timestamp = now;
    statsmatsavename = fullfile(outexpdir,'moving_labels.mat');
    save(statsmatsavename,'flies','t0s','t1s','names','off','timestamp'); 
end
%% compute fraction of time close together
statsperfly.areclose = struct;
statsperflycurr = struct;
behaviorlabel = 'close';
for i = 1:trx.nflies,
    if trx(i).nframes < stats_params.min_nframes_close,
        continue;
    end
    areclose = double(trx(i).dcenter <= stats_params.max_dcenter_close |...
        trx(i).dnose2ell <= stats_params.max_dnose2ell_close | ...
        trx(i).dell2nose <= stats_params.max_dell2nose_close);
    doanalyze = true(size(areclose));
    [statsperflycurr.Z(i),statsperflycurr.mean(i),...
        statsperflycurr.std(i),statsperflycurr.prctiles(:,i)] = ...
        ComputePerFrameStats(areclose,doanalyze,...
        'prctiles_compute',[]);
    statsperflycurr.fracframesanalyzed(i) = nnz(doanalyze) / numel(doanalyze);
    if ~logical(nnz(areclose)),
        continue;
    end
    [flies_all(i), t0s_all{i}, t1s_all{i}, names_all{i}, off_all(i)] = ...
        createlabelfile(behaviorlabel,areclose,trx.firstframes(i),i);
end
statsperexp.areclose = CombinePerFrameStats(statsperflycurr);
statsperfly.areclose = statsperflycurr;
if exist('flies_all','var')
    idx = logical(flies_all);
    flies = flies_all(idx)';
    t0s = t0s_all(~cellfun('isempty', t0s_all));
    t1s = t1s_all(~cellfun('isempty', t1s_all));
    names = names_all(~cellfun('isempty', names_all));
    idx = logical(off_all);
    off = off_all(idx);
    timestamp = now;
    statsmatsavename = fullfile(outexpdir,'close_labels.mat');
    save(statsmatsavename,'flies','t0s','t1s','names','off','timestamp');
end

%% compute fraction of time chasing
%
% statsperfly.chasing  = struct;
% statsperflycurr = struct;
% for i = 1:trx.nflies,
%   if trx(i).nframes < stats_params.min_nframes_chasing,
%     continue;
%   end
%   velmag = trx(i).velmag;
%   dnose2ell = trx(i).dnose2ell(1:end-1);
%   absanglefrom1to2_nose2ell = trx(i).absanglefrom1to2_nose2ell(1:end-1);
%
%   arechasing = velmag >= stats_params.min_velmag_chasing & ...
%     dnose2ell  <= stats_params.max_dnose2ell_chasing & ...
%     absanglefrom1to2_nose2ell <= stats_params.max_absanglefrom1to2_nose2ell_chasing;
% doanalyze = true(size(arechasing));
%   [statsperflycurr.Z(i),statsperflycurr.mean(i),...
%     statsperflycurr.std(i),statsperflycurr.prctiles(:,i)] = ...
%     ComputePerFrameStats(arechasing,doanalyze, ...
%       'prctiles_compute',[]);
%   statsperflycurr.fracframesanalyzed(i) = nnz(doanalyze) / numel(doanalyze);
% end
% statsperexp.chasing = CombinePerFrameStats(statsperflycurr);
% statsperfly.chasing = statsperflycurr;
%% save results to file

% SaveAllPerFrameStatsTxtFile(statstxtsavename,statsperfly,statsperexp);

% save to mat file
statsmatsavename = fullfile(outexpdir,dataloc_params.hackhitstatsmatfilestr);
save(statsmatsavename,'statsperfly','statsperexp','stats_params');
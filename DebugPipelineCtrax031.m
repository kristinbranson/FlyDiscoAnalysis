% test pipeline on new Ctrax data
% DebugPipelineCtrax031

%% set up paths

addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
addpath /groups/branson/bransonlab/projects/olympiad/SAGE/MATLABInterface/Trunk;
addpath /groups/branson/home/bransonk/tracking/code/lds/hmm;
addpath /groups/branson/home/bransonk/tracking/code/Ctrax/matlab/netlab;
settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
rootdatadir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data';
rootoutputdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlCtrax/20111221/ctrax0.3.1_test';


%% parameters

analysis_protocol = '20120210';
ctrax_protocol = '20111221';
expdirfile = 'expdirs_20120210pipelinetest.txt';

%% experiments

tmp = dir(fullfile(rootoutputdir,'*20*'));
experiment_names = {tmp.name};
expdirs = cellfun(@(x) fullfile(rootoutputdir,x),experiment_names,'UniformOutput',false);
inexpdirs = cellfun(@(x) fullfile(rootdatadir,x),experiment_names,'UniformOutput',false);
nexps = numel(expdirs);

%% soft link to stuff

softlink_files = {'FlyBowlDataCaptureParams_EP*.txt',...
  'movie.ufmf',...
  'Log.txt',...
  'Metadata.*',...
  'QuickStats.png',...
  'QuickStats.txt',...
  'SUCCESS',...
  'temperature.txt',...
  'ufmf_diagnostics.txt',...
  'ufmf_log.txt'};

for i = 1:nexps,
  fprintf('%d / %d: %s\n',i,nexps,experiment_names{i});
  for j = 1:numel(softlink_files),
    infilepattern = fullfile(inexpdirs{i},softlink_files{j});
    tmp = dir(infilepattern);
    for k = 1:numel(tmp),
      infile = fullfile(inexpdirs{i},tmp(k).name);
      outfile = fullfile(expdirs{i},tmp(k).name);
      if ~exist(outfile,'file'),
        cmd = sprintf('ln -s %s %s',infile,outfile);
        unix(cmd);
        if ~exist(outfile,'file'),
          error('File %s does not exist after trying to soft-link',outfile);
        end
      end
    end
  end
end

%% automatic checks incoming

% one experiment
i = 1;
expdir = expdirs{i};
FlyBowlAutomaticChecks_Incoming(expdir,'analysis_protocol',analysis_protocol);

% all the rest
fid = fopen(expdirfile,'w');
for i = 1:nexps,
  fprintf(fid,'%s\n',expdirs{i});
end
fclose(fid);

fprintf('Execute:\n');
fprintf('./fork_FlyBowlAutomaticChecks_Incoming.pl %s\n',expdirfile);

%% register trx

% one experiment
i = 1;
expdir = expdirs{i};
FlyBowlRegisterTrx(expdir,'analysis_protocol',analysis_protocol);

% all the rest
fprintf('Execute:\n');
fprintf('./qsub_FlyBowlRegisterTrx.pl %s\n',expdirfile);


%% sex classification

% one experiment
i = 1;
expdir = expdirs{i};
FlyBowlClassifySex2(expdir,'analysis_protocol',analysis_protocol);

% all the rest
fprintf('Execute:\n');
fprintf('./qsub_FlyBowlClassifySex.pl %s\n',expdirfile);

diagnosticsfilestr = 'sexclassifier_diagnostics.txt';

% plot the results
clear diagnostics;
for i = 1:nexps,
  diagnosticsfile = fullfile(expdirs{i},diagnosticsfilestr);
  diagnostics(i) = ReadParams(diagnosticsfile);
end

fns = fieldnames(diagnostics);
nfns = numel(fns);
nr = floor(sqrt(nfns));
nc = ceil(nfns/nr);
hfig = 1;
figure(hfig);
clf;
hax = createsubplots(nr,nc,.02);
if nfns < numel(hax),
  delete(hax(nfns+1:end));
  hax = hax(1:nfns);
end

clear metadata;
for i = 1:nexps,
  metadata(i) = parseExpDir(expdirs{i});
end
timestamps = {metadata.date};
datenums = datenum(timestamps,'yyyymmddTHHMMSS');
idxcontrol = strcmpi({metadata.line},'pBDPGAL4U');

for i = 1:nfns,
  fn = fns{i};
  x = [diagnostics.(fn)];
  plot(hax(i),datenums(idxcontrol),x(idxcontrol),'k.');
  hold(hax(i),'on');
  plot(hax(i),datenums(~idxcontrol),x(~idxcontrol),'r.');
  datetick(hax(i),'x','mmmyy');
  if mod(i,nc) > 0,
    set(hax(i),'XTickLabel',{});
  end
  axisalmosttight([],hax(i));
  title(hax(i),fn,'Interpreter','none');
end
linkaxes(hax,'x');

%% compute per-frame features

% one experiment
i = 1;
expdir = expdirs{i};
FlyBowlComputePerFrameFeatures(expdir,'analysis_protocol',analysis_protocol);

% all the rest
fprintf('Execute:\n');
fprintf('./qsub_FlyBowlComputePerFrameFeatures.pl %s\n',expdirfile);

%% compare results

i = 1;
expdir = expdirs{i};
inexpdir = inexpdirs{i};
experiment_name = experiment_names{i};
trx = Trx('analysis_protocol',analysis_protocol,'settingsdir',settingsdir);
trx.AddExpDir(expdir);

oldtrx = Trx('analysis_protocol',analysis_protocol,'settingsdir',settingsdir);
oldtrx.AddExpDir(inexpdir);

perframefnsfile = fullfile(trx.settingsdir,trx.analysis_protocol,trx.dataloc_params.perframefnsfilestr);
perframefns = importdata(perframefnsfile);

for fly = 1:trx.nflies,
  if fly > oldtrx.nflies,
    break;
  end
  if trx(fly).firstframe ~= oldtrx(fly).firstframe,
    fprintf('Experiment %s: trajectories %d do not match.\n',experiment_name,fly);
    continue;
  end
  n = trx(fly).nframes;
  oldn = oldtrx(fly).nframes;
  n1 = min(n,oldn)-2;
  if n1 < 1,
    continue;
  end
  for j = 1:nfns,
    fn = perframefns{j};
    if ismember(fn,{'sex'}) || ~isempty(regexp(fn,'^closest','once'))
      continue;
    end
    data = trx(fly).(fn)(1:n1);
    olddata = oldtrx(fly).(fn)(1:n1);
    z = std([data,olddata]);
    if ~isempty(regexp(fn,'dangle','once')) || ~isempty(regexp(fn,'dtheta','once')),
      maxdiff = max(abs(modrange(data-olddata,-pi,pi)))/z;
    else
      maxdiff = max(abs(data-olddata))/z;
    end
    if maxdiff > .01
      fprintf('max diff for %s, %d, %s = %f stds\n',experiment_name,fly,fn,maxdiff);
      keyboard;
    end
  end
end

%% compute per-frame stats

% one experiment
i = 1;
expdir = expdirs{i};
FlyBowlComputePerFrameStats(expdir,'analysis_protocol',analysis_protocol);

% all the rest
fprintf('Execute:\n');
fprintf('./qsub_FlyBowlComputePerFrameStats.pl %s\n',expdirfile);

%% plot per-frame stats

% one experiment
i = 1;
expdir = expdirs{i};
FlyBowlPlotPerFrameStats(expdir,'analysis_protocol',analysis_protocol);

% all the rest
fprintf('Execute:\n');
fprintf('./qsub_FlyBowlPlotPerFrameStats.pl %s\n',expdirfile);

%% compute extra diagnostics

% one experiment
i = 1;
expdir = expdirs{i};
FlyBowlExtraDiagnostics(expdir,'analysis_protocol',analysis_protocol);

% all the rest
fprintf('Execute:\n');
fprintf('./qsub_FlyBowlExtraDiagnostics.pl %s\n',expdirfile);

%% ctrax results movie

% one experiment
i = 1;
expdir = expdirs{i};
FlyBowlMakeCtraxResultsMovie(expdir,'analysis_protocol',analysis_protocol);

% all the rest
fprintf('Execute:\n');
fprintf('./qsub_FlyBowlMakeCtraxResultsMovie.pl %s\n',expdirfile);

%% analysis protocol file for current

for i = 1:numel(expdirs),
  expdir = expdirs{i};
  cmd = sprintf('./analysis_protocol.pl %s %s %s',expdir,analysis_protocol,ctrax_protocol);
  unix(cmd);
end


%% completed checks

% one experiment
i = 1;
expdir = expdirs{i};
FlyBowlAutomaticChecks_Complete(expdir,'analysis_protocol',analysis_protocol);

% all the rest
fprintf('Execute:\n');
fprintf('./qsub_FlyBowlAutomaticChecks_Complete.pl %s\n',expdirfile);

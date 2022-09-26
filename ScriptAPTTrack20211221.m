modpath;
apttrkinfofile = '/groups/branson/bransonlab/projects/FlyBubble/APTTrackers/APTTrackerGrone20220329/APTInfo20220712.mat';

%% set up apt tracker info struct
% only need to do this once

tarfile = '/groups/branson/bransonlab/projects/FlyBubble/APTTrackers/APTTrackerGrone20220329.tar';
lObj.projBundleTempDir(tarfile);

aptinfo = struct;
aptinfo.aptpath ='/groups/branson/home/bransonk/tracking/code/APT';
aptinfo.apttracker.name = '20220329T161748';
aptinfo.apttracker.view = 1;
aptinfo.apttracker.cache = '/groups/branson/bransonlab/projects/FlyBubble/APTTrackers/APTTrackerGrone20220329';
aptinfo.apttracker.projname = 'multitarget_bubble_training_20210523_allGT_AR_params20210920';
aptinfo.apttracker.type = 'mdn_joint_fpn';
aptinfo.apttracker.model = fullfile(aptinfo.apttracker.cache,aptinfo.apttracker.projname,aptinfo.apttracker.type,sprintf('view_%d',aptinfo.apttracker.view-1),aptinfo.apttracker.name,'deepnet-60000');
aptinfo.apttracker.strippedlblfile = fullfile(aptinfo.apttracker.cache,aptinfo.apttracker.projname,'20220329T161748_20220329T161749.lbl');
aptinfo.apttracker.trackconfigfile = fullfile(aptinfo.apttracker.cache,aptinfo.apttracker.projname,'track_config.json');
aptinfo.bindpaths = {'/groups/branson','/nrs/branson'};
aptinfo.singularityimg = '/groups/branson/bransonlab/apt/sif/prod.sif';

lObj.tracker.trkCreateConfig(aptinfo.apttracker.trackconfigfile);

assert(exist(aptinfo.apttracker.cache,'dir')>0);
assert(exist(aptinfo.apttracker.model,'file')>0);
assert(exist(aptinfo.apttracker.strippedlblfile,'file')>0);
assert(exist(aptinfo.apttracker.trackconfigfile,'file')>0);

save(apttrkinfofile,'-struct','aptinfo');

%% data to track

rootdatadir = '/groups/branson/bransonlab/flydisco_data';
totrackfile = '/groups/branson/home/bransonk/behavioranalysis/code/FlyDiscoAnalysis/APT_requestlist_20220712.csv';
totrack = importdata(totrackfile);
expnames = totrack;

% expnames = {...
%   'VNC_JRC_SS44225_RigA_20210928T134724'
%   'VNC_JRC_SS47864_RigA_20210426T132018'
%   'VNC_JRC_SS48413_RigA_20210504T144707'
%   'VNC_JRC_SS48619_RigA_20210428T153022'
%   'VNC_JRC_SS67373_RigD_20211027T165117'
%   'VNC_JRC_SS67404_RigA_20211027T153028'
%   'VNC_JRC_SS70427_RigA_20211027T154243'
%   'VNC_JRC_SS71949_RigA_20211027T155921'
%   'VNC_JRC_SS72073_RigA_20211027T161145'
%   'VNC_JRC_SS73717_RigA_20211027T170019'
%   'VNC_JRC_SS73769_RigA_20211007T150616'
%   'VNC_JRC_SS73782_RigA_20210914T155515'
%   
%   };

analysis_protocol = 'current_non_olympiad_dickson_VNC';
settingsdir = '/groups/branson/home/bransonk/behavioranalysis/code/FlyDiscoAnalysis/settings';
datalocparamsfilestr = 'dataloc_params.txt';

%% get APT tracker

load(apttrkinfofile);

assert(exist(apttracker.cache,'dir')>sshcmd0);
assert(exist(apttracker.model,'file')>0);
assert(exist(apttracker.strippedlblfile,'file')>0);
assert(exist(apttracker.trackconfigfile,'file')>0);

timestamp = datestr(now,'yyyymmddTHHMMSS');
rundirstr = 'run';
mintimestamp = 0;

%% exps to track

expdirs = cellfun(@(x) fullfile(rootdatadir,x),expnames,'Uni',0);
assert(all(cellfun(@exist,expdirs)>0));
nexps = numel(expdirs);

dataloc_params = ReadParams(fullfile(settingsdir,analysis_protocol,datalocparamsfilestr));
% adding to dataloc_params cuz not in there yet


%% APT track

binpathstr = sprintf(' -B %s',bindpaths{:});

for i = 1:nexps,

  expdir = expdirs{i};
  [~,expname] = fileparts(expdir);
  
  rundir = fullfile(expdir,rundirstr);
  if ~exist(rundir,'dir'),
    unix(sprintf('ssh bransonlab@login1 "mkdir ''%s''"',rundir));
  end
  %assert(exist(rundir,'dir')>0);

  errfile = fullfile(rundir,sprintf('apt_%s.err',expname,timestamp));
  apttrkfile = fullfile(expdir,dataloc_params.apttrkfilestr);
  movfile = fullfile(expdir,dataloc_params.moviefilestr);
  trxfile = fullfile(expdir,dataloc_params.trxfilestr);
  
  if exist(apttrkfile,'file'),
    res = dir(apttrkfile);
    if res.datenum > mintimestamp,
      fprintf('%s apt tracked, skipping\n',expname);
      continue;
    end
  end
  
  logfile = fullfile(rundir,sprintf('apt_%s.out',timestamp));
  resfile = fullfile(rundir,sprintf('apt_%s_bsub.out',timestamp));
  reserrfile = fullfile(rundir,sprintf('apt_%s_bsub.err',timestamp));
  
  %aptextra = ' -trx_ids 1 -start_frame 1 -end_frame 201';
  aptextra = '';
  aptcmd = sprintf('python %s/deepnet/APT_interface.py -name %s -view %d -cache %s -err_file %s -model_files %s -type %s %s track -out %s -config_file %s -mov %s -trx %s%s > %s 2>&1',...
    aptpath,apttracker.name,apttracker.view,apttracker.cache,errfile,apttracker.model,apttracker.type,apttracker.strippedlblfile,apttrkfile,apttracker.trackconfigfile,movfile,trxfile,aptextra,logfile);
  singcmd = sprintf('singularity exec --nv %s %s bash -c "%s"',binpathstr,singularityimg,aptcmd);
  bsubcmd = sprintf('bsub -n 1 -J apt_%s -gpu "num=1" -q gpu_rtx -o %s -e %s -R"affinity[core(1)]" "%s"',expname,resfile,reserrfile,strrep(singcmd,'"','\"'));
  
  sshcmd = sprintf('ssh bransonlab@login1 "%s"',strrep(strrep(bsubcmd,'\','\\'),'"','\"'));
  
  fprintf('APT tracking %d %s:\n%s\n',i,expname,sshcmd);
  unix(sshcmd);
  
end

%% check that files were created

for i = 1:nexps,

  expdir = expdirs{i};
  [~,expname] = fileparts(expdir);
  
  apttrkfile = fullfile(expdir,dataloc_params.apttrkfilestr);
  fprintf('%s: %d\n',expname,exist(apttrkfile,'file'));
end

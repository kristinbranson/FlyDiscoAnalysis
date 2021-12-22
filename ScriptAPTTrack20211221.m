modpath;
rootdatadir = '/groups/branson/bransonlab/flydisco_data';

expnames = {...
  'VNC_JRC_SS44225_RigA_20210928T134724'
  'VNC_JRC_SS47864_RigA_20210426T132018'
  'VNC_JRC_SS48413_RigA_20210504T144707'
  'VNC_JRC_SS48619_RigA_20210428T153022'
  'VNC_JRC_SS67373_RigD_20211027T165117'
  'VNC_JRC_SS67404_RigA_20211027T153028'
  'VNC_JRC_SS70427_RigA_20211027T154243'
  'VNC_JRC_SS71949_RigA_20211027T155921'
  'VNC_JRC_SS72073_RigA_20211027T161145'
  'VNC_JRC_SS73717_RigA_20211027T170019'
  'VNC_JRC_SS73769_RigA_20211007T150616'
  'VNC_JRC_SS73782_RigA_20210914T155515'
  
  };

analysis_protocol = 'current_non_olympiad_dickson_VNC';
settingsdir = '/groups/branson/home/bransonk/behavioranalysis/code/FlyDiscoAnalysis/settings';
datalocparamsfilestr = 'dataloc_params.txt';

apttrkinfofile = '/groups/branson/home/bransonk/behavioranalysis/code/MABe2022/APTTrackerInfo20211108.mat';

bindpaths = {'/nrs/branson','/groups/branson'};
singularityimg = '/groups/branson/bransonlab/apt/sif/prod.sif';

timestamp = datestr(now,'yyyymmddTHHMMSS');
rundirstr = 'run';
mintimestamp = 0;

%% exps to track

expdirs = cellfun(@(x) fullfile(rootdatadir,x),expnames,'Uni',0);
assert(all(cellfun(@exist,expdirs)>0));
nexps = numel(expdirs);

dataloc_params = ReadParams(fullfile(settingsdir,analysis_protocol,datalocparamsfilestr));
% adding to dataloc_params cuz not in there yet

%% get APT tracker
load(apttrkinfofile);

assert(exist(apttracker.model,'file')>0);
assert(exist(apttracker.strippedlblfile,'file')>0);

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
  aptcmd = sprintf('python %s/deepnet/APT_interface.py -name %s -view %d -cache %s -err_file %s -model_files %s -type %s %s track -out %s -mov %s -trx %s%s > %s 2>&1',...
    aptpath,apttracker.name,apttracker.view,apttracker.cache,errfile,apttracker.model,apttracker.type,apttracker.strippedlblfile,apttrkfile,movfile,trxfile,aptextra,logfile);
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

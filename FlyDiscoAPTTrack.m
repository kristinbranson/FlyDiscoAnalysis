function [success,msgs,jobid] = FlyDiscoAPTTrack(expdir,varargin)

success = true;
msgs = {};
jobid = nan;

version = '0.1';

[analysis_protocol,settingsdir,datalocparamsfilestr,rootoutdir,sshhost,dowait,dryrun,waitchecktime,startframe,endframe,verbose,dooverwrite] = ...
  myparse(varargin,...
  'analysis_protocol','current_bubble',...
  'settingsdir','/groups/branson/home/robiea/Code_versioned/FlyBubbleAnalysis/settings',...
  'datalocparamsfilestr','dataloc_params.txt',...
  'rootoutdir','',...
  'sshhost','',...
  'dowait',true,...
  'dryrun',false,...
  'waitchecktime',10,...
  'startframe',[],...
  'endframe',[],...
  'verbose',1,...
  'dooverwrite',true);

if verbose >= 1,
  fprintf('Starting FlyDiscoAPTTrack version %s\n',version);
end

%% read in the data locations
datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
dataloc_params = ReadParams(datalocparamsfile);

%% read in apt params

% name of parameters file
aptparamsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.aptparamsfilestr);
if ~exist(aptparamsfile,'file'),
  error('APT params file %s does not exist',aptparamsfile);
end
% read
apt_params = ReadParams(aptparamsfile);

%% construct the command

[~,expname] = fileparts(expdir);
movfile = fullfile(expdir,dataloc_params.moviefilestr);
trxfile = fullfile(expdir,dataloc_params.trxfilestr);
if ~isempty(rootoutdir),
  outexpdir = fullfile(rootoutdir,expname);
  if ~exist(outexpdir,'dir'),
    [s,m] = mkdir(outexpdir);
    if ~s,
      error('Could not create output directory %s: %s',outexpdir,m);
    end
  end
else
  outexpdir = expdir;
end
outtrkfile = fullfile(outexpdir,dataloc_params.apttrkfilestr);
if exist(outtrkfile,'file'),
  if ~dooverwrite,
    success = true;
    msg = sprintf('Output trk file %s already exists.',outtrkfile);
    if verbose >= 1,
      fprintf('%s\n',msg);
    end
    msgs{end+1} = msg;
  end

  newloc = rename_file_timestamp(outtrkfile,sshhost);
  if verbose >= 1,
    fprintf('Renamed existing output file %s to %s\n',outtrkfile,newloc);
  end
end

parttrkfile = [outtrkfile,'.part'];
logfile = fullfile(outexpdir,dataloc_params.aptlogfilestr);
errfile = fullfile(outexpdir,dataloc_params.apterrfilestr);
if apt_params.dobsub,
  bsublogfile = fullfile(outexpdir,dataloc_params.aptbsublogfilestr);
  if exist(bsublogfile,'file'),
    newloc = rename_file_timestamp(bsublogfile,sshhost);
    if verbose >= 1,
      fprintf('Renamed existing bsub log file %s to %s\n',bsublogfile,newloc);
    end
  end
else
  bsublogfile = '';
end

if exist(logfile,'file'),
  newloc = rename_file_timestamp(logfile,sshhost);
  if verbose >= 1,
    fprintf('Renamed existing log file %s to %s\n',logfile,newloc);
  end
end

if exist(errfile,'file'),
  newloc = rename_file_timestamp(errfile,sshhost);
  if verbose >= 1,
    fprintf('Renamed existing error file %s to %s\n',errfile,newloc);
  end
end

if ~isempty(endframe) && isempty(startframe),
  startframe = 1;
end
if ~isempty(startframe) && isempty(endframe),
  td = load(trxfile);
  endframe = max([td.trx.endframe]);
end

aptcmd = sprintf('cd "%s"/deepnet; python APT_track.py -lbl_file "%s" -model_ndx %d -mov "%s" -trx "%s" -out "%s" -log_file "%s" -err_file "%s"',...
  apt_params.aptrepopath,apt_params.lbl_file,apt_params.model_ndx,...
  movfile,trxfile,outtrkfile,...
  logfile,errfile);
if ~isempty(startframe),
  aptcmd = sprintf('%s -start_frame %d -end_frame %d',aptcmd,startframe,endframe);
end

bindpathstr = sprintf(' -B %s',apt_params.bindpaths{:});

singcmd = sprintf('singularity exec --nv %s %s bash -c %s',bindpathstr,apt_params.singularityimg,escape_string_for_bash(aptcmd));
if apt_params.dobsub,
  bsubcmd = sprintf('bsub -n 1 -J apt_%s -gpu "num=1" -q %s -o %s -R"affinity[core(1)]" %s 2>&1',expname,apt_params.gpuqueue,bsublogfile,escape_string_for_bash(singcmd));
  if ~isempty(sshhost),
    sshcmd = sprintf('ssh %s %s',sshhost,escape_string_for_bash(bsubcmd));
    cmd = sshcmd;
  else
    cmd = bsubcmd;
  end
else
  cmd = singcmd;
end

if verbose >= 1,
  fprintf('%s\n',cmd);
end

if ~dryrun,
  stdout = system_with_error_handling(cmd) ;
  stdout = strtrim(stdout);
  if verbose >= 1,
    fprintf('%s\n',stdout);
  end
  starttime = tic;
  lastpartupdate = now;
  if apt_params.dobsub,

    jobid = parse_bsub_jobid(stdout);

    if dowait,

      if verbose >= 1,
        fprintf('Waiting for job %d to complete...\n',jobid);
      end


      while true,

        result = get_single_bsub_job_status(jobid,sshhost);

        elapsedtime = toc(starttime);
        if verbose >= 2,
          fprintf('Checking APT track job status after %f seconds...\n',elapsedtime);
        end

        if result == -1, % errored out
          success = false;
          msg = sprintf('APT track command exited without completion.');
          msgs{end+1} = msg;
          if verbose >= 1,
            fprintf('%s\n',msg);
          end
          return;
        elseif result == 1,
          if exist(outtrkfile,'file'),
            if verbose >= 1,
              fprintf('APT track command completed, output trk file found.\n');
            end
            success = true;
            return;
          else
            msg = sprintf('APT track command finished but could not find output file %s',outtrkfile);
            if verbose >= 1,
              fprintf('%s\n',msg);
            end
            success = false;
            msgs{end+1} = msg;
            return;
          end
        end
        % result == 0
        if elapsedtime > apt_params.maxwaittime,
          msg = sprintf('APT track command did not complete within time limit %f s (elapsed time = %f). ',apt_params.maxwaittime,elapsedtime);
          if verbose >= 1,
            fprintf('%s\n',msg);
          end
          success = false;
          msgs{end+1} = msg;
          return;
        end
        if apt_params.maxupdatetime < apt_params.maxwaittime,
          res = dir(parttrkfile);
          if ~isempty(res) && res.datenum > lastpartupdate,
            if verbose >= 2,
              fprintf('Part trk file updated at time %s.\n',res.date);
            end
            lastpartupdate = res.datenum;
          end
          deltaupdate = (now - lastpartupdate)*24*60*60;
          if verbose >= 2,
            fprintf('Time since last part trk file update = %f s\n',deltaupdate);
          end
          if deltaupdate > apt_params.maxupdatetime,

            msg = sprintf('No update to part track file %s in %f > %f seconds',parttrkfile,deltaupdate,apt_params.maxupdatetime);
            if verbose >= 1,
              fprintf('%s\n',msg);
            end

            success = false;
            msgs{end+1} = msg;
            return;
          end
        end
        if verbose >= 2,
          fprintf('Waiting %f seconds to check again...\n',waitchecktime);
        end
        pause(waitchecktime);
      end
    end
  end
end



%python APT_track.py -lbl_file '/groups/branson/home/bransonlab/APTlabelfiles/20220826_Grone_60K_aligned.lbl' -model_ndx 1 -mov '/groups/branson/bransonlab/flydisco_data/VNC2_EXT_VGLUT-GAL4_RigA_20220419T094935/movie.ufmf' -trx '/groups/branson/bransonlab/flydisco_data/VNC2_EXT_VGLUT-GAL4_RigA_20220419T094935/registered_trx.mat' -start_frame 1 -end_frame 1000 -out $PWD/"VNC2_EXT_VGLUT-GAL4_RigA_20220419T094935.trk"

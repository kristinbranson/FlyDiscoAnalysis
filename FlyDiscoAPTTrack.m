function [success,msgs,jobid] = FlyDiscoAPTTrack(expdir,varargin)

success = true;
msgs = {};
jobid = nan;

version = '0.1';

[analysis_protocol,...
 settingsdir, ...
 dataloc_params, ...
 datalocparamsfilestr, ... 
 rootoutdir, ...
 sshhost, ...
 cluster_billing_account_name, ...
 dowait, ...
 dryrun, ...
 waitchecktime, ...
 startframe, ...
 endframe, ...
 verbose, ...
 dooverwrite, ...
 aptrepopath, ...
 docomputemd5s] = ...
  myparse(varargin,...
  'analysis_protocol','current_bubble',...
  'settingsdir',default_settings_folder_path(),...
  'dataloc_params',[], ...
  'datalocparamsfilestr','dataloc_params.txt',...
  'rootoutdir','',...
  'sshhost','',...
  'cluster_billing_account_name', '', ...
  'dowait',true,...
  'dryrun',false,...
  'waitchecktime',10,...
  'startframe',[],...
  'endframe',[],...
  'verbose',1,...
  'dooverwrite',true, ...
  'aptrepopath','', ...
  'docomputemd5s', false);
% Empty sshhost means to bsub locally

if verbose >= 1,
  fprintf('Starting FlyDiscoAPTTrack version %s\n',version);
end

%% read in the data locations if not passed in via dataloc_params
if isempty(dataloc_params) ,
  datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
  dataloc_params = ReadParams(datalocparamsfile);
end

%% read in apt params

% name of parameters file
aptparamsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.aptparamsfilestr);
if ~exist(aptparamsfile,'file'),
  error('APT params file %s does not exist',aptparamsfile);
end
% read
apt_params = ReadParams(aptparamsfile);

%% if aptrepopath was passed as an argument, override the one specified in the APT params file
if ~isempty(aptrepopath),
  apt_params.aptrepopath = aptrepopath ;
end
  
%% Print the md5 hashes of the label file and the singularity image
fprintf('Label file is %s\n', apt_params.lbl_file) ;
[~,~,~,lbl_file_modification_datetime] = simple_dir(apt_params.lbl_file) ;
lbl_file_modification_datetime.TimeZone = 'local' ;
lbl_file_modification_datetime.Format = 'yyyy-MM-dd HH:mm:ss z' ;
fprintf('The modification time of the label file is %s\n', string(lbl_file_modification_datetime)) ;
if docomputemd5s ,
  label_file_md5 = compute_md5_on_local(apt_params.lbl_file) ;
  fprintf('The md5 hash of the label file is %s\n', label_file_md5) ;
end
fprintf('Singularity image file is %s\n', apt_params.singularityimg) ;
[~,~,~,sing_image_file_modification_datetime] = simple_dir(apt_params.singularityimg) ;
sing_image_file_modification_datetime.TimeZone = 'local' ;
sing_image_file_modification_datetime.Format = 'yyyy-MM-dd HH:mm:ss z' ;
fprintf('The modification time of the singularity image file is %s\n', string(sing_image_file_modification_datetime)) ;
if docomputemd5s ,
  singularity_image_file_md5 = compute_md5_on_local(apt_params.singularityimg) ;
  fprintf('The md5 hash of the singularity image file is %s\n', singularity_image_file_md5) ;
end

% %% Print a bunch of of useful info about the .lbl file, APT repo status
% This has issues b/c APT and JAABA have namespace collisions
% APT_Labeler_printInfo_without_adding_APT_to_path(apt_params.lbl_file) ;

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
  P_option = fif(isempty(cluster_billing_account_name), '', sprintf('-P %s', cluster_billing_account_name)) ;  
  bsubcmd = sprintf('bsub %s -n 1 -J apt_%s -gpu "num=1" -q %s -o %s -R"affinity[core(1)]" %s 2>&1', ...
                    P_option, ...
                    expname, ...
                    apt_params.gpuqueue, ...
                    bsublogfile, ...
                    escape_string_for_bash(singcmd));
  if ~isempty(sshhost),
    sshcmd = sprintf('ssh -o BatchMode=yes -o StrictHostKeyChecking=no -o ConnectTimeout=20 %s %s', sshhost, escape_string_for_bash(bsubcmd));
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

        numeric_bsub_job_status = get_bsub_job_status(jobid,sshhost);

        elapsedtime = toc(starttime);
        if verbose >= 2,
          fprintf('Checking APT track job status after %f seconds...\n',elapsedtime);
        end

        if numeric_bsub_job_status == -1, % errored out
          success = false;
          msg = sprintf('APT track command exited without completion.');
          msgs{end+1} = msg; %#ok<AGROW>
          if verbose >= 1,
            fprintf('%s\n',msg);
          end
          return;
        elseif numeric_bsub_job_status == 1,
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
            msgs{end+1} = msg; %#ok<AGROW>
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
          msgs{end+1} = msg; %#ok<AGROW>
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
            msgs{end+1} = msg; %#ok<AGROW>
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

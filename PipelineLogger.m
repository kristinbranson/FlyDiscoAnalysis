classdef PipelineLogger < handle
  
  properties   
    expdir = [];
    stage;
    
    fid = [];
    runInfo;
  end
    
  methods 
    
    function obj = PipelineLogger(expdir,pipelinestage,...
        dataloc_params,logfile_datalocfield,...
        settingsdir,analysis_protocol,...
        varargin)
      % Construct/initialize Logger and print header
      % expdir: char, fullpath
      % pipelinestage: char, eg 'FlyBubbleClassifySex'
      % 
      % Optional PVs:
      % logfid. If provided, Logger takes ownership of this (should be live) file handle.
      % debug. logical scalar.
      % versionstr. Caller version str for print
      
      [logfid,debug,versionstr] = myparse(varargin,'logfid',[],'debug',false,'versionstr','');
      
      obj.expdir = expdir;
      obj.stage = pipelinestage;
      
      if isempty(logfid)
        if isfield(dataloc_params,logfile_datalocfield) && ~debug
          logfile = fullfile(expdir,dataloc_params.(logfile_datalocfield));
          logfid = fopen(logfile,'a');
          if logfid < 1,
            warning('PipelineLogger:log','Could not open log file %s\n',logfile);
            logfid = 1;
          end
        else
          logfid = 1;
        end
      end
      obj.fid = logfid;
      
      if isdeployed
        fprintf(logfid,'Running deployed. CTFRoot=%s\n',ctfroot);
      end
      
      real_analysis_protocol = GetRealAnalysisProtocol(analysis_protocol,settingsdir);
      settingsSS = FlyBubbleBaR.settingssnapshot(settingsdir);
      codeSS = FlyBubbleBaR.codesnapshot();
      timestamp = datestr(now,'yyyymmddTHHMMSS');
      if isempty(versionstr)
        fprintf(logfid,'\n\n***\nRunning %s analysis_protocol %s (real analysis protocol %s) at %s\n',...
          obj.stage,analysis_protocol,real_analysis_protocol,timestamp);
      else
        fprintf(logfid,'\n\n***\nRunning %s ver %s analysis_protocol %s (real analysis protocol %s) at %s\n',...
          obj.stage,versionstr,analysis_protocol,real_analysis_protocol,timestamp);
      end
      fprintf(logfid,'\nCode snapshot: \n');
      fprintf(logfid,'%s\n',codeSS{:});
      fprintf(logfid,'\nSettings snapshot: \n');
      fprintf(logfid,'%s\n',settingsSS{:});
      
      info = struct();
      info.version = versionstr;
      info.analysis_protocol = analysis_protocol;
      info.linked_analysis_protocol = real_analysis_protocol;
      info.settings_snapshot = settingsSS;
      info.code_snapshot = codeSS;
      info.timestamp = timestamp;
      obj.runInfo = info;
    end
    
    function delete(obj)
      obj.closeAbrupt();
    end
    
    function log(obj,varargin)
      fprintf(obj.fid,varargin{:});
    end
    
    function close(obj)
      obj.log('Finished running %s at %s.\n',obj.stage,datestr(now,'yyyymmddTHHMMSS'));
      obj.closeAbrupt();
    end
    
    function closeAbrupt(obj)
      if ~isempty(obj.fid) && obj.fid>1
        fclose(obj.fid);
        obj.fid = [];
      end    
    end
    
  end
  
end
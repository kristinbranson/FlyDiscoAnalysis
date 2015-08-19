classdef FlyBubbleBaR
  % FlyBubble Build and Run
    
  properties (Constant)
    
    Root = fileparts(mfilename('fullpath'));
    
    MANIFESTFILE = 'FlyBubbleManifest.txt';
    
    Manifest = lclReadManifest( fullfile(FlyBubbleBaR.Root,FlyBubbleBaR.MANIFESTFILE) );
    
    SnapshotScript = fullfile(FlyBubbleBaR.Root,'repo_snapshot.sh');    

    BUILDPROJECTS = {
      'FlyBubbleAutomaticChecks_Incoming'
      'FlyBubbleClassifySex'
      'FlyBubbleComputePerFrameFeatures'
      'FlyBubbleDectectIndicatorLedOnOff'
      'FlyBubbleRegisterTrx'
      'FlyBubbleTrackWings'
      };
      
    BUILDSNAPSHOTFILE = 'build.snapshot';    
    BUILDSNAPSHOTFULLFILE = fullfile(FlyBubbleBaR.Root,FlyBubbleBaR.BUILDSNAPSHOTFILE);
    
    BUILDMCCFILE = 'build.mcc';
    BUILDMCCFULLFILE = fullfile(FlyBubbleBaR.Root,FlyBubbleBaR.BUILDMCCFILE);
    
  end
  
  methods (Static)
    
    function p = getpath()
      m = FlyBubbleBaR.Manifest;      
      
      hmmroot = m.hmm;
      jctroot = m.jctrax;
      
      % AL: Order important here. FBA has code that shadows eg hmm
      p = { ...
        FlyBubbleBaR.Root; ...
        fullfile(jctroot,'filehandling'); ...
        fullfile(jctroot,'misc'); ...
        fullfile(jctroot,'simplewing'); ...
        hmmroot;
        };
    end
    
    function setpath()
      p = FlyBubbleBaR.getpath();
      addpath(p{:},'-begin');      
    end
    
    function buildAll()
      projs = FlyBubbleBaR.BUILDPROJECTS;
      for p=projs(:)',p=p{1}; %#ok<FXSET>
        fprintf(' ######\n');
        fprintf(' BUILDING %s...\n',p);
        pause(2);
        FlyBubbleBaR.build(p);
      end
    end
    
    function build(proj)
      % build(proj)
      %
      % FlyBubbleBaR.build FlyBubbleTrackWings
            
      % take snapshot + save it to snapshot file
      codeSSfname = FlyBubbleBaR.BUILDSNAPSHOTFULLFILE;
      fprintf('Taking code snapshot and writing to file: %s...\n',codeSSfname);
      codeSS = FlyBubbleBaR.codesnapshot();
      codeSS = cellstr(codeSS);
      cellstrexport(codeSS,codeSSfname);
      fprintf('... done with snapshot.\n');

      pth = FlyBubbleBaR.getpath();
      pth = pth(:);
      Ipth = [repmat({'-I'},numel(pth),1) pth];
      Ipth = Ipth';
      
      fbroot = FlyBubbleBaR.Root;
      
      % include all compute*.m. For now only used by computePFF but I don't
      % see how it hurts
      cpffs = dir(fullfile(fbroot,'compute*.m'));
      cpffs = {cpffs.name}';
      cpffs = cellfun(@(x)fullfile(fbroot,x),cpffs,'uni',0);
      dashACPFFs = [repmat({'-a'},size(cpffs)) cpffs];
      dashACPFFs = dashACPFFs';
      
      mccargs = {...
       '-o',proj,...
       '-W',['main:' proj],...
       '-T','link:exe',...
       '-d',fbroot,...
       '-w','enable:specified_file_mismatch',...
       '-w','enable:repeated_file',...
       '-w','enable:switch_ignored',...
       '-w','enable:missing_lib_sentinel',...
       '-w','enable:demo_license',...
       '-R','-singleCompThread',...
       '-v',fullfile(fbroot,[proj '.m']),...
       Ipth{:},...
       '-a',fullfile(fbroot,'build.snapshot'),...
       '-a',fullfile(fbroot,'repo_snapshot.sh'),...
       '-a',fullfile(fbroot,'FlyBubbleManifest.txt'),...
       dashACPFFs{:}}; %#ok<CCAT>
     
      fprintf('Writing mcc args to file: %s...\n',FlyBubbleBaR.BUILDMCCFULLFILE);
      cellstrexport(mccargs,FlyBubbleBaR.BUILDMCCFULLFILE);
      
      today = datestr(now,'yyyymmdd');
      fprintf('BEGIN BUILD on %s\n',today);
      pause(2.0);
      mcc(mccargs{:});
      
      % postbuild
      bindir = fullfile(fbroot,'builds','bubble',today);      
      if exist(bindir,'dir')==0
        fprintf('Creating bin dir %s...\n',bindir);
        [succ,msg] = mkdir(bindir);
        if ~succ
          error('FlyBubbleBaR:build','Failed to create bin dir: %s\n',msg);
        end
      end
      fprintf('Moving binaries + build artifacts into: %s\n',bindir);
      % move buildmcc file, buildsnapshot file into bindir with name change
      % move binaries
      binsrc = fullfile(fbroot,proj);
      bindst = fullfile(bindir,proj);
      runsrc = fullfile(fbroot,['run_' proj '.sh']);
      rundst = fullfile(bindir,['run_' proj '.sh']);
      mccsrc = FlyBubbleBaR.BUILDMCCFULLFILE;
      mccdst = fullfile(bindir,[proj '.' FlyBubbleBaR.BUILDMCCFILE]);
      sssrc = FlyBubbleBaR.BUILDSNAPSHOTFULLFILE;
      ssdst = fullfile(bindir,[proj '.' FlyBubbleBaR.BUILDSNAPSHOTFILE]);
      
      if exist(bindst,'file')>0 || exist(rundst,'file')>0 || ...
         exist(mccdst,'file')>0 || exist(ssdst,'file')>0
        warning('FlyBubbleBaR:build','Overwriting existing files in bin dir.');
      end
      FlyBubbleBaR.buildmv(binsrc,bindst);
      FlyBubbleBaR.buildmv(runsrc,rundst);
      FlyBubbleBaR.buildmv(mccsrc,mccdst);
      FlyBubbleBaR.buildmv(sssrc,ssdst);
      fileattrib(bindst,'+x');
      fileattrib(rundst,'+x');
      
      mccExc = fullfile(fbroot,'mccExcludedFiles.log');
      readme = fullfile(fbroot,'readme.txt');
      if exist(mccExc,'file')>0
        delete(mccExc);
      end
      if exist(readme,'file')>0
        delete(readme);
      end
    end
    
    function s = codesnapshot
      % Return status of FlyBubble code and dependencies.
      % When deployed, this returns the state at the time of build.

      if isdeployed
        s = importdata('build.snapshot');        
      elseif isunix
        
        % This method assumes that the user has set their path using
        % FlyBubbleBaR.setpath (so that the Manifest correclty reflects
        % dependencies). Do a quick+dirty check of this assumption.
        twh = which('TrackWingsHelper');
        manifest = FlyBubbleBaR.Manifest;
        if ~isequal(fileparts(twh),fullfile(manifest.jctrax,'simplewing'))
          warning('FlyBubbleBaR:manifest',...
            'Runtime path appears to differ from that specified by Manifest.txt. Code snapshot is likely to be incorrect.');
        end      
        
        script = FlyBubbleBaR.SnapshotScript;
        cmd = sprintf('%s -nocolor -brief %s',script,FlyBubbleBaR.Root);
        [~,s] = system(cmd);
        s = regexp(s,sprintf('\n'),'split');

        cmd = sprintf('%s -nocolor -brief %s',script,manifest.jctrax);
        [~,s2] = system(cmd);
        s2 = regexp(s2,sprintf('\n'),'split');
        
        s = [s(:);s2(:)];
      else
        s = {'No snapshot available, non-unix interactive session.'};
      end
    end
    
    function s = settingssnapshot(settingsdir)
      assert(exist(settingsdir,'dir')>0,'Cannot find dir ''%s''.',settingsdir);
      script = FlyBubbleBaR.SnapshotScript;
      cmd = sprintf('%s -nocolor -brief %s',script,settingsdir);
      [~,s] = system(cmd);
      s = regexp(s,sprintf('\n'),'split');
      s = s(:);
    end
    
    function buildmv(src,dst)
      [succ,msg] = movefile(src,dst);
      if ~succ
        error('FlyBubbleBaR:build','Failed to move file ''%s'' -> ''%s'': %s\n',...
          src,dst,msg);
      end
    end
    
  end  
  
end

function s = lclReadManifest(fname)
tmp = importdata(fname);
tmp = regexp(tmp,',','split');
tmp = cat(1,tmp{:});
s = cell2struct(tmp(:,2),tmp(:,1));
end
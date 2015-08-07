classdef FlyBubbleBaR
  % FlyBubble Build and Run
    
  properties (Constant)
    
    Root = fileparts(mfilename('fullpath'));
    
    MANIFESTFILE = 'FlyBubbleManifest.txt';
    
    Manifest = lclReadManifest( fullfile(FlyBubbleBaR.Root,FlyBubbleBaR.MANIFESTFILE) );
    
    SnapshotScript = fullfile(FlyBubbleBaR.Root,'repo_snapshot.sh');
    
    BUILDSNAPSHOTFILE = 'build.snapshot';    
    BUILDSNAPSHOTFULLFILE = fullfile(FlyBubbleBaR.Root,FlyBubbleBaR.BUILDSNAPSHOTFILE);
    
    BUILDMCCFILE = 'build.mcc';
    BUILDMCCFULLFILE = fullfile(FlyBubbleBaR.Root,FlyBubbleBaR.BUILDMCCFILE);
    
  end
  
  methods (Static)
    
    function p = getpath()
      m = FlyBubbleBaR.Manifest;      
      jctroot = m.jctrax;
      p = { ...
        FlyBubbleBaR.Root; ...
        fullfile(jctroot,'filehandling'); ...
        fullfile(jctroot,'misc'); ...
        fullfile(jctroot,'simplewing')};
    end
    
    function setpath()
      m = FlyBubbleBaR.Manifest;
      
      jctroot = m.jctrax;
      addpath(fullfile(jctroot,'filehandling'));
      addpath(fullfile(jctroot,'misc'));
      addpath(fullfile(jctroot,'simplewing'));
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
       '-a',fullfile(fbroot,'FlyBubbleManifest.txt')}; %#ok<CCAT>
     
      fprintf('Writing mcc args to file: %s...\n',FlyBubbleBaR.BUILDMCCFULLFILE);
      cellstrexport(mccargs,FlyBubbleBaR.BUILDMCCFULLFILE);
      
      fprintf('BEGIN BUILD\n');
      pause(2.0);
      mcc(mccargs{:});
      
      % post
      % move build outputs, snapshot, diary into build dirs
      % change permissions for binaries
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
        
        cmd = sprintf('%s -nocolor -brief %s',script,manifest.jctrax);
        [~,tmp] = system(cmd);
        
        s = [s sprintf('\n') tmp];
        s = cellstr(s);
      else
        s = {'No snapshot available, non-unix interactive session.'};
      end
    end
    
    function s = settingssnapshot(settingsdir)
      assert(exist(settingsdir,'dir')>0,'Cannot find dir ''%s''.',settingsdir);
      script = FlyBubbleBaR.SnapshotScript;
      cmd = sprintf('%s -nocolor -brief %s',script,settingsdir);
      [~,s] = system(cmd);
      s = cellstr(s);
    end
    
  end  
end

function s = lclReadManifest(fname)
tmp = importdata(fname);
tmp = regexp(tmp,',','split');
tmp = tmp{1};
s = cell2struct(tmp(:,2),tmp(:,1));
end
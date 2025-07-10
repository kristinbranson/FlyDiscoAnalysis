function JAABADetectWrapper(expdir, settingsdir, analysis_protocol, forcecompute)
  datalocparamsfilestr = 'dataloc_params.txt';
  datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
  dataloc_params = ReadParams(datalocparamsfile);
  jaabaclassifierparamsfilestrs = fullfile(settingsdir,analysis_protocol,dataloc_params.jaabaclassifierparamsfilestrs);
  jabfiles = read_one_file_name_per_line(jaabaclassifierparamsfilestrs) ;
  fprintf('Running JAABADetect...\n');
  
%   % JAABADetect needs to add things to the path, and those things mess up the rest
%   % of the FlyDiscoAnalysis code, so we'll restore the path when JAABADetect() is
%   % done
%   saved_path = path() ;
%   cleaner = onCleanup(@()(path(saved_path))) ;  % restore the path on exit, whether normally or via exception
% 
%   % Actually call JAABADetect()
%   % For reasons that are unclear to me, the otherwise-useless try/catch
%   % wrapper ensures that errors encounted during loading of user-defined
%   % class objects are still silently ignored, even when "dbstop if error"
%   % is engaged.  ALT, 2021-03-04
%   try
%     JAABADetect(expdir, 'jabfiles', jabfiles, 'forcecompute', forcecompute) ;
%   catch me ,
%     rethrow(me) ;
%   end
  
  % JAABADetect needs its own path, and also seems to (maybe?) hold onto a lot of
  % memory, so we run it in a subshell.
  do_actually_shell_out = false ;
  stdout = run_in_sub_matlab(do_actually_shell_out, @JAABADetect, expdir, 'jabfiles', jabfiles, 'forcecompute', forcecompute) ;
  fprintf('%s', stdout) ; 
end  

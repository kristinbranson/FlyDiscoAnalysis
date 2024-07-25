function FlyDiscoPipelineCore(expdir, ...
                              stage_name_from_stage_index, ...
                              do_run, ...
                              required_file_names_from_stage_name, ...
                              debug, ...
                              settingsdir, ...
                              analysis_protocol, ...
                              stage_function_from_stage_name, ...
                              stage_additional_arguments_from_stage_name, ...
                              dataloc_params)

% Fundamentally, the job of FlyDiscoPipelineCore() is to generate all the
% appropriate analysis files from the `raw' data files. If it generates all the
% appropriate analysis files, that's a successful run. So if the auto-checks
% incoming (ACI) code finds that there are dead flies, then the appropriate
% files to generate for that experiment are just
% automatic_checks_incoming_results.txt, automatic_checks_incoming_info.mat,
% automatic_checks_complete_results.txt, and automatic_checks_complete_info.mat.
% If that is done successfully, then FlyDiscoPipelineCore() should return normally,
% and not throw an error. Thus the throwing of uncaught errors should be
% restricted to cases where it's not even clear to the programmer how to proceed
% to generate the output files for a particular stage. The simplest case would
% be if a write to disk of one of the output files fails, that should clearly
% lead to FlyDiscoPipelineCore() throwing an error. Another simple case is if the
% experiment metadata indicates dead flies. That should stop the core pipeline
% stages from running, but should not cause FlyDiscoPipeline() to throw an
% error.

% Loop over all the stages, calling FlyDiscoPipelineStage() for each.  This
% will cause stage_function_from_stage_name.(stage_name) to be called for each
% stage in turn.
do_run_core_stages = true ;
stage_count = numel(stage_name_from_stage_index) ;
for stage_index = 1 : stage_count ,
  stage_name = stage_name_from_stage_index{stage_index} ;
  is_core = ~(strcmp(stage_name, 'automaticchecksincoming') || strcmp(stage_name, 'automaticcheckscomplete')) ;
  % We run the stage if it is non-core (ACI or ACC), OR if do_run_core_stages is
  % true.
  do_run_stage =  ~is_core || do_run_core_stages ;
  if do_run_stage ,
    stage_function = stage_function_from_stage_name.(stage_name) ;
    FlyDiscoPipelineStage(expdir, ...
                          stage_name, ...
                          do_run, ...
                          required_file_names_from_stage_name, ...
                          debug, ...
                          settingsdir, ...
                          analysis_protocol, ...
                          stage_function, ...
                          stage_additional_arguments_from_stage_name) ;
  
    % Special check after ACI stage.  This will disable running of the core stages
    % in some cases
    if strcmp(stage_name, 'automaticchecksincoming') ,
      do_run_core_stages = CheckACIFile(expdir, dataloc_params, do_run) ;
    end  
  end  % if do_run_stage
end  % for loop

end  % function





function do_run_core_stages = CheckACIFile(expdir, dataloc_params, do_run)

if is_on_or_force(do_run.automaticchecksincoming) ,
  % Make sure the file usually named automatic_checks_incoming_results.txt
  % contains either "automated_pf,P" or "automated_pf,U", but not
  % "automated_pf,F" (nor anything else).
  % This ensures the the core pipeline only runs if the current run or a previous run of the ACI stage
  % indicated that it should do so.
  do_run_core_stages = CheckACIResultsFileContents(expdir, dataloc_params) ;
  if ~do_run_core_stages ,
    fprintf('FlyDisco pipeline encountered issues with experiment during automaticchecksincoming stage.  Skipping core pipeline stages.\n');
  end
else
  do_run_core_stages = true ;
end

end  % function
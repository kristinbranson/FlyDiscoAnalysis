function FlyDiscoPipelineCore(expdir, ...
                              stage_name_from_stage_index, ...
                              do_run, ...
                              required_file_names_from_stage_name, ...
                              debug, ...
                              settingsdir, ...
                              analysis_protocol, ...
                              stage_function_from_stage_name, ...
                              stage_additional_arguments_from_stage_name, ...
                              dataloc_params, ...
                              dotrystages)

% Fundamentally, the job of FlyDiscoPipelineCore() is to generate all the
% appropriate analysis files from the `raw' data files. If it generates all
% the appropriate analysis files, that's a successful run. So if the
% auto-checks incoming (ACI) code finds that there are dead flies, then the
% appropriate files to generate for that experiment are just
% automatic_checks_incoming_results.txt, automatic_checks_incoming_info.mat,
% automatic_checks_complete_results.txt, and
% automatic_checks_complete_info.mat. If that is done successfully, then
% FlyDiscoPipelineCore() should return normally, and not throw an error. Thus
% the throwing of uncaught errors should be restricted to cases where one of
% the individual stage functions throws an error.  The simplest case would be
% if a write to disk of one of the output files fails, that should clearly
% lead to FlyDiscoPipelineCore() throwing an error. Another simple case (as
% described above) is if the experiment metadata indicates dead flies. That
% should stop the core pipeline stages from running, but should not cause
% FlyDiscoPipelineCore() to throw an error.

% Deal with optional args
if ~exist('dotrystages', 'var') || isempty(dotrystages) ,
  dotrystages = true ;
end

% Loop over all the stages, calling FlyDiscoPipelineStage() for each.  This
% will cause stage_function_from_stage_name.(stage_name) to be called for each
% stage in turn.
do_run_core_stages = true ;
core_pipeline_exception_maybe = [] ;
exception_stage_name_maybe = {} ;
stage_count = numel(stage_name_from_stage_index) ;
for stage_index = 1 : stage_count ,
  stage_name = stage_name_from_stage_index{stage_index} ;
  is_core = ~(strcmp(stage_name, 'automaticchecksincoming') || strcmp(stage_name, 'automaticcheckscomplete') || strcmp(stage_name, 'cleanup')) ;
  % Decide whether or not to run the stage.
  % Note that logic inside FlyDiscoPipelineStage() can decide not to *actually*
  % run the stage, even if do_run_stage is true.
  do_run_stage =  ~is_core || (do_run_core_stages && isempty(core_pipeline_exception_maybe));
  if do_run_stage ,
    stage_function = stage_function_from_stage_name.(stage_name) ;
    if dotrystages ,
      try
        FlyDiscoPipelineStage(expdir, ...
                              stage_name, ...
                              do_run, ...
                              required_file_names_from_stage_name, ...
                              debug, ...
                              settingsdir, ...
                              analysis_protocol, ...
                              stage_function, ...
                              stage_additional_arguments_from_stage_name) ;
      catch pipeline_exception
        core_pipeline_exception_maybe = pipeline_exception ;
        exception_stage_name_maybe = { stage_name } ;
        fprintf('\nAn exception has occured during stage %s of the pipeline:\n', stage_name) ;
        fprintf('%s\n', pipeline_exception.getReport())
        fprintf('\nThis exception has been caught, but will be rethrown after completion of the ACC phase.\n') ;
      end
    else
      FlyDiscoPipelineStage(expdir, ...
                            stage_name, ...
                            do_run, ...
                            required_file_names_from_stage_name, ...
                            debug, ...
                            settingsdir, ...
                            analysis_protocol, ...
                            stage_function, ...
                            stage_additional_arguments_from_stage_name) ;
    end
  
    % Special check after ACI stage.  This will disable running of the core stages
    % in some cases
    if strcmp(stage_name, 'automaticchecksincoming') ,
      do_run_core_stages = CheckACIFile(expdir, dataloc_params, do_run) ;
    end  
  end  % if do_run_stage
end  % for loop

% If there was an exception thrown during the pipeline, rethrow it so that
% we exit with an error return code when running in batch mode.
if ~isempty(core_pipeline_exception_maybe) ,
  exception_stage_name = exception_stage_name_maybe{1} ; 
  fprintf('\nRethrowing an error that occurred during stage %s of the pipeline:\n', exception_stage_name) ;
  core_pipeline_exception = core_pipeline_exception_maybe(1) ;
  rethrow(core_pipeline_exception) ;
end

% If we get here, analysis has completed successfully
fprintf('Analysis pipeline completed!\n');

end  % function

function do_run_cleanup_stage = CheckACCFile(expdir, dataloc_params, do_run)
% Check if the ACC stage is on, the results file is present and indicates the pipeline passed.

if is_on_or_force(do_run.automaticcheckscomplete) ,
  % Make sure the file usually named automatic_checks_complete_results.txt
  % contains "automated_pf,P", but not
  % "automated_pf,F" (nor anything else).
  % This ensures the the cleanup stage only runs if the current run or a previous run of the ACC stage
  % indicated that it should do so.
  do_run_cleanup_stage = CheckACCResultsFileContents(expdir, dataloc_params) ;
  if ~do_run_cleanup_stage ,
    fprintf('FlyDisco pipeline encountered issues with experiment during automaticcheckscomplete stage.  Skipping cleanup stage.\n');
  end
else
  do_run_cleanup_stage = false ;
end

end  % function

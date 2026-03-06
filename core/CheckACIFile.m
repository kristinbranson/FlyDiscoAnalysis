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

function result = determine_downstream_trx_file_path(expdir, dataloc_params, do_run)
% Determines the path to the trx file for stages downstream of addpflies

if is_on_or_force(do_run.addpflies) ,
  result = fullfile(expdir, determine_addpflies_trx_output_file_name(dataloc_params)) ;
else
  result = fullfile(expdir,dataloc_params.ctraxfilestr) ;
end

end

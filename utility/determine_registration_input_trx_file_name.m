function result = determine_registration_input_trx_file_name(expdir, dataloc_params)

% Determine the relative path of the trx file the registraion stage should use
% as trx input.  Uses addpflies trx output if it exists, otherwise used
% FlyTracker trx output.

ft_trx_output_file_name = dataloc_params.ctraxfilestr ;
ft_trx_output_file_path = fullfile(expdir, ft_trx_output_file_name) ;
addpflies_trx_output_file_name = determine_addpflies_trx_output_file_name(dataloc_params) ;
addpflies_trx_output_file_path = fullfile(expdir, addpflies_trx_output_file_name) ;
if exist(addpflies_trx_output_file_path, 'file') ,
  result = addpflies_trx_output_file_name ;
elseif exist(ft_trx_output_file_path, 'file') ,
  result = ft_trx_output_file_name ;
else
  error('Unable to find any of the possible registration trx input files') ;
end

end



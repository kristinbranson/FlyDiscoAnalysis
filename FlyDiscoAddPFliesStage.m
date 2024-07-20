function FlyDiscoAddPFliesStage(expdir, dataloc_params, settingsdir, analysis_protocol, doaddpflies, debug)

stage = 'addpflies';
addpflies_trx_output_file_name = determine_addpflies_trx_output_file_name(dataloc_params) ;
addpflies_trx_output_file_path = fullfile(expdir, addpflies_trx_output_file_name) ;
if isequal(doaddpflies, 'force') ,
    do_run_stage = true ;
    FlyDiscoAddPFlies(expdir,...
                      'settingsdir', settingsdir, ...
                      'analysis_protocol', analysis_protocol, ...
                      'dataloc_params', dataloc_params, ...
                      'debug', debug) ;
elseif isequal(doaddpflies, 'on') ,
    does_trx_output_file_exist = logical(exist(addpflies_trx_output_file_path, 'file')) ;
    do_run_stage = ~does_trx_output_file_exist ;
    if do_run_stage ,
        FlyDiscoAddPFlies(expdir,...
                          'settingsdir', settingsdir, ...
                          'analysis_protocol', analysis_protocol, ...
                          'dataloc_params', dataloc_params, ...
                          'debug', debug) ;
    end
else
    do_run_stage = false ;
end
if do_run_stage ,            
    does_trx_output_file_exist = logical(exist(addpflies_trx_output_file_path, 'file')) ;
    if ~does_trx_output_file_exist ,
        msgs = { sprintf('Missing file %s', addpflies_trx_output_file_name) } ;
        flydisco_pipeline_error(stage, msgs) ;
    end
end

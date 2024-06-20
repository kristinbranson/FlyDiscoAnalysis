function FlyDiscoAddPFliesStage(expdir, dataloc_params, settingsdir, analysis_protocol, doaddpflies, debug)

stage = 'addpflies';
unregistered_trx_file_path = fullfile(expdir, dataloc_params.ctraxfilestr) ;
if isequal(doaddpflies, 'force') ,
    do_run_stage = true ;
    FlyDiscoAddPFlies(expdir,...
                      'settingsdir', settingsdir, ...
                      'analysis_protocol', analysis_protocol, ...
                      'dataloc_params', dataloc_params, ...
                      'debug', debug) ;
elseif isequal(doaddpflies, 'on') ,
    does_trx_file_have_pflies = determine_does_trx_file_have_pflies(unregistered_trx_file_path) ;
    do_run_stage = ~does_trx_file_have_pflies ;
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
    does_trx_file_have_pflies = determine_does_trx_file_have_pflies(unregistered_trx_file_path) ;
    if ~does_trx_file_have_pflies ,
        msgs = { sprintf('Missing pflies in file %s',unregistered_trx_file_path) } ;
        flydisco_pipeline_error(stage, msgs) ;
    end
end

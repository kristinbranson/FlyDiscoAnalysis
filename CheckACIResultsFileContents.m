function CheckACIResultsFileContents(expdir, dataloc_params, stage)
    aci_results_as_text_file_name = dataloc_params.automaticchecksincomingresultsfilestr ;
    aci_results_as_text_file_path = fullfile(expdir, aci_results_as_text_file_name) ;
    aci_results = ReadParams(aci_results_as_text_file_path) ;
    if isfield(aci_results, 'automated_pf') ,
        if strcmp(aci_results.automated_pf, 'P')  || strcmp(aci_results.automated_pf, 'U') ,
            % do nothing, all is well
        elseif strcmp(aci_results.automated_pf, 'F') ,
            msgs = { sprintf('ACI results file %s has automated_pf==F', aci_results_as_text_file_name) } ;
            flydisco_pipeline_error(stage, msgs) ;
        else
            msgs = { sprintf('ACI results file %s has automated_pf value that is not P nor U', aci_results_as_text_file_name) } ;
            flydisco_pipeline_error(stage, msgs) ;
        end
    else
        msgs = { sprintf('ACI results file %s has no automated_pf field', aci_results_as_text_file_name) } ;
        flydisco_pipeline_error(stage, msgs) ;
    end
end

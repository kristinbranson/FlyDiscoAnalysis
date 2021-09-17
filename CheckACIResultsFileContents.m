function do_continue_pipeline = CheckACIResultsFileContents(expdir, dataloc_params, stage)  %#ok<INUSD>
    aci_results_as_text_file_name = dataloc_params.automaticchecksincomingresultsfilestr ;
    aci_results_as_text_file_path = fullfile(expdir, aci_results_as_text_file_name) ;
    aci_results = ReadParams(aci_results_as_text_file_path) ;
    if isfield(aci_results, 'automated_pf') ,
        if strcmp(aci_results.automated_pf, 'P')  || strcmp(aci_results.automated_pf, 'U') ,
            do_continue_pipeline = true ;
        elseif strcmp(aci_results.automated_pf, 'F') ,
            fprintf('ACI results file %s has automated_pf==F\n', aci_results_as_text_file_name) ;
            do_continue_pipeline = false ;
        else
            fprintf('ACI results file %s has automated_pf value that is not P nor U\n', aci_results_as_text_file_name) ;
            do_continue_pipeline = false ;
        end
    else
        fprintf('ACI results file %s has no automated_pf field\n', aci_results_as_text_file_name) ;
        do_continue_pipeline = false ;
    end
end
